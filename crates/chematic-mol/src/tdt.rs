//! Streaming Daylight TDT (Tagged Data) file I/O — `.tdt`.
//!
//! A TDT record is a sequence of `TAGNAME<value>` lines terminated by a line
//! starting with `|`. This is chematic's counterpart to RDKit's
//! `TDTMolSupplier`/`TDTWriter`
//! (`Code/GraphMol/FileParsers/TDTMolSupplier.cpp`/`TDTWriter.cpp`, RDKit
//! commit `8afba32ec539dcb2369bc84549d802aca3f7eb39`, the true resolution of
//! tag `Release_2026_03_4`).
//!
//! ```text
//! $SMI<CCO>
//! NAME<ethanol>
//! ACTIVITY<7.2>
//! |
//! ```
//!
//! **Grammar, source-confirmed (not guessed):** a record starts with the
//! literal 5 characters `$SMI<` at column 0. A generic tag line is
//! `TAGNAME<value>`: the tag name is the substring before the *first* `<`
//! (trimmed of leading/trailing space/tab only); the value is the substring
//! between that first `<` and the *last* `>` on the line — both `<` and `>`
//! must appear on the same physical line for a generic tag, or it's a parse
//! error. Only the `2D`/`3D` coordinate tags support multi-line values
//! (comma-separated numbers terminated by a token containing `;>`). The
//! record terminator `|` must be the first character of its line; nothing
//! after it is inspected. There is no escaping mechanism for `<`/`>` inside
//! a value anywhere in this grammar (RDKit has none either).
//!
//! **Deliberate, documented divergences from RDKit, all found via this
//! session's source audit and confirmed against a live RDKit oracle:**
//!
//! - **Coordinate-block parsing bug, fixed, not replicated.** RDKit's own
//!   `TDTMolSupplier` silently drops the *last* atom's position when
//!   reading a `2D`/`3D` coordinate tag (its comma-tokenizer treats the
//!   token containing the trailing `;>` as "found the terminator," and
//!   never pushes that token's own numeric value) — confirmed against a
//!   live RDKit 2026.03.3 run, not merely read from source. This chematic
//!   port parses the full coordinate list correctly. Precision-first
//!   library, real bug found in the reference implementation, not
//!   reproduced.
//! - **Reader/writer name-tag symmetry, fixed, not replicated.** RDKit's
//!   `TDTWriter` hard-codes its name-tag output as literally `"NAME"`
//!   (unconfigurable), while `TDTMolSupplier`'s own `nameRecord` parameter
//!   defaults to `""` (meaning *no* tag populates `_Name`) — so a bare
//!   `TDTWriter`/`TDTMolSupplier()` round trip silently loses the molecule
//!   name by default, confirmed empirically against real RDKit. This is a
//!   genuine, low-value RDKit design wart, not a deliberate feature;
//!   [`TdtReaderOptions::name_tag`]/[`TdtWriterOptions::name_tag`] both
//!   default to `Some("NAME")`, so a default reader+writer pair round-trips
//!   the name correctly out of the box.
//! - **Malformed-record recovery, hardened, not replicated.** In real
//!   RDKit, a missing `>` on a generic tag line throws an exception that is
//!   *not* caught inside `TDTMolSupplier::next()`'s own position-advance
//!   bookkeeping — confirmed against source and empirically: naively
//!   retrying `next()` after catching the exception re-seeks to the exact
//!   same failed record and throws again, indefinitely. This reader instead
//!   scans forward to the next `$SMI<`/record terminator internally before
//!   returning `Err` for the broken record, so the *next* call to
//!   [`Iterator::next`] always makes progress — matching this crate's
//!   established "malformed input never causes non-termination" discipline
//!   (see [`crate::smiles_table`]'s own `strict_parsing` design).
//! - **`2D`/`3D` tags require an explicit opt-in to be interpreted as
//!   coordinates** ([`TdtReaderOptions::read_2d`]/`read_3d`, both default
//!   `false`) — matching RDKit's own default of *not* reading a conformer
//!   (RDKit's Python binding explicitly pins `confId3D=-1` regardless of
//!   the underlying C++ class's own differing default, a cross-entry-point
//!   inconsistency this port does not reproduce). When disabled, a `2D`/`3D`
//!   tag is stored as a generic string property holding the raw,
//!   unparsed coordinate text — never silently discarded.
//! - **A file with no trailing newline drops its final tag line in real
//!   RDKit, discovered against a live oracle this session (not predicted by
//!   the initial source audit).** A record whose last generic tag line is
//!   also the file's last byte, with no trailing `\n`, has that tag line
//!   silently dropped by RDKit's `TDTMolSupplier` — confirmed empirically
//!   (`$SMI<CC>\nNAME<ethane>` with no trailing newline yields `_Name == ""`
//!   in real RDKit, not `"ethane"`). This reader's line-reading loop
//!   (`BufRead::read_line`) returns a final unterminated line correctly, so
//!   chematic does not reproduce this gap.

use std::io::{BufRead, Write};

use chematic_core::Molecule;
use chematic_smiles::SmilesError;

use crate::record::MoleculeRecord;

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

/// Configuration for [`TdtRecordReader`].
#[derive(Debug, Clone)]
pub struct TdtReaderOptions {
    /// Which tag's value populates [`MoleculeRecord::name`]. `None` means no
    /// tag does (RDKit's own `nameRecord=""` sentinel).
    pub name_tag: Option<String>,
    /// Interpret a `2D` tag as a 2D coordinate list (see the module docs for
    /// why this defaults to `false`, matching RDKit's real Python default).
    pub read_2d: bool,
    /// Interpret a `3D` tag as a 3D coordinate list.
    pub read_3d: bool,
    /// `true`: a malformed record yields one `Err`, then the reader stops.
    /// `false`: the same `Err` is yielded, but the reader recovers and
    /// continues to the next record (always, regardless of this flag --
    /// see the module docs' malformed-record-recovery note -- this flag
    /// only controls whether the *caller* sees further items afterward).
    pub strict_parsing: bool,
    /// Hard cap on a single line's byte length.
    pub max_line_bytes: usize,
    /// Maximum records yielded by the reader.
    pub max_records: usize,
    /// Maximum tags retained in one record.
    pub max_tags_per_record: usize,
}

impl Default for TdtReaderOptions {
    fn default() -> Self {
        Self {
            name_tag: Some("NAME".to_string()),
            read_2d: false,
            read_3d: false,
            strict_parsing: false,
            max_line_bytes: 1 << 20,
            max_records: 100_000,
            max_tags_per_record: 10_000,
        }
    }
}

/// Configuration for [`TdtRecordWriter`].
#[derive(Debug, Clone)]
pub struct TdtWriterOptions {
    /// Tag name under which [`MoleculeRecord::name`] is written, or `None`
    /// to omit the name entirely. Matches RDKit's `TDTWriter` default of
    /// writing a `"NAME"` tag when a name is present.
    pub name_tag: Option<String>,
    /// Write a `2D` coordinate tag from [`MoleculeRecord::coordinates_2d`],
    /// if present.
    pub write_2d: bool,
    /// Write a `3D` coordinate tag from [`MoleculeRecord::coordinates_3d`],
    /// if present.
    pub write_3d: bool,
    /// Which properties to write, and in what order. `None` writes every
    /// property on the record (matches RDKit's own `TDTWriter` default of
    /// writing all non-computed properties); `Some(keys)` restricts to
    /// exactly those keys, in that order (matches `SetProps`).
    pub properties: Option<Vec<String>>,
    /// Number of significant digits for written coordinate values (matches
    /// `TDTWriter::SetNumDigits`, default 4).
    pub precision: usize,
}

impl Default for TdtWriterOptions {
    fn default() -> Self {
        Self {
            name_tag: Some("NAME".to_string()),
            write_2d: false,
            write_3d: true, // matches RDKit's TDTWriter default (3D, not 2D)
            properties: None,
            precision: 4,
        }
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors from [`TdtRecordReader`]/[`TdtRecordWriter`].
#[derive(Debug, Clone, PartialEq)]
pub enum TdtError {
    /// The record did not start with `$SMI<...>`.
    MissingSmilesTag { line: usize, record_index: usize },
    /// The `$SMI` tag's value could not be parsed as SMILES.
    InvalidSmiles {
        line: usize,
        record_index: usize,
        detail: SmilesError,
    },
    /// A generic tag line had `<` but no closing `>` on the same line.
    UnterminatedTag {
        line: usize,
        record_index: usize,
        tag_name: String,
    },
    /// A `2D`/`3D` coordinate list could not be parsed as comma-separated
    /// numbers, or its length didn't match `3 * atom_count`/`2 * atom_count`.
    MalformedCoordinateList {
        line: usize,
        record_index: usize,
        tag_name: String,
        detail: String,
    },
    /// A line exceeded `max_line_bytes`.
    LineTooLong {
        line: usize,
        record_index: usize,
        limit: usize,
    },
    /// A record contained more tags than allowed.
    TooManyTags {
        line: usize,
        record_index: usize,
        actual: usize,
        limit: usize,
    },
    /// The reader reached its maximum record budget.
    TooManyRecords { limit: usize },
    /// An IO error occurred while reading the input stream.
    Io(String),
}

impl std::fmt::Display for TdtError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingSmilesTag { line, record_index } => {
                write!(
                    f,
                    "record {record_index} at line {line}: missing $SMI<...> tag"
                )
            }
            Self::InvalidSmiles {
                line,
                record_index,
                detail,
            } => write!(
                f,
                "invalid SMILES at line {line} (record {record_index}): {detail}"
            ),
            Self::UnterminatedTag {
                line,
                record_index,
                tag_name,
            } => write!(
                f,
                "no closing '>' for tag {tag_name:?} at line {line} (record {record_index})"
            ),
            Self::MalformedCoordinateList {
                line,
                record_index,
                tag_name,
                detail,
            } => write!(
                f,
                "malformed {tag_name} coordinate list at line {line} (record {record_index}): {detail}"
            ),
            Self::LineTooLong {
                line,
                record_index,
                limit,
            } => write!(
                f,
                "line {line} (record {record_index}) exceeds the {limit}-byte limit"
            ),
            Self::TooManyTags {
                line,
                record_index,
                actual,
                limit,
            } => write!(
                f,
                "record {record_index} at line {line} has {actual} tags, exceeding the {limit}-tag limit"
            ),
            Self::TooManyRecords { limit } => {
                write!(f, "TDT input exceeds the {limit}-record limit")
            }
            Self::Io(msg) => write!(f, "IO error: {msg}"),
        }
    }
}

impl std::error::Error for TdtError {}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

enum RawTag {
    /// A generic `TAGNAME<value>` tag, fully resolved on one line.
    Generic { name: String, value: String },
    /// A `2D`/`3D` coordinate tag, already collected across however many
    /// physical lines its `;>`-terminated number list spanned.
    Coordinates { name: String, raw: String },
}

/// Streaming iterator over a TDT file.
pub struct TdtRecordReader<R: BufRead> {
    reader: R,
    options: TdtReaderOptions,
    line_number: usize,
    record_index: usize,
    stopped: bool,
    at_eof: bool,
}

impl<R: BufRead> TdtRecordReader<R> {
    pub fn new(reader: R, options: TdtReaderOptions) -> Self {
        Self {
            reader,
            options,
            line_number: 0,
            record_index: 0,
            stopped: false,
            at_eof: false,
        }
    }

    fn read_raw_line(&mut self) -> Result<Option<String>, TdtError> {
        let mut buf = String::new();
        let n = self
            .reader
            .read_line(&mut buf)
            .map_err(|e| TdtError::Io(e.to_string()))?;
        if n == 0 {
            self.at_eof = true;
            return Ok(None);
        }
        self.line_number += 1;
        if buf.len() > self.options.max_line_bytes {
            return Err(TdtError::LineTooLong {
                line: self.line_number,
                record_index: self.record_index,
                limit: self.options.max_line_bytes,
            });
        }
        while buf.ends_with('\n') || buf.ends_with('\r') {
            buf.pop();
        }
        Ok(Some(buf))
    }

    /// Skip forward to (but not past) the next `$SMI<` line, recovering
    /// from a malformed record. Consumes lines up to and including a `|`
    /// terminator if one is seen first, or up to (excluding) the next
    /// `$SMI<` line.
    fn recover_to_next_record(&mut self) -> Result<(), TdtError> {
        loop {
            match self.reader.fill_buf() {
                Ok([]) => return Ok(()), // EOF
                Ok(_) => {}
                Err(e) => return Err(TdtError::Io(e.to_string())),
            }
            // Peek one line without consuming past a `$SMI<` boundary.
            let mut probe = String::new();
            let consumed = {
                let buf = self
                    .reader
                    .fill_buf()
                    .map_err(|e| TdtError::Io(e.to_string()))?;
                let text = String::from_utf8_lossy(buf);
                if text.starts_with("$SMI<") {
                    0
                } else if let Some(nl) = text.find('\n') {
                    probe.push_str(&text[..=nl]);
                    nl + 1
                } else {
                    let len = text.len();
                    probe.push_str(&text);
                    len
                }
            };
            if consumed == 0 {
                return Ok(());
            }
            self.reader.consume(consumed);
            self.line_number += 1;
            let is_terminator = probe.trim_end_matches(['\r', '\n']).starts_with('|');
            if is_terminator {
                return Ok(());
            }
        }
    }

    fn read_record_tags(&mut self) -> Result<Option<(usize, Vec<RawTag>)>, TdtError> {
        let start_line = loop {
            let line = match self.read_raw_line()? {
                None => return Ok(None),
                Some(l) => l,
            };
            if line.trim().is_empty() {
                continue;
            }
            break line;
        };

        let record_line = self.line_number;

        if !start_line.starts_with("$SMI<") {
            // A real record's first non-blank line must be `$SMI<...>`. This
            // is a malformed record (missing SMILES entirely), not noise to
            // skip -- recover to the next record boundary so the *next*
            // `next()` call still makes progress.
            self.recover_to_next_record()?;
            return Err(TdtError::MissingSmilesTag {
                line: record_line,
                record_index: self.record_index,
            });
        }

        // Any error from here on (malformed tag, oversized line, malformed
        // coordinate list, IO error mid-record) must recover to the next
        // record boundary before propagating -- otherwise a leftover
        // fragment of THIS record (e.g. its `|` terminator, already read
        // into a buffer that got rejected for being too long) is
        // misinterpreted as the start of a phantom next record. Centralized
        // here rather than at each individual error site inside
        // `read_record_body`, so no future error path can forget it.
        match self.read_record_body(&start_line, record_line) {
            Ok(tags) => Ok(Some((record_line, tags))),
            Err(e) => {
                self.recover_to_next_record()?;
                Err(e)
            }
        }
    }

    fn push_tag(&self, tags: &mut Vec<RawTag>, tag: RawTag) -> Result<(), TdtError> {
        if tags.len() >= self.options.max_tags_per_record {
            return Err(TdtError::TooManyTags {
                line: self.line_number,
                record_index: self.record_index,
                actual: tags.len() + 1,
                limit: self.options.max_tags_per_record,
            });
        }
        tags.push(tag);
        Ok(())
    }

    fn read_record_body(
        &mut self,
        start_line: &str,
        record_line: usize,
    ) -> Result<Vec<RawTag>, TdtError> {
        let mut tags = Vec::new();

        let smi_value = extract_tag_value(start_line).ok_or_else(|| TdtError::UnterminatedTag {
            line: record_line,
            record_index: self.record_index,
            tag_name: "$SMI".to_string(),
        })?;
        self.push_tag(
            &mut tags,
            RawTag::Generic {
                name: "$SMI".to_string(),
                value: smi_value,
            },
        )?;

        loop {
            let line = match self.read_raw_line()? {
                None => break, // EOF mid-record: return what we have (matches RDKit)
                Some(l) => l,
            };
            if line.starts_with('|') {
                break;
            }
            if line.trim().is_empty() {
                continue;
            }

            let tag_name = tag_name_of(&line);
            if tag_name == "2D" || tag_name == "3D" {
                let raw = self.read_coordinate_list(&line, &tag_name)?;
                self.push_tag(
                    &mut tags,
                    RawTag::Coordinates {
                        name: tag_name,
                        raw,
                    },
                )?;
            } else {
                match extract_tag_value(&line) {
                    Some(value) => self.push_tag(
                        &mut tags,
                        RawTag::Generic {
                            name: tag_name,
                            value,
                        },
                    )?,
                    None => {
                        return Err(TdtError::UnterminatedTag {
                            line: self.line_number,
                            record_index: self.record_index,
                            tag_name,
                        });
                    }
                }
            }
        }

        Ok(tags)
    }

    /// Read a coordinate tag's value, following continuation lines until a
    /// token containing `;>` is seen (matches RDKit's real multi-line
    /// number-list grammar for `2D`/`3D` only).
    fn read_coordinate_list(
        &mut self,
        first_line: &str,
        tag_name: &str,
    ) -> Result<String, TdtError> {
        let start = first_line
            .find('<')
            .ok_or_else(|| TdtError::UnterminatedTag {
                line: self.line_number,
                record_index: self.record_index,
                tag_name: tag_name.to_string(),
            })?;
        let mut buf = first_line[start + 1..].to_string();
        while !buf.contains(";>") {
            match self.read_raw_line()? {
                None => {
                    return Err(TdtError::MalformedCoordinateList {
                        line: self.line_number,
                        record_index: self.record_index,
                        tag_name: tag_name.to_string(),
                        detail: "unterminated coordinate list at EOF".to_string(),
                    });
                }
                Some(more) => {
                    buf.push(',');
                    buf.push_str(&more);
                }
            }
        }
        let end = buf.find(";>").unwrap();
        Ok(buf[..end].to_string())
    }
}

fn tag_name_of(line: &str) -> String {
    match line.find('<') {
        Some(pos) => line[..pos].trim().to_string(),
        None => line.trim().to_string(),
    }
}

/// Extract a generic tag's value: substring between the first `<` and the
/// last `>` on the line. `None` if either delimiter is missing.
fn extract_tag_value(line: &str) -> Option<String> {
    let start = line.find('<')?;
    let end = line.rfind('>')?;
    if end < start {
        return None;
    }
    Some(line[start + 1..end].to_string())
}

fn parse_coordinate_list(raw: &str) -> Result<Vec<f64>, String> {
    raw.split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.trim().parse::<f64>().map_err(|e| e.to_string()))
        .collect()
}

impl<R: BufRead> Iterator for TdtRecordReader<R> {
    type Item = Result<MoleculeRecord, TdtError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.stopped {
            return None;
        }
        if self.record_index >= self.options.max_records {
            self.stopped = true;
            return Some(Err(TdtError::TooManyRecords {
                limit: self.options.max_records,
            }));
        }

        let (record_line, tags) = match self.read_record_tags() {
            Ok(None) => return None,
            Ok(Some(t)) => t,
            Err(e) => {
                let record_index = self.record_index;
                self.record_index += 1;
                if self.options.strict_parsing {
                    self.stopped = true;
                }
                return Some(Err(match e {
                    TdtError::UnterminatedTag { line, tag_name, .. } => TdtError::UnterminatedTag {
                        line,
                        record_index,
                        tag_name,
                    },
                    other => other,
                }));
            }
        };

        let record_index = self.record_index;
        self.record_index += 1;

        let mut smi: Option<String> = None;
        let mut properties: Vec<(String, String)> = Vec::new();
        let mut name = String::new();
        let mut coordinates_2d: Option<Vec<[f64; 2]>> = None;
        let mut coordinates_3d: Option<Vec<[f64; 3]>> = None;
        let mut raw_2d: Option<String> = None;
        let mut raw_3d: Option<String> = None;

        for tag in tags {
            match tag {
                RawTag::Generic {
                    name: tag_name,
                    value,
                } if tag_name == "$SMI" => {
                    smi = Some(value);
                }
                RawTag::Generic {
                    name: tag_name,
                    value,
                } => {
                    if self.options.name_tag.as_deref() == Some(tag_name.as_str()) {
                        name = value.clone();
                    }
                    set_property(&mut properties, tag_name, value);
                }
                RawTag::Coordinates {
                    name: tag_name,
                    raw,
                } => {
                    if tag_name == "2D" {
                        raw_2d = Some(raw);
                    } else {
                        raw_3d = Some(raw);
                    }
                }
            }
        }

        let Some(smi) = smi else {
            return Some(Err(TdtError::MissingSmilesTag {
                line: record_line,
                record_index,
            }));
        };

        let mol: Molecule = match chematic_smiles::parse(&smi) {
            Ok(m) => m,
            Err(detail) => {
                let err = TdtError::InvalidSmiles {
                    line: record_line,
                    record_index,
                    detail,
                };
                if self.options.strict_parsing {
                    self.stopped = true;
                }
                return Some(Err(err));
            }
        };

        if let Some(raw) = raw_2d {
            if self.options.read_2d {
                match parse_coordinate_list(&raw) {
                    Ok(nums) if nums.len() == 2 * mol.atom_count() => {
                        coordinates_2d = Some(
                            nums.as_chunks::<2>()
                                .0
                                .iter()
                                .map(|c| [c[0], c[1]])
                                .collect(),
                        );
                    }
                    Ok(nums) => {
                        return Some(Err(TdtError::MalformedCoordinateList {
                            line: record_line,
                            record_index,
                            tag_name: "2D".to_string(),
                            detail: format!(
                                "expected {} numbers (2 * {} atoms), found {}",
                                2 * mol.atom_count(),
                                mol.atom_count(),
                                nums.len()
                            ),
                        }));
                    }
                    Err(detail) => {
                        return Some(Err(TdtError::MalformedCoordinateList {
                            line: record_line,
                            record_index,
                            tag_name: "2D".to_string(),
                            detail,
                        }));
                    }
                }
            } else {
                set_property(&mut properties, "2D".to_string(), raw);
            }
        }
        if let Some(raw) = raw_3d {
            if self.options.read_3d {
                match parse_coordinate_list(&raw) {
                    Ok(nums) if nums.len() == 3 * mol.atom_count() => {
                        coordinates_3d = Some(
                            nums.as_chunks::<3>()
                                .0
                                .iter()
                                .map(|c| [c[0], c[1], c[2]])
                                .collect(),
                        );
                    }
                    Ok(nums) => {
                        return Some(Err(TdtError::MalformedCoordinateList {
                            line: record_line,
                            record_index,
                            tag_name: "3D".to_string(),
                            detail: format!(
                                "expected {} numbers (3 * {} atoms), found {}",
                                3 * mol.atom_count(),
                                mol.atom_count(),
                                nums.len()
                            ),
                        }));
                    }
                    Err(detail) => {
                        return Some(Err(TdtError::MalformedCoordinateList {
                            line: record_line,
                            record_index,
                            tag_name: "3D".to_string(),
                            detail,
                        }));
                    }
                }
            } else {
                set_property(&mut properties, "3D".to_string(), raw);
            }
        }

        Some(Ok(MoleculeRecord {
            mol,
            name,
            properties,
            coordinates_2d,
            coordinates_3d,
        }))
    }
}

/// Insert-or-overwrite-in-place, matching RDKit's own `Dict::setVal`
/// "last wins, same position" semantics for repeated tags.
fn set_property(properties: &mut Vec<(String, String)>, key: String, value: String) {
    if let Some(existing) = properties.iter_mut().find(|(k, _)| *k == key) {
        existing.1 = value;
    } else {
        properties.push((key, value));
    }
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Streaming writer for a TDT file.
pub struct TdtRecordWriter<W: Write> {
    writer: W,
    options: TdtWriterOptions,
    wrote_any: bool,
    closed: bool,
}

impl<W: Write> TdtRecordWriter<W> {
    pub fn new(writer: W, options: TdtWriterOptions) -> Self {
        Self {
            writer,
            options,
            wrote_any: false,
            closed: false,
        }
    }

    /// Write one record. The molecule's SMILES is serialized via
    /// [`chematic_smiles::write`].
    pub fn write_record(&mut self, record: &MoleculeRecord) -> std::io::Result<()> {
        if self.wrote_any {
            writeln!(self.writer, "|")?;
        }
        self.wrote_any = true;

        let smi = chematic_smiles::write(&record.mol);
        writeln!(self.writer, "$SMI<{smi}>")?;

        if let Some(tag) = &self.options.name_tag
            && !record.name.is_empty()
        {
            writeln!(self.writer, "{tag}<{}>", record.name)?;
        }

        if self.options.write_2d
            && let Some(coords) = &record.coordinates_2d
        {
            self.write_coordinate_tag("2D", coords.iter().flat_map(|c| c.iter().copied()))?;
        }
        if self.options.write_3d
            && let Some(coords) = &record.coordinates_3d
        {
            self.write_coordinate_tag("3D", coords.iter().flat_map(|c| c.iter().copied()))?;
        }

        match &self.options.properties {
            None => {
                for (k, v) in &record.properties {
                    writeln!(self.writer, "{k}<{}>", v.replace('\n', " "))?;
                }
            }
            Some(keys) => {
                for k in keys {
                    if let Some(v) = record.get_property(k) {
                        writeln!(self.writer, "{k}<{}>", v.replace('\n', " "))?;
                    }
                }
            }
        }

        Ok(())
    }

    fn write_coordinate_tag(
        &mut self,
        tag: &str,
        values: impl Iterator<Item = f64>,
    ) -> std::io::Result<()> {
        let precision = self.options.precision;
        let formatted: Vec<String> = values.map(|v| format!("{v:.precision$}")).collect();
        writeln!(self.writer, "{tag}<{};>", formatted.join(","))
    }

    /// Terminate the final record (matches RDKit's `TDTWriter::close`,
    /// which writes the trailing `|` only on close, not on every `write`).
    /// Idempotent -- a second call is a no-op, so an explicit `close()`
    /// followed by `Drop` never double-writes the terminator.
    pub fn close(&mut self) -> std::io::Result<()> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        if self.wrote_any {
            writeln!(self.writer, "|")?;
        }
        self.writer.flush()
    }

    /// Consume the writer, returning the underlying `W` without finalizing
    /// (call [`Self::close`] first if the trailing `|` terminator matters --
    /// this writer does not finalize on drop, matching this crate's other
    /// streaming writers).
    pub fn into_inner(self) -> W {
        self.writer
    }

    /// Replace the set (and order) of properties written per record.
    /// `None` writes every property (RDKit's own `TDTWriter` default).
    pub fn set_properties(&mut self, properties: Option<Vec<String>>) {
        self.options.properties = properties;
    }

    /// Replace the tag under which the record name is written, or `None`
    /// to omit it entirely (matches `SetWriteNames`).
    pub fn set_name_tag(&mut self, name_tag: Option<String>) {
        self.options.name_tag = name_tag;
    }

    /// Replace the number of significant digits used for written coordinate
    /// values (matches `SetNumDigits`).
    pub fn set_precision(&mut self, precision: usize) {
        self.options.precision = precision;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn reader_over(
        input: &str,
        options: TdtReaderOptions,
    ) -> TdtRecordReader<BufReader<Cursor<Vec<u8>>>> {
        TdtRecordReader::new(
            BufReader::new(Cursor::new(input.as_bytes().to_vec())),
            options,
        )
    }

    #[test]
    fn basic_record_with_name_and_property() {
        let input = "$SMI<CCO>\nNAME<ethanol>\nACTIVITY<7.2>\n|\n";
        let rec = reader_over(input, TdtReaderOptions::default())
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(rec.mol.atom_count(), 3);
        assert_eq!(rec.name, "ethanol");
        assert_eq!(rec.get_property("ACTIVITY"), Some("7.2"));
    }

    #[test]
    fn multiple_records() {
        let input = "$SMI<CC>\nNAME<ethane>\n|\n$SMI<CCO>\nNAME<ethanol>\n|\n";
        let recs: Vec<_> = reader_over(input, TdtReaderOptions::default())
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "ethane");
        assert_eq!(recs[1].name, "ethanol");
    }

    #[test]
    fn missing_smi_tag_is_explicit_error() {
        let input = "NAME<no_smiles>\n|\n";
        let mut reader = reader_over(input, TdtReaderOptions::default());
        assert!(matches!(
            reader.next(),
            Some(Err(TdtError::MissingSmilesTag { .. }))
        ));
    }

    #[test]
    fn empty_value_tag_is_stored_as_empty_string() {
        let input = "$SMI<CC>\nNOTE<>\n|\n";
        let rec = reader_over(input, TdtReaderOptions::default())
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(rec.get_property("NOTE"), Some(""));
    }

    #[test]
    fn repeated_tag_last_wins_same_position() {
        let input = "$SMI<CC>\nFOO<first>\nBAR<x>\nFOO<second>\n|\n";
        let rec = reader_over(input, TdtReaderOptions::default())
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(rec.get_property("FOO"), Some("second"));
        assert_eq!(rec.properties[0].0, "FOO"); // position preserved, not moved to end
    }

    #[test]
    fn unknown_tags_stored_as_generic_properties() {
        let input = "$SMI<CC>\nMFCD<12345>\nCAS<64-17-5>\n|\n";
        let rec = reader_over(input, TdtReaderOptions::default())
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(rec.get_property("MFCD"), Some("12345"));
        assert_eq!(rec.get_property("CAS"), Some("64-17-5"));
    }

    #[test]
    fn unterminated_generic_tag_recovers_to_next_record() {
        let input = "$SMI<CC>\nBROKEN<no_close\n|\n$SMI<CCO>\nNAME<ethanol>\n|\n";
        let mut reader = reader_over(input, TdtReaderOptions::default());
        assert!(matches!(
            reader.next(),
            Some(Err(TdtError::UnterminatedTag { .. }))
        ));
        let second = reader.next().unwrap().unwrap();
        assert_eq!(second.name, "ethanol");
        assert!(reader.next().is_none());
    }

    #[test]
    fn strict_parsing_stops_after_first_error() {
        let opts = TdtReaderOptions {
            strict_parsing: true,
            ..Default::default()
        };
        let input = "NAME<no_smiles>\n|\n$SMI<CC>\n|\n";
        let mut reader = reader_over(input, opts);
        assert!(matches!(
            reader.next(),
            Some(Err(TdtError::MissingSmilesTag { .. }))
        ));
        assert!(reader.next().is_none());
    }

    #[test]
    fn eof_mid_record_no_terminator_still_succeeds() {
        let input = "$SMI<CC>\nNAME<ethane>";
        let rec = reader_over(input, TdtReaderOptions::default())
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(rec.name, "ethane");
    }

    #[test]
    fn coordinate_2d_opt_in_parses_full_list_no_last_atom_drop() {
        let opts = TdtReaderOptions {
            read_2d: true,
            ..Default::default()
        };
        // 3-atom molecule: 6 numbers.
        let input = "$SMI<CCO>\n2D<0.0,0.0,1.0,0.0,2.0,1.0;>\n|\n";
        let rec = reader_over(input, opts).next().unwrap().unwrap();
        let coords = rec.coordinates_2d.expect("2D coords present");
        assert_eq!(coords.len(), 3);
        assert_eq!(coords[2], [2.0, 1.0]); // last atom NOT dropped (RDKit's own bug, not replicated)
    }

    #[test]
    fn coordinate_3d_opt_in_parses_full_list() {
        let opts = TdtReaderOptions {
            read_3d: true,
            ..Default::default()
        };
        let input = "$SMI<CCO>\n3D<0.0,0.0,0.0,1.0,0.0,0.0,2.0,1.0,0.5;>\n|\n";
        let rec = reader_over(input, opts).next().unwrap().unwrap();
        let coords = rec.coordinates_3d.expect("3D coords present");
        assert_eq!(coords.len(), 3);
        assert_eq!(coords[2], [2.0, 1.0, 0.5]);
    }

    #[test]
    fn coordinate_tag_without_opt_in_stored_as_raw_property() {
        let input = "$SMI<CCO>\n2D<0.0,0.0,1.0,0.0,2.0,1.0;>\n|\n";
        let rec = reader_over(input, TdtReaderOptions::default())
            .next()
            .unwrap()
            .unwrap();
        assert!(rec.coordinates_2d.is_none());
        assert_eq!(rec.get_property("2D"), Some("0.0,0.0,1.0,0.0,2.0,1.0"));
    }

    #[test]
    fn multiline_coordinate_list_spanning_lines() {
        let opts = TdtReaderOptions {
            read_2d: true,
            ..Default::default()
        };
        let input = "$SMI<CCO>\n2D<0.0,0.0,\n1.0,0.0,\n2.0,1.0;>\n|\n";
        let rec = reader_over(input, opts).next().unwrap().unwrap();
        let coords = rec.coordinates_2d.expect("2D coords present");
        assert_eq!(coords, vec![[0.0, 0.0], [1.0, 0.0], [2.0, 1.0]]);
    }

    #[test]
    fn writer_default_roundtrip_preserves_name() {
        let mol = chematic_smiles::parse("CCO").unwrap();
        let mut record = MoleculeRecord::new(mol);
        record.name = "ethanol".to_string();
        record
            .properties
            .push(("Activity".to_string(), "7.2".to_string()));

        let mut out = Vec::new();
        {
            let mut writer = TdtRecordWriter::new(&mut out, TdtWriterOptions::default());
            writer.write_record(&record).unwrap();
            writer.close().unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with("$SMI<"));
        assert!(text.contains("NAME<ethanol>"));
        assert!(text.contains("Activity<7.2>"));
        assert!(text.trim_end().ends_with('|'));

        let mut reader = reader_over(&text, TdtReaderOptions::default());
        let back = reader.next().unwrap().unwrap();
        assert_eq!(back.mol.atom_count(), 3);
        assert_eq!(back.name, "ethanol"); // default reader+writer round-trips the name
    }

    #[test]
    fn writer_multiple_records_terminator_placement() {
        let mol1 = chematic_smiles::parse("CC").unwrap();
        let mol2 = chematic_smiles::parse("CCO").unwrap();
        let mut out = Vec::new();
        {
            let mut writer = TdtRecordWriter::new(&mut out, TdtWriterOptions::default());
            writer.write_record(&MoleculeRecord::new(mol1)).unwrap();
            writer.write_record(&MoleculeRecord::new(mol2)).unwrap();
            writer.close().unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let terminator_count = text.lines().filter(|l| *l == "|").count();
        assert_eq!(terminator_count, 2);
    }

    #[test]
    fn writer_properties_none_writes_all() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let mut record = MoleculeRecord::new(mol);
        record.properties = vec![
            ("A".to_string(), "1".to_string()),
            ("B".to_string(), "2".to_string()),
        ];
        let mut out = Vec::new();
        {
            let mut writer = TdtRecordWriter::new(
                &mut out,
                TdtWriterOptions {
                    name_tag: None,
                    ..Default::default()
                },
            );
            writer.write_record(&record).unwrap();
            writer.close().unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("A<1>"));
        assert!(text.contains("B<2>"));
    }

    #[test]
    fn writer_properties_selected_subset_and_order() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let mut record = MoleculeRecord::new(mol);
        record.properties = vec![
            ("A".to_string(), "1".to_string()),
            ("B".to_string(), "2".to_string()),
        ];
        let mut out = Vec::new();
        {
            let mut writer = TdtRecordWriter::new(
                &mut out,
                TdtWriterOptions {
                    name_tag: None,
                    properties: Some(vec!["B".to_string()]),
                    ..Default::default()
                },
            );
            writer.write_record(&record).unwrap();
            writer.close().unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("B<2>"));
        assert!(!text.contains("A<1>"));
    }

    #[test]
    fn writer_embedded_newline_in_property_replaced_with_space() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let mut record = MoleculeRecord::new(mol);
        record
            .properties
            .push(("Note".to_string(), "line1\nline2".to_string()));
        let mut out = Vec::new();
        {
            let mut writer = TdtRecordWriter::new(
                &mut out,
                TdtWriterOptions {
                    name_tag: None,
                    ..Default::default()
                },
            );
            writer.write_record(&record).unwrap();
            writer.close().unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("Note<line1 line2>"));
    }
}

// ---------------------------------------------------------------------------
// Adversarial / fuzz-style tests
// ---------------------------------------------------------------------------
//
// Same rationale as `crate::smiles_table::adversarial_tests`: no
// cargo-fuzz/libfuzzer harness exists anywhere in this workspace, and
// introducing one for this module was judged disproportionate. Deterministic
// adversarial unit tests instead, all asserting only "no panic/hang/OOM".

#[cfg(test)]
mod adversarial_tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn drain(input: &[u8], options: TdtReaderOptions) -> usize {
        TdtRecordReader::new(BufReader::new(Cursor::new(input.to_vec())), options).count()
    }

    #[test]
    fn empty_input_no_panic() {
        assert_eq!(drain(b"", TdtReaderOptions::default()), 0);
    }

    #[test]
    fn truncated_input_mid_smi_tag_no_panic() {
        assert_eq!(drain(b"$SMI<CC", TdtReaderOptions::default()), 1);
    }

    #[test]
    fn truncated_input_at_dollar_sign_only_no_panic() {
        assert_eq!(drain(b"$", TdtReaderOptions::default()), 1);
    }

    #[test]
    fn huge_line_is_explicit_error_not_panic_or_oom() {
        let mut line = b"$SMI<CC>\nNOTE<".to_vec();
        line.extend(std::iter::repeat_n(b'a', 10_000_000));
        line.push(b'>');
        line.push(b'\n');
        line.extend_from_slice(b"|\n");
        let opts = TdtReaderOptions {
            max_line_bytes: 1 << 20,
            ..Default::default()
        };
        let reader = TdtRecordReader::new(BufReader::new(Cursor::new(line)), opts);
        let results: Vec<_> = reader.collect();
        assert_eq!(results.len(), 1);
        assert!(matches!(results[0], Err(TdtError::LineTooLong { .. })));
    }

    #[test]
    fn record_and_tag_limits_are_explicit_errors() {
        let opts = TdtReaderOptions {
            max_records: 1,
            ..Default::default()
        };
        let input = b"$SMI<CC>\n|\n$SMI<CCO>\n|\n";
        let mut reader = TdtRecordReader::new(BufReader::new(Cursor::new(input)), opts);
        assert!(reader.next().unwrap().is_ok());
        assert!(matches!(
            reader.next(),
            Some(Err(TdtError::TooManyRecords { limit: 1 }))
        ));

        let opts = TdtReaderOptions {
            max_tags_per_record: 1,
            ..Default::default()
        };
        let mut reader = TdtRecordReader::new(
            BufReader::new(Cursor::new(b"$SMI<CC>\nNAME<ethane>\n|\n")),
            opts,
        );
        assert!(matches!(
            reader.next(),
            Some(Err(TdtError::TooManyTags { limit: 1, .. }))
        ));
    }

    #[test]
    fn huge_property_value_within_limit_does_not_panic() {
        let mut line = b"$SMI<CC>\nNOTE<".to_vec();
        line.extend(std::iter::repeat_n(b'x', 500_000));
        line.extend_from_slice(b">\n|\n");
        let opts = TdtReaderOptions {
            max_line_bytes: 1 << 21,
            ..Default::default()
        };
        let results: Vec<_> =
            TdtRecordReader::new(BufReader::new(Cursor::new(line)), opts).collect();
        assert_eq!(results.len(), 1);
        assert!(results[0].is_ok());
    }

    #[test]
    fn invalid_utf8_byte_path_yields_io_error_not_panic() {
        let input: &[u8] = b"$SMI<CC>\nNOTE<\xFF\xFE>\n|\n";
        let results: Vec<_> = TdtRecordReader::new(
            BufReader::new(Cursor::new(input.to_vec())),
            TdtReaderOptions::default(),
        )
        .collect();
        assert_eq!(results.len(), 1);
        assert!(matches!(results[0], Err(TdtError::Io(_))));
    }

    #[test]
    fn excessive_property_count_does_not_panic() {
        let mut line = String::from("$SMI<CC>\n");
        for i in 0..5000 {
            line.push_str(&format!("P{i}<v{i}>\n"));
        }
        line.push_str("|\n");
        let results: Vec<_> = TdtRecordReader::new(
            BufReader::new(Cursor::new(line.into_bytes())),
            TdtReaderOptions::default(),
        )
        .collect();
        assert_eq!(results.len(), 1);
        assert!(results[0].is_ok());
        assert_eq!(results[0].as_ref().unwrap().properties.len(), 5000);
    }

    #[test]
    fn very_large_molecule_smiles_does_not_panic() {
        let smi = "C".repeat(3000);
        let line = format!("$SMI<{smi}>\n|\n");
        let results: Vec<_> = TdtRecordReader::new(
            BufReader::new(Cursor::new(line.into_bytes())),
            TdtReaderOptions::default(),
        )
        .collect();
        assert_eq!(results.len(), 1);
        assert!(results[0].is_ok());
        assert_eq!(results[0].as_ref().unwrap().mol.atom_count(), 3000);
    }

    #[test]
    fn malformed_coordinate_list_does_not_panic() {
        let opts = TdtReaderOptions {
            read_2d: true,
            ..Default::default()
        };
        let input = b"$SMI<CCO>\n2D<not,numbers,here;>\n|\n";
        let results: Vec<_> =
            TdtRecordReader::new(BufReader::new(Cursor::new(input.to_vec())), opts).collect();
        assert_eq!(results.len(), 1);
        assert!(matches!(
            results[0],
            Err(TdtError::MalformedCoordinateList { .. })
        ));
    }

    #[test]
    fn coordinate_list_never_terminated_does_not_hang() {
        let opts = TdtReaderOptions {
            read_2d: true,
            max_line_bytes: 1 << 16,
            ..Default::default()
        };
        // No ";>" anywhere -- must not loop forever waiting for one.
        let input = "$SMI<CC>\n2D<1.0,2.0,3.0\n".repeat(50) + "|\n";
        let results: Vec<_> =
            TdtRecordReader::new(BufReader::new(Cursor::new(input.into_bytes())), opts).collect();
        // Whatever the outcome, this line finishing at all (not hanging) is the test.
        assert!(results.len() <= 1);
    }

    /// Tiny deterministic PRNG (splitmix64), matching
    /// `crate::smiles_table::adversarial_tests`' own approach -- no `rand`
    /// dependency needed for a seeded, reproducible mutation corpus.
    struct SplitMix64(u64);
    impl SplitMix64 {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        }
    }

    #[test]
    fn seeded_random_mutation_corpus_never_panics() {
        let base: &[u8] = b"$SMI<CCO>\nNAME<ethanol>\n2D<0.0,0.0,1.0,0.0,2.0,1.0;>\nNOTE<x>\n|\n";
        let mut rng = SplitMix64(0xC0FFEE1234567890);

        for _ in 0..2000 {
            let mut mutated = base.to_vec();
            let n_mutations = 1 + (rng.next() % 5) as usize;
            for _ in 0..n_mutations {
                if mutated.is_empty() {
                    break;
                }
                let idx = (rng.next() as usize) % mutated.len();
                let op = rng.next() % 3;
                match op {
                    0 => mutated[idx] = (rng.next() % 256) as u8,
                    1 => {
                        mutated.insert(idx, (rng.next() % 256) as u8);
                    }
                    _ => {
                        mutated.remove(idx);
                    }
                }
            }

            let opts = TdtReaderOptions {
                read_2d: true,
                max_line_bytes: 1 << 16,
                ..Default::default()
            };
            let _ = drain(&mutated, opts);
        }
    }
}
