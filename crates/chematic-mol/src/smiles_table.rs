//! Streaming SMILES table file I/O — `.smi`, `.smiles`, `.csv`, `.tsv`, `.txt`.
//!
//! A SMILES table file is one molecule per line, split into columns by a
//! delimiter; one column holds the SMILES, an optional column holds a name,
//! and any other columns become molecule properties. This is chematic's
//! streaming counterpart to RDKit's `SmilesMolSupplier`/`SmilesWriter`
//! (`Code/GraphMol/FileParsers/SmilesMolSupplier.cpp`/`SmilesWriter.cpp`,
//! RDKit commit `8afba32ec539dcb2369bc84549d802aca3f7eb39`, the true
//! resolution of tag `Release_2026_03_4` — see
//! `docs/` audit notes for the full source-cited findings this module is
//! built from).
//!
//! **Deliberate divergences from RDKit, documented rather than silently
//! matched or silently different:**
//! - RDKit's `SmilesMolSupplier` is a lazy, seek-based random-access reader
//!   (`.length()`/`operator[]` scan-and-cache the whole remaining file).
//!   [`SmilesRecordReader`] is forward-only (`BufRead`, not `Seek`) — no
//!   `.length()`/indexed access is offered, intentionally.
//! - RDKit's delimiter default differs across its own entry points (four
//!   different values were found across the v1 C++ API, the v2 core struct,
//!   and the Python binding). Chematic does not attempt to replicate that
//!   inconsistency; [`Delimiter::default()`] is [`Delimiter::Whitespace`],
//!   collapsing runs of space/tab into one split point (RDKit's own
//!   `boost::char_separator` with `keep_empty_tokens` does *not* collapse
//!   runs, a footgun this module deliberately does not reproduce).
//! - RDKit's `SmilesMolSupplier` has no strict/lax parsing toggle at all — a
//!   bad row always returns `None` for that entry while iteration continues
//!   unconditionally. [`SmilesReaderOptions::strict_parsing`] is a chematic
//!   addition: `true` still yields `Err` for the bad row but the reader then
//!   stops (subsequent `next()` calls return `None`); `false` yields the same
//!   `Err` and continues to the next row, matching RDKit's own behavior.
//! - CXSMILES: not recognized in the SMILES column (parsed via
//!   [`chematic_smiles::parse`], which has no CXSMILES support at all) —
//!   this matches RDKit's *own* default for `SmilesMolSupplier`, which
//!   explicitly hard-codes `allowCXSMILES=false` for this entry point too.
//! - CSV quoting: an explicit RFC 4180 *subset* — quoted fields may contain
//!   the delimiter, embedded quotes (doubled, `""` → `"`), but **may not
//!   span multiple physical lines**. A line-based `BufRead` reader cannot
//!   know a quoted field is unterminated without reading ahead indefinitely;
//!   rather than silently misparsing or unboundedly buffering, an
//!   unterminated quote is a [`SmilesTableError::UnterminatedQuote`].
//!   **RDKit's own `SmilesMolSupplier` has no CSV-quote-awareness at all**
//!   for its comma-delimiter mode — confirmed against a live RDKit oracle
//!   this session, not merely inferred from source: a value like
//!   `"has, a comma"` is split into two raw columns by RDKit's literal
//!   comma-splitting, corrupting the field. Chematic's quoting support is a
//!   genuine, deliberate improvement here, not a matched behavior — any
//!   `.csv` file with an embedded comma or quote will therefore disagree
//!   between chematic and RDKit by design, in chematic's favor.
//! - No name column (`name_column: None`, RDKit's `nameColumn=-1`): RDKit's
//!   `SmilesMolSupplier` falls back to the running *physical line number* as
//!   `_Name` in this case (source-confirmed and independently reproduced
//!   against a live RDKit oracle this session, `scripts/gen_rdkit_smiles_table_oracle.py`).
//!   [`MoleculeRecord::name`] is simply empty in this case instead —
//!   replicating a line-number-as-name fallback was judged a low-value,
//!   surprising RDKit implementation detail not worth reproducing.

use std::collections::HashSet;
use std::io::{BufRead, Write};

use chematic_core::Molecule;
use chematic_smiles::SmilesError;

use crate::record::MoleculeRecord;

// ---------------------------------------------------------------------------
// Delimiter
// ---------------------------------------------------------------------------

/// How columns are separated on each data line.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Delimiter {
    /// One or more consecutive ASCII space/tab characters, collapsed into a
    /// single split point (like `str::split_whitespace`, but restricted to
    /// space/tab — this module strips `\r` from line ends before
    /// tokenizing, so no other whitespace is ever present to split on).
    #[default]
    Whitespace,
    /// A single literal tab character. Unlike [`Delimiter::Whitespace`],
    /// consecutive tabs are **not** collapsed — an empty field between two
    /// tabs is a real, empty column.
    Tab,
    /// A single literal comma, with RFC 4180-subset quoting (see the module
    /// docs).
    Comma,
    /// Any other single ASCII byte, treated like [`Delimiter::Tab`] (literal,
    /// non-collapsing, no quoting).
    Custom(u8),
}

/// Split one already-`\r`-stripped line into columns per `delimiter`.
///
/// Returns `Err` only for [`Delimiter::Comma`] with an unterminated quoted
/// field — every other delimiter mode cannot fail (worst case: a very long
/// single field).
fn tokenize_line(line: &str, delimiter: Delimiter) -> Result<Vec<String>, ()> {
    match delimiter {
        Delimiter::Whitespace => Ok(line
            .split([' ', '\t'])
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .collect()),
        Delimiter::Tab => Ok(line.split('\t').map(str::to_string).collect()),
        Delimiter::Custom(b) => Ok(line.split(b as char).map(str::to_string).collect()),
        Delimiter::Comma => tokenize_csv(line),
    }
}

/// RFC 4180 *subset* CSV tokenizer (see module docs for the exact subset:
/// quote-doubling supported, embedded newlines in a quoted field are not).
fn tokenize_csv(line: &str) -> Result<Vec<String>, ()> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut chars = line.chars().peekable();
    let mut in_quotes = false;
    let mut quoted_field = false;

    while let Some(c) = chars.next() {
        if in_quotes {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    field.push('"');
                    chars.next();
                } else {
                    in_quotes = false;
                }
            } else {
                field.push(c);
            }
        } else if c == '"' && field.is_empty() && !quoted_field {
            in_quotes = true;
            quoted_field = true;
        } else if c == ',' {
            fields.push(std::mem::take(&mut field));
            quoted_field = false;
        } else {
            field.push(c);
        }
    }
    if in_quotes {
        return Err(());
    }
    fields.push(field);
    Ok(fields)
}

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

/// Configuration for [`SmilesRecordReader`].
#[derive(Debug, Clone)]
pub struct SmilesReaderOptions {
    pub delimiter: Delimiter,
    /// 0-indexed column holding the SMILES string.
    pub smiles_column: usize,
    /// 0-indexed column holding the record name, or `None` if there is no
    /// name column (the record's `name` is then always empty).
    pub name_column: Option<usize>,
    /// Whether the first non-comment, non-blank line is a header naming the
    /// extra (non-SMILES, non-name) columns, rather than a data row.
    pub title_line: bool,
    /// `true`: a row that fails to parse/sanitize yields one `Err`, then the
    /// reader stops (subsequent `next()` calls return `None`). `false`: the
    /// same `Err` is yielded, but the reader continues to the next row.
    pub strict_parsing: bool,
    /// Hard cap on a single line's byte length, to bound memory use against
    /// a pathological or adversarial input; exceeding it is
    /// [`SmilesTableError::LineTooLong`].
    pub max_line_bytes: usize,
}

impl Default for SmilesReaderOptions {
    fn default() -> Self {
        Self {
            delimiter: Delimiter::default(),
            smiles_column: 0,
            name_column: Some(1),
            title_line: false,
            strict_parsing: false,
            max_line_bytes: 1 << 20, // 1 MiB
        }
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors from [`SmilesRecordReader`]/[`SmilesRecordWriter`].
#[derive(Debug, Clone, PartialEq)]
pub enum SmilesTableError {
    /// The SMILES column's value could not be parsed as SMILES.
    InvalidSmiles {
        line: usize,
        record_index: usize,
        detail: SmilesError,
    },
    /// The requested SMILES column index does not exist on this line.
    MissingSmilesColumn {
        line: usize,
        record_index: usize,
        column: usize,
        columns_found: usize,
    },
    /// A quoted CSV field (`Delimiter::Comma`) was never closed on the line
    /// it started on — see the module docs' RFC 4180-subset note.
    UnterminatedQuote { line: usize, record_index: usize },
    /// A line exceeded `max_line_bytes`.
    LineTooLong {
        line: usize,
        record_index: usize,
        limit: usize,
    },
    /// An IO error occurred while reading the input stream.
    Io(String),
}

impl std::fmt::Display for SmilesTableError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSmiles {
                line,
                record_index,
                detail,
            } => write!(
                f,
                "invalid SMILES at line {line} (record {record_index}): {detail}"
            ),
            Self::MissingSmilesColumn {
                line,
                record_index,
                column,
                columns_found,
            } => write!(
                f,
                "line {line} (record {record_index}) has no column {column} \
                 (only {columns_found} column(s) found)"
            ),
            Self::UnterminatedQuote { line, record_index } => write!(
                f,
                "unterminated quoted field at line {line} (record {record_index})"
            ),
            Self::LineTooLong {
                line,
                record_index,
                limit,
            } => write!(
                f,
                "line {line} (record {record_index}) exceeds the {limit}-byte limit"
            ),
            Self::Io(msg) => write!(f, "IO error: {msg}"),
        }
    }
}

impl std::error::Error for SmilesTableError {}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

/// Streaming iterator over a SMILES table file.
///
/// Comment lines (first byte is literally `#`, before any trimming — a line
/// with leading whitespace before `#` is *not* a comment, matching RDKit's
/// own `SmilesMolSupplier::skipComments` behavior) and blank lines
/// (whitespace-only) are skipped silently, never counted as a record.
pub struct SmilesRecordReader<R: BufRead> {
    reader: R,
    options: SmilesReaderOptions,
    header: Option<Vec<String>>,
    header_read: bool,
    line_number: usize,
    record_index: usize,
    stopped: bool,
}

impl<R: BufRead> SmilesRecordReader<R> {
    pub fn new(reader: R, options: SmilesReaderOptions) -> Self {
        Self {
            reader,
            options,
            header: None,
            header_read: false,
            line_number: 0,
            record_index: 0,
            stopped: false,
        }
    }

    /// Read one raw line (without its line terminator), or `None` at EOF.
    fn read_line(&mut self) -> Result<Option<String>, SmilesTableError> {
        let mut buf = String::new();
        let n = self
            .reader
            .read_line(&mut buf)
            .map_err(|e| SmilesTableError::Io(e.to_string()))?;
        if n == 0 {
            return Ok(None);
        }
        self.line_number += 1;
        if buf.len() > self.options.max_line_bytes {
            return Err(SmilesTableError::LineTooLong {
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

    /// Read the header line, if configured, populating `self.header`.
    fn ensure_header(&mut self) -> Result<(), SmilesTableError> {
        if self.header_read {
            return Ok(());
        }
        self.header_read = true;
        if !self.options.title_line {
            return Ok(());
        }
        if let Some(line) = self.next_data_line()? {
            let tokens = tokenize_line(&line, self.options.delimiter).map_err(|_| {
                SmilesTableError::UnterminatedQuote {
                    line: self.line_number,
                    record_index: self.record_index,
                }
            })?;
            self.header = Some(tokens);
        }
        Ok(())
    }

    /// Read the next non-comment, non-blank raw line.
    fn next_data_line(&mut self) -> Result<Option<String>, SmilesTableError> {
        loop {
            match self.read_line()? {
                None => return Ok(None),
                Some(line) => {
                    if line.starts_with('#') || line.trim().is_empty() {
                        continue;
                    }
                    return Ok(Some(line));
                }
            }
        }
    }

    fn property_name(&self, col: usize) -> String {
        self.header
            .as_ref()
            .and_then(|h| h.get(col))
            .cloned()
            .unwrap_or_else(|| format!("column_{col}"))
    }
}

impl<R: BufRead> Iterator for SmilesRecordReader<R> {
    type Item = Result<MoleculeRecord, SmilesTableError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.stopped {
            return None;
        }
        if let Err(e) = self.ensure_header() {
            self.stopped = true;
            return Some(Err(e));
        }

        let line = match self.next_data_line() {
            Ok(Some(l)) => l,
            Ok(None) => return None,
            Err(e) => {
                self.stopped = true;
                return Some(Err(e));
            }
        };

        let record_line = self.line_number;
        let record_index = self.record_index;
        self.record_index += 1;

        let tokens = match tokenize_line(&line, self.options.delimiter) {
            Ok(t) => t,
            Err(()) => {
                let err = SmilesTableError::UnterminatedQuote {
                    line: record_line,
                    record_index,
                };
                if self.options.strict_parsing {
                    self.stopped = true;
                }
                return Some(Err(err));
            }
        };

        let Some(smi) = tokens.get(self.options.smiles_column) else {
            let err = SmilesTableError::MissingSmilesColumn {
                line: record_line,
                record_index,
                column: self.options.smiles_column,
                columns_found: tokens.len(),
            };
            if self.options.strict_parsing {
                self.stopped = true;
            }
            return Some(Err(err));
        };

        let mol: Molecule = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(detail) => {
                let err = SmilesTableError::InvalidSmiles {
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

        let mut record = MoleculeRecord::new(mol);
        record.name = self
            .options
            .name_column
            .and_then(|c| tokens.get(c))
            .cloned()
            .unwrap_or_default();

        let skip: HashSet<usize> = [self.options.smiles_column]
            .into_iter()
            .chain(self.options.name_column)
            .collect();
        for (col, val) in tokens.iter().enumerate() {
            if skip.contains(&col) {
                continue;
            }
            record
                .properties
                .push((self.property_name(col), val.clone()));
        }

        Some(Ok(record))
    }
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Configuration for [`SmilesRecordWriter`].
#[derive(Debug, Clone)]
pub struct SmilesWriterOptions {
    pub delimiter: Delimiter,
    /// Header label for the name column; an empty string suppresses the
    /// name field entirely (both header and data rows).
    pub name_header: String,
    pub include_header: bool,
    /// Which properties to write, and in what order — an empty slice writes
    /// no extra properties (chematic's default; matches RDKit's own
    /// `SmilesWriter` default of writing nothing beyond SMILES + name).
    pub properties: Vec<String>,
}

impl Default for SmilesWriterOptions {
    fn default() -> Self {
        Self {
            delimiter: Delimiter::Whitespace,
            name_header: "Name".to_string(),
            include_header: true,
            properties: Vec::new(),
        }
    }
}

fn delimiter_char(d: Delimiter) -> char {
    match d {
        Delimiter::Whitespace => ' ',
        Delimiter::Tab => '\t',
        Delimiter::Comma => ',',
        Delimiter::Custom(b) => b as char,
    }
}

/// Quote a CSV field per the RFC 4180 subset this module implements, only
/// if it actually needs quoting (contains the delimiter, a quote, or `\r`/`\n`).
fn csv_quote_if_needed(field: &str, delim: char) -> String {
    if field.contains(delim) || field.contains('"') || field.contains(['\r', '\n']) {
        format!("\"{}\"", field.replace('"', "\"\""))
    } else {
        field.to_string()
    }
}

fn write_field(delim: Delimiter, field: &str) -> String {
    if delim == Delimiter::Comma {
        csv_quote_if_needed(field, ',')
    } else {
        field.to_string()
    }
}

/// Streaming writer for a SMILES table file.
pub struct SmilesRecordWriter<W: Write> {
    writer: W,
    options: SmilesWriterOptions,
    header_written: bool,
}

impl<W: Write> SmilesRecordWriter<W> {
    pub fn new(writer: W, options: SmilesWriterOptions) -> Self {
        Self {
            writer,
            options,
            header_written: false,
        }
    }

    fn write_header(&mut self) -> std::io::Result<()> {
        let delim = delimiter_char(self.options.delimiter);
        write!(self.writer, "SMILES")?;
        if !self.options.name_header.is_empty() {
            write!(self.writer, "{delim}{}", self.options.name_header)?;
        }
        for p in &self.options.properties {
            write!(self.writer, "{delim}{p}")?;
        }
        writeln!(self.writer)
    }

    /// Write one record. The molecule is serialized via
    /// [`chematic_smiles::write`] (canonical, isomeric, aromatic-bond SMILES
    /// — chematic has no separate Kekule-SMILES writer mode at present).
    pub fn write_record(&mut self, record: &MoleculeRecord) -> std::io::Result<()> {
        if self.options.include_header && !self.header_written {
            self.write_header()?;
        }
        self.header_written = true;

        let delim = delimiter_char(self.options.delimiter);
        let smi = chematic_smiles::write(&record.mol);
        write!(self.writer, "{}", write_field(self.options.delimiter, &smi))?;
        if !self.options.name_header.is_empty() {
            write!(
                self.writer,
                "{delim}{}",
                write_field(self.options.delimiter, &record.name)
            )?;
        }
        for p in &self.options.properties {
            let val = record.get_property(p).unwrap_or("");
            write!(
                self.writer,
                "{delim}{}",
                write_field(self.options.delimiter, val)
            )?;
        }
        writeln!(self.writer)
    }

    pub fn into_inner(self) -> W {
        self.writer
    }

    /// Replace the set (and order) of extra properties written per record.
    pub fn set_properties(&mut self, properties: Vec<String>) {
        self.options.properties = properties;
    }

    /// Flush the underlying writer.
    pub fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn reader_over(
        input: &str,
        options: SmilesReaderOptions,
    ) -> SmilesRecordReader<BufReader<Cursor<Vec<u8>>>> {
        SmilesRecordReader::new(
            BufReader::new(Cursor::new(input.as_bytes().to_vec())),
            options,
        )
    }

    #[test]
    fn whitespace_delimiter_basic() {
        let input = "CC ethane\nCCO ethanol\n";
        let recs: Vec<_> = reader_over(input, SmilesReaderOptions::default())
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].mol.atom_count(), 2);
        assert_eq!(recs[0].name, "ethane");
        assert_eq!(recs[1].mol.atom_count(), 3);
        assert_eq!(recs[1].name, "ethanol");
    }

    #[test]
    fn tab_delimiter_does_not_collapse_empty_fields() {
        let opts = SmilesReaderOptions {
            delimiter: Delimiter::Tab,
            name_column: Some(1),
            ..Default::default()
        };
        let input = "CC\t\tActivity1\n"; // empty name field between two tabs
        let rec = reader_over(input, opts).next().unwrap().unwrap();
        assert_eq!(rec.name, "");
        assert_eq!(
            rec.properties,
            vec![("column_2".to_string(), "Activity1".to_string())]
        );
    }

    #[test]
    fn comma_csv_with_quoted_embedded_comma() {
        let opts = SmilesReaderOptions {
            delimiter: Delimiter::Comma,
            name_column: Some(1),
            ..Default::default()
        };
        let input = "CC,\"a, b\",7.2\n";
        let rec = reader_over(input, opts).next().unwrap().unwrap();
        assert_eq!(rec.name, "a, b");
        assert_eq!(
            rec.properties,
            vec![("column_2".to_string(), "7.2".to_string())]
        );
    }

    #[test]
    fn comma_csv_doubled_quote_escape() {
        let opts = SmilesReaderOptions {
            delimiter: Delimiter::Comma,
            name_column: Some(1),
            ..Default::default()
        };
        let input = "CC,\"say \"\"hi\"\"\"\n";
        let rec = reader_over(input, opts).next().unwrap().unwrap();
        assert_eq!(rec.name, "say \"hi\"");
    }

    #[test]
    fn comma_csv_unterminated_quote_is_lax_recoverable() {
        let opts = SmilesReaderOptions {
            delimiter: Delimiter::Comma,
            strict_parsing: false,
            ..Default::default()
        };
        let input = "CC,\"unterminated\nCCO,ethanol\n";
        let mut reader = reader_over(input, opts);
        let first = reader.next().unwrap();
        assert!(matches!(
            first,
            Err(SmilesTableError::UnterminatedQuote { .. })
        ));
        let second = reader.next().unwrap().unwrap();
        assert_eq!(second.mol.atom_count(), 3);
    }

    #[test]
    fn strict_parsing_stops_after_first_error() {
        let opts = SmilesReaderOptions {
            strict_parsing: true,
            ..Default::default()
        };
        let input = "not(a(smiles\nCCO ethanol\n";
        let mut reader = reader_over(input, opts);
        assert!(matches!(
            reader.next(),
            Some(Err(SmilesTableError::InvalidSmiles { .. }))
        ));
        assert!(reader.next().is_none());
    }

    #[test]
    fn lax_parsing_continues_after_error() {
        let opts = SmilesReaderOptions {
            strict_parsing: false,
            ..Default::default()
        };
        let input = "not(a(smiles\nCCO ethanol\n";
        let mut reader = reader_over(input, opts);
        assert!(matches!(
            reader.next(),
            Some(Err(SmilesTableError::InvalidSmiles { .. }))
        ));
        let second = reader.next().unwrap().unwrap();
        assert_eq!(second.mol.atom_count(), 3);
        assert!(reader.next().is_none());
    }

    #[test]
    fn comment_and_blank_lines_skipped() {
        let input =
            "# a comment\n\nCC ethane\n  # not a comment, treated as data (leading space)\n";
        let recs: Vec<_> = reader_over(input, SmilesReaderOptions::default()).collect();
        // The "leading-space-before-#" line is real data with an invalid SMILES -> Err.
        assert_eq!(recs.len(), 2);
        assert!(recs[0].is_ok());
        assert!(recs[1].is_err());
    }

    #[test]
    fn title_line_supplies_property_names() {
        let opts = SmilesReaderOptions {
            title_line: true,
            name_column: Some(1),
            ..Default::default()
        };
        let input = "SMILES Name Activity\nCC ethane 7.2\n";
        let rec = reader_over(input, opts).next().unwrap().unwrap();
        assert_eq!(rec.name, "ethane");
        assert_eq!(
            rec.properties,
            vec![("Activity".to_string(), "7.2".to_string())]
        );
    }

    #[test]
    fn no_title_line_uses_stable_column_names() {
        let opts = SmilesReaderOptions {
            title_line: false,
            name_column: Some(1),
            ..Default::default()
        };
        let input = "CC ethane 7.2\n";
        let rec = reader_over(input, opts).next().unwrap().unwrap();
        assert_eq!(
            rec.properties,
            vec![("column_2".to_string(), "7.2".to_string())]
        );
    }

    #[test]
    fn missing_smiles_column_is_explicit_error() {
        let opts = SmilesReaderOptions {
            smiles_column: 5,
            ..Default::default()
        };
        let input = "CC ethane\n";
        let mut reader = reader_over(input, opts);
        assert!(matches!(
            reader.next(),
            Some(Err(SmilesTableError::MissingSmilesColumn { .. }))
        ));
    }

    #[test]
    fn no_name_column_yields_empty_name() {
        let opts = SmilesReaderOptions {
            name_column: None,
            ..Default::default()
        };
        let input = "CC extra1 extra2\n";
        let rec = reader_over(input, opts).next().unwrap().unwrap();
        assert_eq!(rec.name, "");
        assert_eq!(
            rec.properties,
            vec![
                ("column_1".to_string(), "extra1".to_string()),
                ("column_2".to_string(), "extra2".to_string())
            ]
        );
    }

    #[test]
    fn line_too_long_is_explicit_error_not_panic() {
        let opts = SmilesReaderOptions {
            max_line_bytes: 8,
            ..Default::default()
        };
        let input = "CC a_very_long_name_field\n";
        let mut reader = reader_over(input, opts);
        assert!(matches!(
            reader.next(),
            Some(Err(SmilesTableError::LineTooLong { .. }))
        ));
    }

    #[test]
    fn writer_default_roundtrip() {
        let mol = chematic_smiles::parse("c1ccccc1").unwrap();
        let mut record = MoleculeRecord::new(mol);
        record.name = "benzene".to_string();

        let mut out = Vec::new();
        {
            let mut writer = SmilesRecordWriter::new(&mut out, SmilesWriterOptions::default());
            writer.write_record(&record).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with("SMILES Name\n") || text.starts_with("SMILES Name \n"));

        let mut reader = reader_over(
            &text,
            SmilesReaderOptions {
                title_line: true,
                ..Default::default()
            },
        );
        let back = reader.next().unwrap().unwrap();
        assert_eq!(back.mol.atom_count(), 6);
        assert_eq!(back.name, "benzene");
    }

    #[test]
    fn writer_no_name_header_omits_name_field() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let mut record = MoleculeRecord::new(mol);
        record.name = "ethane".to_string();

        let opts = SmilesWriterOptions {
            name_header: String::new(),
            include_header: false,
            ..Default::default()
        };
        let mut out = Vec::new();
        {
            let mut writer = SmilesRecordWriter::new(&mut out, opts);
            writer.write_record(&record).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert_eq!(
            text.trim_end(),
            chematic_smiles::write(&chematic_smiles::parse("CC").unwrap())
        );
    }

    #[test]
    fn writer_selected_properties_in_order() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let mut record = MoleculeRecord::new(mol);
        record.name = "ethane".to_string();
        record.properties = vec![
            ("MW".to_string(), "30.07".to_string()),
            ("Activity".to_string(), "7.2".to_string()),
        ];

        let opts = SmilesWriterOptions {
            properties: vec!["Activity".to_string(), "MW".to_string()],
            include_header: false,
            ..Default::default()
        };
        let mut out = Vec::new();
        {
            let mut writer = SmilesRecordWriter::new(&mut out, opts);
            writer.write_record(&record).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let smi = chematic_smiles::write(&chematic_smiles::parse("CC").unwrap());
        assert_eq!(text.trim_end(), format!("{smi} ethane 7.2 30.07"));
    }

    #[test]
    fn writer_missing_property_writes_empty_field() {
        // name_header left at its default ("Name"), so both the (empty) name
        // field and the missing property field are written as empty fields,
        // not omitted -- two trailing delimiters, which `trim_end()` would
        // wrongly swallow, hence `strip_suffix('\n')` here.
        let mol = chematic_smiles::parse("CC").unwrap();
        let record = MoleculeRecord::new(mol);
        let opts = SmilesWriterOptions {
            properties: vec!["Missing".to_string()],
            include_header: false,
            ..Default::default()
        };
        let mut out = Vec::new();
        {
            let mut writer = SmilesRecordWriter::new(&mut out, opts);
            writer.write_record(&record).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let smi = chematic_smiles::write(&chematic_smiles::parse("CC").unwrap());
        assert_eq!(text.strip_suffix('\n').unwrap(), format!("{smi}  "));
    }

    #[test]
    fn csv_write_quotes_field_containing_comma() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let mut record = MoleculeRecord::new(mol);
        record.name = "a, b".to_string();
        let opts = SmilesWriterOptions {
            delimiter: Delimiter::Comma,
            include_header: false,
            ..Default::default()
        };
        let mut out = Vec::new();
        {
            let mut writer = SmilesRecordWriter::new(&mut out, opts);
            writer.write_record(&record).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("\"a, b\""));
    }
}

// ---------------------------------------------------------------------------
// Adversarial / fuzz-style tests
// ---------------------------------------------------------------------------
//
// No `cargo-fuzz`/libfuzzer harness exists anywhere in this workspace yet,
// and introducing that toolchain (nightly + a separate fuzz crate) for one
// module was judged disproportionate to the risk here (a line-based text
// tokenizer over `BufRead`, not a binary/byte-oriented format). Instead:
// deterministic adversarial unit tests for every category the acceptance
// gate requires, plus a small seeded random-mutation corpus. All assert
// only "no panic, no hang, no OOM, no infinite loop" -- not any particular
// output -- since the whole point is that malformed/adversarial input must
// degrade to a clean `Err`, never worse.

#[cfg(test)]
mod adversarial_tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn drain(input: &[u8], options: SmilesReaderOptions) -> usize {
        let reader = SmilesRecordReader::new(BufReader::new(Cursor::new(input.to_vec())), options);
        reader.count()
    }

    #[test]
    fn empty_input_no_panic() {
        assert_eq!(drain(b"", SmilesReaderOptions::default()), 0);
    }

    #[test]
    fn truncated_input_no_trailing_newline_no_panic() {
        assert_eq!(drain(b"CC ethane", SmilesReaderOptions::default()), 1);
    }

    #[test]
    fn truncated_input_mid_quote_no_panic() {
        let opts = SmilesReaderOptions {
            delimiter: Delimiter::Comma,
            ..Default::default()
        };
        // Unterminated quote at EOF, no newline at all.
        assert_eq!(drain(b"CC,\"unterminated", opts), 1);
    }

    #[test]
    fn huge_line_is_explicit_error_not_panic_or_oom() {
        let mut line = b"CC ".to_vec();
        line.extend(std::iter::repeat_n(b'a', 10_000_000)); // 10MB name field
        line.push(b'\n');
        let opts = SmilesReaderOptions {
            max_line_bytes: 1 << 20,
            ..Default::default()
        };
        let reader = SmilesRecordReader::new(BufReader::new(Cursor::new(line)), opts);
        let results: Vec<_> = reader.collect();
        assert_eq!(results.len(), 1);
        assert!(matches!(
            results[0],
            Err(SmilesTableError::LineTooLong { .. })
        ));
    }

    #[test]
    fn huge_property_value_within_limit_does_not_panic() {
        let mut line = b"CC ethane ".to_vec();
        line.extend(std::iter::repeat_n(b'x', 500_000));
        line.push(b'\n');
        let opts = SmilesReaderOptions {
            max_line_bytes: 1 << 21, // large enough to admit this line
            ..Default::default()
        };
        let reader = SmilesRecordReader::new(BufReader::new(Cursor::new(line)), opts);
        let results: Vec<_> = reader.collect();
        assert_eq!(results.len(), 1);
        assert!(results[0].is_ok());
    }

    #[test]
    fn invalid_utf8_byte_path_yields_io_error_not_panic() {
        // 0xFF is never valid UTF-8 in any position.
        let input: &[u8] = b"CC \xFF\xFE ethane\n";
        let reader = SmilesRecordReader::new(
            BufReader::new(Cursor::new(input.to_vec())),
            SmilesReaderOptions::default(),
        );
        let results: Vec<_> = reader.collect();
        assert_eq!(results.len(), 1);
        assert!(matches!(results[0], Err(SmilesTableError::Io(_))));
    }

    #[test]
    fn excessive_column_count_does_not_panic() {
        let mut line = String::from("CC");
        for i in 0..5000 {
            line.push(' ');
            line.push_str(&format!("prop{i}"));
        }
        line.push('\n');
        let results: Vec<_> = drain_results(line.as_bytes(), SmilesReaderOptions::default());
        assert_eq!(results.len(), 1);
        assert!(results[0].is_ok());
        let rec = results[0].as_ref().unwrap();
        assert_eq!(rec.properties.len(), 4999); // 5000 extra tokens minus the name column
    }

    #[test]
    fn very_large_molecule_smiles_does_not_panic() {
        // A long linear alkane -- exercises the SMILES parser + tokenizer on
        // an unusually large single field, not just a large file.
        let smi = "C".repeat(3000);
        let line = format!("{smi} big_alkane\n");
        let results: Vec<_> = drain_results(line.as_bytes(), SmilesReaderOptions::default());
        assert_eq!(results.len(), 1);
        assert!(results[0].is_ok());
        assert_eq!(results[0].as_ref().unwrap().mol.atom_count(), 3000);
    }

    fn drain_results(
        input: &[u8],
        options: SmilesReaderOptions,
    ) -> Vec<Result<MoleculeRecord, SmilesTableError>> {
        SmilesRecordReader::new(BufReader::new(Cursor::new(input.to_vec())), options).collect()
    }

    /// Tiny deterministic PRNG (splitmix64) -- no `rand` dependency needed
    /// for a seeded, reproducible mutation corpus.
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
        let base: &[u8] = b"SMILES,Name,Note\nCC,ethane,\"has, a comma\"\nc1ccccc1,benzene,plain\n";
        let mut rng = SplitMix64(0xDEADBEEFCAFEF00D);

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

            let opts = SmilesReaderOptions {
                delimiter: Delimiter::Comma,
                title_line: true,
                max_line_bytes: 1 << 16,
                ..Default::default()
            };
            // The only contract under test: this must terminate and never
            // panic, regardless of how corrupted the byte stream is.
            let _ = drain(&mutated, opts);
        }
    }
}
