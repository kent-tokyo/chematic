//! SDF (Structure-Data File) reader.
//!
//! An SDF file contains one or more MOL V2000 blocks separated by `$$$$`
//! delimiter lines.  Data-field sections between `M  END` and `$$$$` are
//! accepted but ignored.

use chematic_core::{Coords3D, Molecule};
use chematic_perception::{EzDirectionDiagnostic, StereoDiagnostic};

use crate::error::MolParseError;
use crate::mol2000::{
    CoordinateDimension, GeometryRank, MolMetadata, Stereo3DDiagnostic, parse_mol,
    read_mol_with_diagnostics,
};

/// Iterator over molecules in an SDF string.
///
/// Each call to `next()` returns the next `(Molecule, MolMetadata)` pair
/// parsed from the string, or the first `MolParseError` encountered.
/// Returns `None` when the entire input has been consumed.
pub struct SdfReader<'a> {
    remaining: &'a str,
    current_mol_num: usize,
}

impl<'a> SdfReader<'a> {
    /// Create a new `SdfReader` over the given SDF string.
    pub fn new(input: &'a str) -> Self {
        Self {
            remaining: input,
            current_mol_num: 0,
        }
    }
}

impl<'a> Iterator for SdfReader<'a> {
    type Item = Result<(Molecule, MolMetadata), MolParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        // A blank first line is a legal MOL name line (issue #171) -- do not
        // skip it. A genuinely empty gap between/after `$$$$` delimiters is
        // already handled below via `mol_block.trim().is_empty()`.
        if self.remaining.is_empty() {
            return None;
        }

        self.current_mol_num += 1;

        // Scan line by line so that a `$$$$` substring inside a data value
        // does not trigger a false match.  When the delimiter is found, the
        // mol block runs up to (but excluding) it, and the rest continues
        // after the delimiter line.  When EOF is reached without a delimiter,
        // the entire remainder is treated as a single mol block.
        let mut byte_offset = 0usize;
        let (end_byte, after_delim) = loop {
            let rest = &self.remaining[byte_offset..];
            match rest.find('\n') {
                Some(nl) => {
                    let line = rest[..nl].trim_end_matches('\r');
                    if line == "$$$$" {
                        break (byte_offset, &self.remaining[byte_offset + nl + 1..]);
                    }
                    byte_offset += nl + 1;
                }
                None => {
                    // Last line, no trailing newline.
                    if rest.trim_end_matches('\r') == "$$$$" {
                        break (byte_offset, "");
                    }
                    break (self.remaining.len(), "");
                }
            }
        };

        let mol_block = &self.remaining[..end_byte];
        self.remaining = after_delim;

        if mol_block.trim().is_empty() {
            // Empty block between two `$$$$` lines — skip and try next.
            return self.next();
        }

        Some(parse_mol(mol_block))
    }
}

/// Parse all molecules from an SDF string.
///
/// Stops and returns an error on the first parse failure.
pub fn parse_sdf(input: &str) -> Result<Vec<(Molecule, MolMetadata)>, MolParseError> {
    parse_sdf_with_limits(input, SdfParseLimits::default())
}

/// Parse all molecules from an SDF string with explicit resource limits.
///
/// This uses the same bounded record path as [`SdfFileReader`], while keeping
/// the compact `(Molecule, MolMetadata)` result of [`parse_sdf`].
pub fn parse_sdf_with_limits(
    input: &str,
    limits: SdfParseLimits,
) -> Result<Vec<(Molecule, MolMetadata)>, MolParseError> {
    use std::io::{BufReader, Cursor};

    SdfFileReader::with_limits(BufReader::new(Cursor::new(input.as_bytes())), limits)
        .map(|result| result.map(|record| (record.mol, record.meta)))
        .collect()
}

// ---------------------------------------------------------------------------
// SdfRecord — molecule + SD data fields
// ---------------------------------------------------------------------------

/// A fully-parsed SDF record: molecule, metadata, 2D coordinates, and SD properties.
///
/// `SdfRecordReader` yields one `SdfRecord` per molecule in an SDF string.
/// The `properties` map corresponds to `> <FieldName>` data blocks.
pub struct SdfRecord {
    /// Parsed molecule (heavy atoms only, no explicit H).
    pub mol: Molecule,
    /// Metadata from the three-line MOL header (name, comment).
    pub meta: MolMetadata,
    /// 2D atom coordinates in Å, indexed by atom position.
    /// Empty when the MOL block contains no coordinate data.
    pub coords: Vec<(f64, f64)>,
    /// SD data fields.  Keys are field names; values are field content.
    /// Multi-line values are joined with `\n`.
    pub properties: std::collections::HashMap<String, String>,
    /// Rejected wedge/hash stereocenters from this record (see
    /// [`chematic_perception::StereoDiagnostic`]). Empty unless a wedge/hash
    /// bond was present at some center and got rejected.
    pub stereo_diagnostics: Vec<StereoDiagnostic>,
    /// Rejected stereogenic double bonds from this record (see
    /// [`chematic_perception::EzDirectionDiagnostic`]). Empty unless a
    /// stereogenic double bond's direction was rejected.
    pub ez_diagnostics: Vec<EzDirectionDiagnostic>,
    /// Real 3D coordinates, `Some` exactly when the record's atom block has
    /// non-(near-)zero z values -- see [`crate::mol2000::MolReadReport::conformer`].
    pub conformer: Option<Coords3D>,
    /// The record's own header-declared dimensionality (line 2's "2D"/"3D"
    /// tag). `Unknown` is the common case -- most writers never populate it.
    pub coordinate_dimension: CoordinateDimension,
    /// What the record's actual coordinates look like, independent of
    /// `coordinate_dimension` -- see [`crate::mol2000::GeometryRank`].
    pub geometry_rank: GeometryRank,
    /// 3D-geometry-related diagnostics for this record (dimension-vs-geometry
    /// mismatches, wedge-vs-3D-geometry parity conflicts). See
    /// [`crate::mol2000::Stereo3DDiagnostic`] -- empty does not mean
    /// "verified correct", it means nothing was declared to check.
    pub stereo3d_diagnostics: Vec<Stereo3DDiagnostic>,
}

/// Iterator over SDF records that also captures SD data fields.
///
/// Unlike [`SdfReader`], this iterator yields [`SdfRecord`] values so that
/// callers can access per-molecule properties (e.g. activity values, MW, etc.).
pub struct SdfRecordReader<'a> {
    remaining: &'a str,
    input_bytes: usize,
    records_read: usize,
    limits: SdfParseLimits,
}

impl<'a> SdfRecordReader<'a> {
    /// Create a new reader over the given SDF string.
    pub fn new(input: &'a str) -> Self {
        Self::with_limits(input, SdfParseLimits::default())
    }

    /// Create a record reader with explicit input, line, record-size, and
    /// record-count limits.
    pub fn with_limits(input: &'a str, limits: SdfParseLimits) -> Self {
        Self {
            remaining: input,
            input_bytes: input.len(),
            records_read: 0,
            limits,
        }
    }
}

impl<'a> Iterator for SdfRecordReader<'a> {
    type Item = Result<SdfRecord, MolParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.input_bytes > self.limits.max_input_bytes {
                self.remaining = "";
                return Some(Err(MolParseError::ResourceLimit {
                    resource: "input bytes",
                    actual: self.input_bytes,
                    limit: self.limits.max_input_bytes,
                }));
            }
            // A blank first line is a legal MOL name line (issue #171) -- do not
            // skip it. A genuinely empty gap between/after `$$$$` delimiters is
            // already handled below via `block.trim().is_empty()`.
            if self.remaining.is_empty() {
                return None;
            }

            // Scan to find the $$$$ delimiter (line-by-line to avoid false matches).
            let mut byte_offset = 0usize;
            let (end_byte, after_delim) = loop {
                let rest = &self.remaining[byte_offset..];
                match rest.find('\n') {
                    Some(nl) => {
                        let line = rest[..nl].trim_end_matches('\r');
                        if line.len() > self.limits.max_line_bytes {
                            self.remaining = "";
                            return Some(Err(MolParseError::ResourceLimit {
                                resource: "line bytes",
                                actual: line.len(),
                                limit: self.limits.max_line_bytes,
                            }));
                        }
                        if line == "$$$$" {
                            break (byte_offset, &self.remaining[byte_offset + nl + 1..]);
                        }
                        byte_offset += nl + 1;
                    }
                    None => {
                        if rest.trim_end_matches('\r') == "$$$$" {
                            break (byte_offset, "");
                        }
                        if rest.trim_end_matches('\r').len() > self.limits.max_line_bytes {
                            self.remaining = "";
                            return Some(Err(MolParseError::ResourceLimit {
                                resource: "line bytes",
                                actual: rest.trim_end_matches('\r').len(),
                                limit: self.limits.max_line_bytes,
                            }));
                        }
                        break (self.remaining.len(), "");
                    }
                }
            };

            let block = &self.remaining[..end_byte];
            self.remaining = after_delim;

            if block.trim().is_empty() {
                continue;
            }

            if block.len() > self.limits.max_record_bytes {
                self.remaining = "";
                return Some(Err(MolParseError::ResourceLimit {
                    resource: "record bytes",
                    actual: block.len(),
                    limit: self.limits.max_record_bytes,
                }));
            }
            if self.records_read >= self.limits.max_records {
                self.remaining = "";
                return Some(Err(MolParseError::ResourceLimit {
                    resource: "records",
                    actual: self.records_read.saturating_add(1),
                    limit: self.limits.max_records,
                }));
            }
            self.records_read += 1;

            // Parse molecule + 2D coordinates + stereo diagnostics.
            let report = match read_mol_with_diagnostics(block) {
                Ok(report) => report,
                Err(e) => return Some(Err(e)),
            };

            // Extract data fields from the part after "M  END".
            let data_part = block
                .find("M  END")
                .map(|pos| &block[pos + 6..])
                .unwrap_or("");
            let properties: std::collections::HashMap<String, String> =
                parse_sd_fields(data_part).into_iter().collect();

            return Some(Ok(SdfRecord {
                mol: report.mol,
                meta: report.metadata,
                coords: report.coords,
                properties,
                stereo_diagnostics: report.stereo_diagnostics,
                ez_diagnostics: report.ez_diagnostics,
                conformer: report.conformer,
                coordinate_dimension: report.coordinate_dimension,
                geometry_rank: report.geometry_rank,
                stereo3d_diagnostics: report.stereo3d_diagnostics,
            }));
        }
    }
}

/// Parse SD data fields from the section after `M  END`.
///
/// Each field starts with `> <FieldName>` on its own line.  The value is
/// everything on subsequent lines until a blank line (or end of input).
fn parse_sd_fields(data: &str) -> Vec<(String, String)> {
    let mut fields = Vec::new();
    let mut current_key: Option<String> = None;
    let mut current_value_lines: Vec<&str> = Vec::new();

    for raw_line in data.lines() {
        let line = raw_line.trim_end_matches('\r');

        if let Some(key) = parse_sd_field_header(line) {
            // Flush previous field.
            if let Some(k) = current_key.take() {
                fields.push((k, current_value_lines.join("\n")));
                current_value_lines.clear();
            }
            current_key = Some(key);
        } else if line.is_empty() {
            // Blank line ends the current field's value.
            if let Some(k) = current_key.take() {
                fields.push((k, current_value_lines.join("\n")));
                current_value_lines.clear();
            }
        } else if current_key.is_some() {
            current_value_lines.push(line);
        }
    }
    // Flush trailing field with no blank line.
    if let Some(k) = current_key {
        fields.push((k, current_value_lines.join("\n")));
    }

    fields
}

/// Parse `> <FieldName>` header lines, returning the field name or `None`.
fn parse_sd_field_header(line: &str) -> Option<String> {
    // SDF spec: field headers start with "> " and contain the name in `<...>`.
    // We accept "> <Name>" and also ">  <Name>" (extra spaces).
    let rest = line.strip_prefix('>')?;
    let rest = rest.trim();
    let inner = rest.strip_prefix('<')?.strip_suffix('>')?;
    Some(inner.to_string())
}

// ---------------------------------------------------------------------------
// Conformer ensembles — bundle records sharing the same molecular graph
// ---------------------------------------------------------------------------

/// One group of SDF records that share the same molecular graph (see
/// [`same_graph_identity`]) and each carry a real 3D conformer.
///
/// Records with no 3D conformer (`conformer.is_none()`, i.e. flagged 2D or
/// flat -- see [`crate::mol2000::MolReadReport::conformer`]) are not part of
/// any ensemble: this supplier is specifically about bundling repeated 3D
/// conformers of "the same" molecule (e.g. `EmbedMultipleConfs` +
/// `MolToMolBlock` per conformer, written as consecutive SDF records).
pub struct ConformerEnsemble {
    /// The molecular graph shared by every member (the first record's
    /// `mol`; every other member's `mol` is checked equal by
    /// [`same_graph_identity`] and then discarded, not stored again).
    pub mol: Molecule,
    /// Metadata from the first record in the group.
    pub metadata: MolMetadata,
    /// One conformer per member record, in file order.
    pub conformers: Vec<Coords3D>,
}

/// Deliberately simple identity check: same atom count/order, same
/// per-index `(element, charge, isotope)`, same bond list `(atom1, atom2,
/// order)` in index order.
///
/// This is order-sensitive structural equality, NOT graph isomorphism or
/// canonical-SMILES equality -- and that is intentional, not a shortcut
/// taken for lack of time: it is correct for the actual use case this
/// supplier targets (the same molecule written N times with different
/// coordinates by the same embedding tool, which never reorders atoms or
/// bonds between calls -- confirmed against RDKit's own
/// `EmbedMultipleConfs` + per-conformer `MolToMolBlock`, see this crate's
/// fixture tests). A more rigorous, reordering-tolerant identity layer is a
/// separate, parallel workstream in the 3D Breakthrough Program (Agent B);
/// this one is intentionally not that, and the Coordinator is expected to
/// reconcile the two later.
fn same_graph_identity(a: &Molecule, b: &Molecule) -> bool {
    if a.atom_count() != b.atom_count() || a.bond_count() != b.bond_count() {
        return false;
    }
    for ((_, atom_a), (_, atom_b)) in a.atoms().zip(b.atoms()) {
        if atom_a.element != atom_b.element
            || atom_a.charge != atom_b.charge
            || atom_a.isotope != atom_b.isotope
        {
            return false;
        }
    }
    for ((_, bond_a), (_, bond_b)) in a.bonds().zip(b.bonds()) {
        if bond_a.atom1 != bond_b.atom1
            || bond_a.atom2 != bond_b.atom2
            || bond_a.order != bond_b.order
        {
            return false;
        }
    }
    true
}

/// Parse an SDF string and group records that share the same molecular graph
/// (see [`same_graph_identity`]) and each carry a real 3D conformer into
/// [`ConformerEnsemble`]s.
///
/// Records with no 3D conformer are silently omitted from the result (this
/// is specifically a *3D conformer* ensemble supplier, not a general SDF
/// grouping utility). Grouping is order-independent within the file (a
/// record is compared against every existing group's representative, not
/// just its immediate predecessor), matching how real multi-conformer SDF
/// exports are usually -- but not always -- written consecutively.
///
/// Stops and returns an error on the first parse failure, same as
/// [`read_sdf_with_diagnostics`].
pub fn read_sdf_conformer_ensembles(input: &str) -> Result<Vec<ConformerEnsemble>, MolParseError> {
    let reports = crate::mol2000::read_sdf_with_diagnostics(input)?;
    let mut ensembles: Vec<ConformerEnsemble> = Vec::new();

    for report in reports {
        let Some(conformer) = report.conformer else {
            continue;
        };
        if let Some(existing) = ensembles
            .iter_mut()
            .find(|e| same_graph_identity(&e.mol, &report.mol))
        {
            existing.conformers.push(conformer);
        } else {
            ensembles.push(ConformerEnsemble {
                mol: report.mol,
                metadata: report.metadata,
                conformers: vec![conformer],
            });
        }
    }

    Ok(ensembles)
}

// ---------------------------------------------------------------------------
// SdfFileReader — streaming reader for file-backed SDF
// ---------------------------------------------------------------------------

/// Streaming SDF iterator over any [`std::io::BufRead`] source.
///
/// Unlike [`SdfRecordReader`] (which requires the entire file in memory as a
/// `&str`), `SdfFileReader` reads one MOL block at a time, suitable for large
/// SDF files.  IO errors are surfaced as [`MolParseError::Io`]; molecule parse
/// errors are returned as `Err` items so the caller can decide to skip or stop.
pub struct SdfFileReader<R: std::io::BufRead> {
    reader: R,
    done: bool,
    limits: SdfParseLimits,
    bytes_read: usize,
    records_read: usize,
}

/// Resource limits for streaming SDF input.
#[derive(Debug, Clone, Copy)]
pub struct SdfParseLimits {
    /// Maximum bytes read from the source.
    pub max_input_bytes: usize,
    /// Maximum bytes in one MOL/data record.
    pub max_record_bytes: usize,
    /// Maximum bytes in one physical input line.
    pub max_line_bytes: usize,
    /// Maximum non-empty records yielded.
    pub max_records: usize,
}

impl Default for SdfParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 1 << 30,
            max_record_bytes: 16 << 20,
            max_line_bytes: 16 << 20,
            max_records: 100_000,
        }
    }
}

impl<R: std::io::BufRead> SdfFileReader<R> {
    /// Wrap any `BufRead` source (e.g. `BufReader<File>`).
    pub fn new(reader: R) -> Self {
        Self::with_limits(reader, SdfParseLimits::default())
    }

    /// Wrap a `BufRead` source and enforce input, record-size, and record-count limits.
    pub fn with_limits(reader: R, limits: SdfParseLimits) -> Self {
        Self {
            reader,
            done: false,
            limits,
            bytes_read: 0,
            records_read: 0,
        }
    }
}

impl<R: std::io::BufRead> Iterator for SdfFileReader<R> {
    type Item = Result<SdfRecord, MolParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let mut block = String::with_capacity(2048);
        let mut line = String::new();

        loop {
            line.clear();
            match self.reader.read_line(&mut line) {
                Err(e) => {
                    self.done = true;
                    return Some(Err(MolParseError::Io(e.to_string())));
                }
                Ok(0) => {
                    // EOF
                    self.done = true;
                    if block.trim().is_empty() {
                        return None;
                    }
                    // Trailing block without $$$$ delimiter — parse it.
                    break;
                }
                Ok(_) => {
                    if line.len() > self.limits.max_line_bytes {
                        self.done = true;
                        return Some(Err(MolParseError::ResourceLimit {
                            resource: "line bytes",
                            actual: line.len(),
                            limit: self.limits.max_line_bytes,
                        }));
                    }
                    self.bytes_read = self.bytes_read.saturating_add(line.len());
                    if self.bytes_read > self.limits.max_input_bytes {
                        self.done = true;
                        return Some(Err(MolParseError::ResourceLimit {
                            resource: "input bytes",
                            actual: self.bytes_read,
                            limit: self.limits.max_input_bytes,
                        }));
                    }
                    let trimmed = line.trim_end_matches(['\r', '\n']);
                    if trimmed == "$$$$" {
                        break;
                    }
                    block.push_str(&line);
                }
            }
        }

        if block.trim().is_empty() {
            // Empty block (e.g. two consecutive $$$$) — advance to next.
            return self.next();
        }

        if block.len() > self.limits.max_record_bytes {
            self.done = true;
            return Some(Err(MolParseError::ResourceLimit {
                resource: "record bytes",
                actual: block.len(),
                limit: self.limits.max_record_bytes,
            }));
        }
        if self.records_read >= self.limits.max_records {
            self.done = true;
            return Some(Err(MolParseError::ResourceLimit {
                resource: "records",
                actual: self.records_read.saturating_add(1),
                limit: self.limits.max_records,
            }));
        }
        self.records_read += 1;

        // Reuse the same parse path as SdfRecordReader.
        let report = match read_mol_with_diagnostics(&block) {
            Ok(report) => report,
            Err(e) => return Some(Err(e)),
        };

        let data_part = block
            .find("M  END")
            .map(|pos| &block[pos + 6..])
            .unwrap_or("");
        let properties: std::collections::HashMap<String, String> =
            parse_sd_fields(data_part).into_iter().collect();

        Some(Ok(SdfRecord {
            mol: report.mol,
            meta: report.metadata,
            coords: report.coords,
            properties,
            stereo_diagnostics: report.stereo_diagnostics,
            ez_diagnostics: report.ez_diagnostics,
            conformer: report.conformer,
            coordinate_dimension: report.coordinate_dimension,
            geometry_rank: report.geometry_rank,
            stereo3d_diagnostics: report.stereo3d_diagnostics,
        }))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const MOL_A: &str = "\
mol_a
  chematic

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
";

    const MOL_B: &str = "\
mol_b
  chematic

  3  2  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
M  END
";

    fn two_mol_sdf() -> String {
        format!("{MOL_A}$$$$\n{MOL_B}$$$$\n")
    }

    #[test]
    fn test_sdf_reader_two_molecules() {
        let sdf = two_mol_sdf();
        let results: Vec<_> = SdfReader::new(&sdf).collect();
        assert_eq!(results.len(), 2);
        let (mol_a, meta_a) = results[0].as_ref().expect("mol_a parse");
        let (mol_b, meta_b) = results[1].as_ref().expect("mol_b parse");
        assert_eq!(mol_a.atom_count(), 2);
        assert_eq!(mol_a.bond_count(), 1);
        assert_eq!(meta_a.name, "mol_a");
        assert_eq!(mol_b.atom_count(), 3);
        assert_eq!(mol_b.bond_count(), 2);
        assert_eq!(meta_b.name, "mol_b");
    }

    #[test]
    fn test_parse_sdf_all() {
        let sdf = two_mol_sdf();
        let mols = parse_sdf(&sdf).expect("parse_sdf");
        assert_eq!(mols.len(), 2);
    }

    #[test]
    fn test_sdf_reader_single_molecule_no_delimiter() {
        // An SDF with a single molecule that has no trailing $$$$ is still valid.
        let results: Vec<_> = SdfReader::new(MOL_A).collect();
        assert_eq!(results.len(), 1);
        let (mol, _) = results[0].as_ref().expect("parse");
        assert_eq!(mol.atom_count(), 2);
    }

    #[test]
    fn test_sdf_reader_stops_on_error() {
        // Second molecule has a bad counts line; parse_sdf should return Err.
        let bad_sdf = format!("{MOL_A}$$$$\nbad\n  prog\n\n  X  Y\nM  END\n$$$$\n");
        let result = parse_sdf(&bad_sdf);
        assert!(result.is_err());
    }

    #[test]
    fn test_sdf_reader_empty_input() {
        let results: Vec<_> = SdfReader::new("").collect();
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_sdf_reader_names_preserved() {
        let sdf = two_mol_sdf();
        let mols = parse_sdf(&sdf).expect("parse");
        assert_eq!(mols[0].1.name, "mol_a");
        assert_eq!(mols[1].1.name, "mol_b");
    }

    #[test]
    fn test_sdf_with_data_fields() {
        // SDF with data fields between M  END and $$$$ — should be ignored.
        let sdf_with_data = format!("{MOL_A}> <MW>\n44.0\n\n$$$$\n");
        let results: Vec<_> = SdfReader::new(&sdf_with_data).collect();
        assert_eq!(results.len(), 1);
        let (mol, _) = results[0].as_ref().expect("parse");
        assert_eq!(mol.atom_count(), 2);
    }

    #[test]
    fn test_sdf_reader_reports_truncated_large_count_record() {
        let bad_sdf = "\
max_atoms
  chematic

999  0  0  0  0  0  0  0  0  0  0 V2000
$$$$
";
        let mut reader = SdfReader::new(bad_sdf);
        assert!(matches!(
            reader.next(),
            Some(Err(MolParseError::UnexpectedEnd))
        ));
        assert!(reader.next().is_none());
        assert!(matches!(
            parse_sdf(bad_sdf),
            Err(MolParseError::UnexpectedEnd)
        ));
    }

    #[test]
    fn test_sdf_file_reader_streaming() {
        // Use Cursor<Vec<u8>> as an in-memory BufRead to avoid touching the filesystem.
        use std::io::{BufReader, Cursor};
        let sdf = two_mol_sdf();
        let cursor = Cursor::new(sdf.into_bytes());
        let records: Vec<_> = SdfFileReader::new(BufReader::new(cursor))
            .filter_map(|r| r.ok())
            .collect();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].mol.atom_count(), 2); // mol_a: 2 C atoms
        assert_eq!(records[1].mol.atom_count(), 3); // mol_b: C, N, O
    }

    #[test]
    fn test_sdf_file_reader_skips_empty_block() {
        use std::io::{BufReader, Cursor};
        let sdf = format!("$$$$\n{MOL_A}$$$$\n");
        let cursor = Cursor::new(sdf.into_bytes());
        let records: Vec<_> = SdfFileReader::new(BufReader::new(cursor))
            .filter_map(|r| r.ok())
            .collect();
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn test_sdf_file_reader_enforces_resource_limits() {
        use std::io::{BufReader, Cursor};

        let limits = SdfParseLimits {
            max_record_bytes: 8,
            ..SdfParseLimits::default()
        };
        let result =
            SdfFileReader::with_limits(BufReader::new(Cursor::new(MOL_A.as_bytes())), limits)
                .next()
                .unwrap();
        assert!(matches!(
            result,
            Err(MolParseError::ResourceLimit {
                resource: "record bytes",
                ..
            })
        ));

        let limits = SdfParseLimits {
            max_records: 1,
            ..SdfParseLimits::default()
        };
        let input = format!("{MOL_A}$$$$\n{MOL_B}$$$$\n");
        let mut reader =
            SdfFileReader::with_limits(BufReader::new(Cursor::new(input.into_bytes())), limits);
        assert!(reader.next().unwrap().is_ok());
        assert!(matches!(
            reader.next().unwrap(),
            Err(MolParseError::ResourceLimit {
                resource: "records",
                ..
            })
        ));

        let limits = SdfParseLimits {
            max_line_bytes: 8,
            ..SdfParseLimits::default()
        };
        let result =
            SdfFileReader::with_limits(BufReader::new(Cursor::new(MOL_A.as_bytes())), limits)
                .next()
                .unwrap();
        assert!(matches!(
            result,
            Err(MolParseError::ResourceLimit {
                resource: "line bytes",
                ..
            })
        ));
    }

    #[test]
    fn test_sdf_record_reader_enforces_limits_without_recursive_skip() {
        let input = format!("$$$$\n{}$$$$\n", MOL_A);
        let mut reader = SdfRecordReader::with_limits(
            &input,
            SdfParseLimits {
                max_records: 0,
                ..SdfParseLimits::default()
            },
        );
        assert!(matches!(
            reader.next(),
            Some(Err(MolParseError::ResourceLimit {
                resource: "records",
                ..
            }))
        ));

        let mut reader = SdfRecordReader::with_limits(
            &input,
            SdfParseLimits {
                max_input_bytes: 8,
                ..SdfParseLimits::default()
            },
        );
        assert!(matches!(
            reader.next(),
            Some(Err(MolParseError::ResourceLimit {
                resource: "input bytes",
                ..
            }))
        ));
    }

    #[test]
    fn test_parse_sdf_with_limits_uses_bounded_record_path() {
        let input = two_mol_sdf();
        let result = parse_sdf_with_limits(
            &input,
            SdfParseLimits {
                max_records: 1,
                ..SdfParseLimits::default()
            },
        );
        assert!(matches!(
            result,
            Err(MolParseError::ResourceLimit {
                resource: "records",
                ..
            })
        ));

        let result = parse_sdf_with_limits(
            &input,
            SdfParseLimits {
                max_record_bytes: 8,
                ..SdfParseLimits::default()
            },
        );
        assert!(matches!(
            result,
            Err(MolParseError::ResourceLimit {
                resource: "record bytes",
                ..
            })
        ));
    }

    #[test]
    fn test_sdf_file_reader_malformed_record_yields_err() {
        use std::io::{BufReader, Cursor};

        // 3-record SDF: valid, malformed (bad atom count line), valid
        let malformed = "broken\n  prog\n\n  NOTNUM  0  0 V2000\nM  END\n";
        let sdf = format!("{MOL_A}$$$$\n{malformed}$$$$\n{MOL_B}$$$$\n");
        let cursor = Cursor::new(sdf.into_bytes());
        let results: Vec<_> = SdfFileReader::new(BufReader::new(cursor)).collect();

        // Three items: Ok, Err, Ok
        assert_eq!(results.len(), 3);
        assert!(results[0].is_ok(), "first record should be ok");
        assert!(
            results[1].is_err(),
            "second record should be err (malformed)"
        );
        assert!(results[2].is_ok(), "third record should be ok");
    }

    // ── Issue #171: blank MOL name line must not be eaten as inter-record
    // padding ──────────────────────────────────────────────────────────────
    // Fixture generated by a live RDKit 2026.03.3 oracle:
    //   AllChem.Compute2DCoords(Chem.MolFromSmiles("CC")); Chem.MolToMolBlock(mol)
    // RDKit's own MolToMolBlock leaves the name line blank for an unnamed
    // molecule -- this is the literal, spec-legal shape that previously broke
    // SdfReader/SdfRecordReader (see issue #171's repro).
    const MOL_BLANK_NAME: &str = "\
\n     RDKit          2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n";

    fn crlf(s: &str) -> String {
        s.replace('\n', "\r\n")
    }

    #[test]
    fn test_sdf_reader_blank_name_first_record() {
        let sdf = format!("{MOL_BLANK_NAME}$$$$\n{MOL_A}$$$$\n");
        let results: Vec<_> = SdfReader::new(&sdf).collect();
        assert_eq!(results.len(), 2);
        let (mol0, meta0) = results[0].as_ref().expect("blank-name record parse");
        assert_eq!(meta0.name, "");
        assert_eq!(mol0.atom_count(), 2);
        let (mol1, meta1) = results[1].as_ref().expect("mol_a parse");
        assert_eq!(meta1.name, "mol_a");
        assert_eq!(mol1.atom_count(), 2);
    }

    #[test]
    fn test_sdf_reader_blank_name_middle_record() {
        let sdf = format!("{MOL_A}$$$$\n{MOL_BLANK_NAME}$$$$\n{MOL_B}$$$$\n");
        let results: Vec<_> = SdfReader::new(&sdf).collect();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].as_ref().expect("mol_a").1.name, "mol_a");
        let (mol1, meta1) = results[1].as_ref().expect("blank-name record parse");
        assert_eq!(meta1.name, "");
        assert_eq!(mol1.atom_count(), 2);
        assert_eq!(results[2].as_ref().expect("mol_b").1.name, "mol_b");
    }

    #[test]
    fn test_sdf_reader_blank_name_last_record() {
        let sdf = format!("{MOL_A}$$$$\n{MOL_BLANK_NAME}$$$$\n");
        let results: Vec<_> = SdfReader::new(&sdf).collect();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].as_ref().expect("mol_a").1.name, "mol_a");
        let (mol1, meta1) = results[1].as_ref().expect("blank-name record parse");
        assert_eq!(meta1.name, "");
        assert_eq!(mol1.atom_count(), 2);
    }

    #[test]
    fn test_sdf_reader_blank_name_crlf() {
        let sdf = crlf(&format!("{MOL_A}$$$$\n{MOL_BLANK_NAME}$$$$\n{MOL_B}$$$$\n"));
        let results: Vec<_> = SdfReader::new(&sdf).collect();
        assert_eq!(results.len(), 3);
        let (mol1, meta1) = results[1].as_ref().expect("blank-name record parse (CRLF)");
        assert_eq!(meta1.name, "");
        assert_eq!(mol1.atom_count(), 2);
    }

    #[test]
    fn test_sdf_record_reader_blank_name_first_middle_last() {
        // first
        let sdf_first = format!("{MOL_BLANK_NAME}$$$$\n{MOL_A}$$$$\n");
        let recs: Vec<_> = SdfRecordReader::new(&sdf_first)
            .map(|r| r.expect("parse"))
            .collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].meta.name, "");
        assert_eq!(recs[0].mol.atom_count(), 2);

        // middle
        let sdf_middle = format!("{MOL_A}$$$$\n{MOL_BLANK_NAME}$$$$\n{MOL_B}$$$$\n");
        let recs: Vec<_> = SdfRecordReader::new(&sdf_middle)
            .map(|r| r.expect("parse"))
            .collect();
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[1].meta.name, "");
        assert_eq!(recs[1].mol.atom_count(), 2);

        // last
        let sdf_last = format!("{MOL_A}$$$$\n{MOL_BLANK_NAME}$$$$\n");
        let recs: Vec<_> = SdfRecordReader::new(&sdf_last)
            .map(|r| r.expect("parse"))
            .collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[1].meta.name, "");
        assert_eq!(recs[1].mol.atom_count(), 2);
    }

    #[test]
    fn test_sdf_reader_malformed_input_recovery_unaffected() {
        // Existing malformed-input behavior (bad counts line) must still
        // error, confirming the blank-line-skip removal didn't loosen error
        // recovery. Same shape as test_sdf_reader_stops_on_error.
        let bad_sdf = format!("{MOL_A}$$$$\nbad\n  prog\n\n  X  Y\nM  END\n$$$$\n");
        assert!(parse_sdf(&bad_sdf).is_err());
    }

    // ── stereo diagnostics: direct parse vs. SDF supplier ────────────────

    /// A valid CHFClBr wedge block, generated via the crate's own writer
    /// (not hand-typed fixed-width text) so column layout can't drift from
    /// what the parser actually expects.
    fn wedge_mol_block() -> String {
        use crate::mol2000::{MolMetadata, write_mol_with_coords};
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (-1.0, 0.4),
            (0.9, 0.7),
            (-0.5, -1.1),
            (0.8, -0.6),
        ];
        write_mol_with_coords(&mol, &MolMetadata::default().with_name("wedge"), &coords)
    }

    /// A contradictory-wedge CHFClBr block (two disagreeing wedges) --
    /// same shape used by `chematic_perception::stereo2d_local`'s own
    /// `contradictory_wedges_no_assignment` fixture.
    fn contradictory_mol_block() -> String {
        use crate::mol2000::{MolMetadata, write_mol_with_coords};
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Up).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            (0.0, 0.0),
            (-1.0, 0.4),
            (0.9, 0.7),
            (-0.5, -1.1),
            (0.8, -0.6),
        ];
        write_mol_with_coords(&mol, &MolMetadata::default().with_name("bad"), &coords)
    }

    #[test]
    fn sdf_record_reader_valid_wedge_matches_direct_parse() {
        let block = wedge_mol_block();
        let direct = read_mol_with_diagnostics(&block).expect("direct parse");
        assert!(direct.stereo_diagnostics.is_empty());

        let sdf = format!("{block}$$$$\n");
        let rec = SdfRecordReader::new(&sdf)
            .next()
            .expect("one record")
            .expect("parse ok");
        assert!(rec.stereo_diagnostics.is_empty());
        assert_eq!(
            rec.mol.atom(chematic_core::AtomIdx(0)).chirality,
            direct.mol.atom(chematic_core::AtomIdx(0)).chirality
        );
    }

    #[test]
    fn sdf_record_reader_contradictory_wedge_diagnostic_matches_direct_parse() {
        let block = contradictory_mol_block();
        let direct = read_mol_with_diagnostics(&block).expect("direct parse");
        assert_eq!(direct.stereo_diagnostics.len(), 1);

        let sdf = format!("{block}$$$$\n");
        let rec = SdfRecordReader::new(&sdf)
            .next()
            .expect("one record")
            .expect("parse ok");
        assert_eq!(rec.stereo_diagnostics, direct.stereo_diagnostics);
    }

    #[test]
    fn sdf_file_reader_diagnostics_match_direct_parse() {
        use std::io::{BufReader, Cursor};

        let block = contradictory_mol_block();
        let direct = read_mol_with_diagnostics(&block).expect("direct parse");

        let sdf = format!("{block}$$$$\n");
        let cursor = Cursor::new(sdf.into_bytes());
        let rec = SdfFileReader::new(BufReader::new(cursor))
            .next()
            .expect("one record")
            .expect("parse ok");
        assert_eq!(rec.stereo_diagnostics, direct.stereo_diagnostics);
    }
}
