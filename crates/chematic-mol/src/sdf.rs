//! SDF (Structure-Data File) reader.
//!
//! An SDF file contains one or more MOL V2000 blocks separated by `$$$$`
//! delimiter lines.  Data-field sections between `M  END` and `$$$$` are
//! accepted but ignored.

use chematic_core::Molecule;

use crate::error::MolParseError;
use crate::mol2000::{MolMetadata, parse_mol};

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
        // Skip leading blank lines between records (defensive; well-formed SDF
        // should not have them, but some writers emit a trailing blank).
        while self.remaining.starts_with('\n') || self.remaining.starts_with("\r\n") {
            let skip = if self.remaining.starts_with("\r\n") { 2 } else { 1 };
            self.remaining = &self.remaining[skip..];
        }

        if self.remaining.is_empty() {
            return None;
        }

        self.current_mol_num += 1;

        // Locate the `$$$$` delimiter line. We search line by line so that a
        // `$$$$` substring inside a data value does not trigger a false match.
        let mut end_byte = self.remaining.len(); // default: consume everything
        let mut after_delim = "";
        let mut scan = self.remaining;
        let mut byte_offset = 0usize;

        loop {
            // Find the next newline.
            let nl_pos = scan.find('\n');
            let (line, rest_after_nl) = match nl_pos {
                Some(pos) => {
                    let line = &scan[..pos];
                    // Strip a possible leading \r from the line.
                    let line = line.trim_end_matches('\r');
                    (line, &scan[pos + 1..])
                }
                None => {
                    // Last line with no trailing newline.
                    let line = scan.trim_end_matches('\r');
                    (line, "")
                }
            };

            if line == "$$$$" {
                end_byte = byte_offset;
                // `after_delim` starts right after the newline that follows `$$$$`.
                let delim_line_len = if nl_pos.is_some() {
                    line.len() + 1 + 1 // line content + possible \r + \n
                } else {
                    line.len()
                };
                // More precise: skip from byte_offset to past the delimiter line.
                let delim_start = byte_offset;
                let delim_end = match self.remaining[delim_start..].find('\n') {
                    Some(pos) => delim_start + pos + 1,
                    None => self.remaining.len(),
                };
                let _ = delim_line_len; // suppress warning
                after_delim = &self.remaining[delim_end..];
                break;
            }

            if nl_pos.is_none() {
                // Reached the last line without finding `$$$$`.
                break;
            }

            let advance = nl_pos.unwrap() + 1;
            byte_offset += advance;
            scan = rest_after_nl;
        }

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
    SdfReader::new(input).collect()
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
        let sdf_with_data = format!(
            "{MOL_A}> <MW>\n44.0\n\n$$$$\n"
        );
        let results: Vec<_> = SdfReader::new(&sdf_with_data).collect();
        assert_eq!(results.len(), 1);
        let (mol, _) = results[0].as_ref().expect("parse");
        assert_eq!(mol.atom_count(), 2);
    }
}
