//! Multi-molecule SMILES file (.smi) reader and writer.
//!
//! The `.smi` format is one molecule per line:
//!
//! ```text
//! CC\tethane
//! CCO\tethanol
//! c1ccccc1\tbenzene
//! ```
//!
//! The separator between SMILES and name is a tab **or** one or more spaces.
//! Lines starting with `#` and blank lines are silently skipped.

use crate::error::SmilesError;
use crate::parser::parse;
use crate::writer::write;
use chematic_core::Molecule;

/// Resource limits for multi-molecule `.smi` input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmiFileParseLimits {
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_records: usize,
}

impl Default for SmiFileParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_bytes: 16 << 20,
            max_records: 100_000,
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse a multi-molecule `.smi` string.
///
/// Each element is `Ok((molecule, name))` for a successfully parsed line, or
/// `Err(SmilesError)` for a malformed SMILES.  Blank lines and `#` comments
/// are skipped and do not produce entries.
///
/// If a line has no name field, the name is an empty string.
pub fn parse_smi_file(s: &str) -> Vec<Result<(Molecule, String), SmilesError>> {
    parse_smi_file_with_limits(s, &SmiFileParseLimits::default())
}

/// Parse a multi-molecule `.smi` string with explicit resource limits.
pub fn parse_smi_file_with_limits(
    s: &str,
    limits: &SmiFileParseLimits,
) -> Vec<Result<(Molecule, String), SmilesError>> {
    if s.len() > limits.max_input_bytes {
        return vec![Err(SmilesError::ResourceLimit {
            resource: "smi file input bytes",
            actual: s.len(),
            limit: limits.max_input_bytes,
        })];
    }
    let mut results = Vec::new();
    for line in s.lines() {
        if line.len() > limits.max_line_bytes {
            results.push(Err(SmilesError::ResourceLimit {
                resource: "smi file line bytes",
                actual: line.len(),
                limit: limits.max_line_bytes,
            }));
            break;
        }
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // Split on first whitespace run: first token = SMILES, rest = name.
        let mut parts = line.splitn(2, |c: char| c.is_whitespace());
        let smiles = parts.next().unwrap_or("").trim();
        let name = parts.next().unwrap_or("").trim().to_string();
        if smiles.is_empty() {
            continue;
        }
        if results.len() >= limits.max_records {
            results.push(Err(SmilesError::ResourceLimit {
                resource: "smi file records",
                actual: results.len() + 1,
                limit: limits.max_records,
            }));
            break;
        }
        results.push(parse(smiles).map(|mol| (mol, name)));
    }
    results
}

/// Write a list of `(molecule, name)` pairs to `.smi` format.
///
/// Each molecule is written as `SMILES<TAB>name\n`.
/// If `name` is empty, the tab and name are omitted.
pub fn write_smi_file(records: &[(Molecule, &str)]) -> String {
    let mut out = String::new();
    for (mol, name) in records {
        let smiles = write(mol);
        if name.is_empty() {
            out.push_str(&smiles);
        } else {
            out.push_str(&smiles);
            out.push('\t');
            out.push_str(name);
        }
        out.push('\n');
    }
    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tab_separated() {
        let s = "CC\tethane\nCCO\tethanol\n";
        let results = parse_smi_file(s);
        assert_eq!(results.len(), 2);
        let (mol0, name0) = results[0].as_ref().unwrap();
        assert_eq!(mol0.atom_count(), 2);
        assert_eq!(name0, "ethane");
        let (mol1, name1) = results[1].as_ref().unwrap();
        assert_eq!(mol1.atom_count(), 3);
        assert_eq!(name1, "ethanol");
    }

    #[test]
    fn test_parse_space_separated() {
        let s = "CC ethane\nc1ccccc1 benzene\n";
        let results = parse_smi_file(s);
        assert_eq!(results.len(), 2);
        assert!(results[0].is_ok());
        assert_eq!(results[0].as_ref().unwrap().1, "ethane");
    }

    #[test]
    fn test_parse_skips_comments_and_blanks() {
        let s = "# comment\n\nCC\tethane\n# another\n";
        let results = parse_smi_file(s);
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn test_parse_no_name() {
        let s = "CC\n";
        let results = parse_smi_file(s);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].as_ref().unwrap().1, "");
    }

    #[test]
    fn test_parse_invalid_smiles_is_err() {
        // Unclosed ring closure is invalid SMILES.
        let s = "C1CC\tbad\n";
        let results = parse_smi_file(s);
        assert_eq!(results.len(), 1);
        assert!(results[0].is_err(), "unclosed ring should be Err");
    }

    #[test]
    fn test_write_roundtrip() {
        use crate::parser::parse as parse_smiles;
        let benzene = parse_smiles("c1ccccc1").unwrap();
        let ethane = parse_smiles("CC").unwrap();
        let records: Vec<(Molecule, &str)> = vec![(benzene, "benzene"), (ethane, "ethane")];
        let s = write_smi_file(&records);
        let back = parse_smi_file(&s);
        assert_eq!(back.len(), 2);
        assert_eq!(back[0].as_ref().unwrap().1, "benzene");
        assert_eq!(back[1].as_ref().unwrap().1, "ethane");
    }

    #[test]
    fn test_parse_limits_report_file_resource_errors() {
        assert!(matches!(
            parse_smi_file_with_limits(
                "CC\n",
                &SmiFileParseLimits {
                    max_input_bytes: 1,
                    ..Default::default()
                }
            )[0],
            Err(SmilesError::ResourceLimit {
                resource: "smi file input bytes",
                ..
            })
        ));
        assert!(matches!(
            parse_smi_file_with_limits(
                "CC name\n",
                &SmiFileParseLimits {
                    max_line_bytes: 3,
                    ..Default::default()
                }
            )[0],
            Err(SmilesError::ResourceLimit {
                resource: "smi file line bytes",
                ..
            })
        ));
        assert!(matches!(
            parse_smi_file_with_limits(
                "CC\nCCO\n",
                &SmiFileParseLimits {
                    max_records: 1,
                    ..Default::default()
                }
            )[1],
            Err(SmilesError::ResourceLimit {
                resource: "smi file records",
                ..
            })
        ));
    }
}
