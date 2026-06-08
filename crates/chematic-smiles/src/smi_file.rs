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
    let mut results = Vec::new();
    for line in s.lines() {
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
}
