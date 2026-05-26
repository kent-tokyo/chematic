//! MOL V2000 (Ctab) parser and writer.
//!
//! Reference format:
//!   Line 1  — molecule name (may be blank)
//!   Line 2  — program/date info (may be blank)
//!   Line 3  — comment (may be blank)
//!   Line 4  — counts line: fixed-width fields for atom count, bond count, version tag
//!   Lines 5..5+natoms — atom block (one line per atom)
//!   Lines 5+natoms..5+natoms+nbonds — bond block (one line per bond)
//!   "M  END" — molecule terminator

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

use crate::error::MolParseError;

/// Metadata extracted from the three-line MOL header.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MolMetadata {
    /// Molecule name from header line 1.
    pub name: String,
    /// Comment string from header line 3.
    pub comment: String,
}

// ---------------------------------------------------------------------------
// Charge encoding table (V2000 ccc field → formal charge)
// ---------------------------------------------------------------------------

/// Decode a V2000 charge code into a formal charge value.
fn decode_charge(code: i8) -> i8 {
    match code {
        1 => 3,
        2 => 2,
        3 => 1,
        4 => 0, // doublet radical — treated as neutral
        5 => -1,
        6 => -2,
        7 => -3,
        _ => 0,
    }
}

/// Encode a formal charge into a V2000 charge code.
fn encode_charge(charge: i8) -> u8 {
    match charge {
        3 => 1,
        2 => 2,
        1 => 3,
        -1 => 5,
        -2 => 6,
        -3 => 7,
        _ => 0,
    }
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a fixed-width 3-character integer field from a string slice.
///
/// Returns an error using `make_err` when the slice is missing or the text
/// cannot be parsed as an integer.
fn parse_field3(
    line: &str,
    start: usize,
    line_num: usize,
    make_err: impl Fn(usize, String) -> MolParseError,
) -> Result<usize, MolParseError> {
    let field = line.get(start..start + 3).ok_or_else(|| {
        make_err(line_num, format!("line too short at column {start}"))
    })?;
    field.trim().parse::<usize>().map_err(|_| {
        make_err(line_num, format!("cannot parse integer from '{field}'"))
    })
}

/// Parse a MOL V2000 string into a `(Molecule, MolMetadata)` pair.
///
/// The parser follows the MDL/CTfile fixed-width column layout.
/// Coordinates are parsed but not stored (the core `Molecule` type has no
/// coordinate fields yet).
pub fn parse_mol(input: &str) -> Result<(Molecule, MolMetadata), MolParseError> {
    let mut lines = input.lines().enumerate().peekable();

    // -- Header block: lines 1–3 -------------------------------------------

    let name = match lines.next() {
        Some((_, l)) => l.to_string(),
        None => return Err(MolParseError::UnexpectedEnd),
    };

    // Line 2: program/date — ignored but must be present
    match lines.next() {
        Some(_) => {}
        None => return Err(MolParseError::UnexpectedEnd),
    }

    let comment = match lines.next() {
        Some((_, l)) => l.to_string(),
        None => return Err(MolParseError::UnexpectedEnd),
    };

    let metadata = MolMetadata { name, comment };

    // -- Counts line (line 4) -----------------------------------------------

    let (counts_lineno, counts_line) = match lines.next() {
        Some((i, l)) => (i + 1, l), // 1-based for error messages
        None => return Err(MolParseError::UnexpectedEnd),
    };

    // Verify V2000 version tag (bytes 33–38 in the spec; be lenient with
    // shorter lines — just check that "V2000" appears somewhere in the line).
    if !counts_line.contains("V2000") {
        return Err(MolParseError::InvalidCountLine {
            line: counts_lineno,
            detail: "missing V2000 version tag".to_string(),
        });
    }

    let make_count_err = |ln: usize, d: String| MolParseError::InvalidCountLine { line: ln, detail: d };

    let natoms = parse_field3(counts_line, 0, counts_lineno, make_count_err)?;
    let nbonds = parse_field3(counts_line, 3, counts_lineno, make_count_err)?;

    // -- Atom block ---------------------------------------------------------

    let mut builder = MoleculeBuilder::new();

    for atom_i in 0..natoms {
        let (raw_lineno, atom_line) = match lines.next() {
            Some((i, l)) => (i + 1, l),
            None => return Err(MolParseError::UnexpectedEnd),
        };

        let make_atom_err = |ln: usize, d: String| MolParseError::InvalidAtomLine { line: ln, detail: d };

        // Element symbol: bytes 31–33 (3 chars, left-padded with a space in
        // the spec, but writers vary; trim both ends).
        let sym = atom_line.get(31..34).ok_or_else(|| {
            make_atom_err(raw_lineno, format!("atom line {atom_i} too short for element field"))
        })?;
        let sym = sym.trim();

        let element = Element::from_symbol(sym).ok_or_else(|| MolParseError::UnknownElement {
            symbol: sym.to_string(),
            line: raw_lineno,
        })?;

        // Charge code: bytes 36–38 (3 chars).
        let charge = if let Some(ccc) = atom_line.get(36..39) {
            let code: i8 = ccc.trim().parse().unwrap_or(0);
            decode_charge(code)
        } else {
            0
        };

        let mut atom = Atom::new(element);
        atom.charge = charge;
        builder.add_atom(atom);
    }

    // -- Bond block ---------------------------------------------------------

    for bond_i in 0..nbonds {
        let (raw_lineno, bond_line) = match lines.next() {
            Some((i, l)) => (i + 1, l),
            None => return Err(MolParseError::UnexpectedEnd),
        };

        let make_bond_err = |ln: usize, d: String| MolParseError::InvalidBondLine { line: ln, detail: d };

        let a1_raw = parse_field3(bond_line, 0, raw_lineno, make_bond_err)?;
        let a2_raw = parse_field3(bond_line, 3, raw_lineno, make_bond_err)?;
        let btype_raw = parse_field3(bond_line, 6, raw_lineno, make_bond_err)?;

        if a1_raw == 0 || a2_raw == 0 {
            return Err(MolParseError::InvalidBondLine {
                line: raw_lineno,
                detail: format!("bond {bond_i}: atom indices are 1-based; got {a1_raw}/{a2_raw}"),
            });
        }

        let a1 = AtomIdx((a1_raw - 1) as u32);
        let a2 = AtomIdx((a2_raw - 1) as u32);

        // Stereo field (columns 9-11, 0-indexed): only meaningful for single bonds.
        let stereo_raw: usize = if bond_line.len() >= 12 {
            parse_field3(bond_line, 9, raw_lineno, make_bond_err).unwrap_or(0)
        } else {
            0
        };

        let order = match btype_raw {
            1 => match stereo_raw {
                1 | 4 => BondOrder::Up,
                6 => BondOrder::Down,
                _ => BondOrder::Single,
            },
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Aromatic,
            // 5 = single-or-double, 6 = single-or-aromatic, 7 = double-or-aromatic, 8 = any
            // Fall back to Single for query bond types not representable in BondOrder.
            _ => BondOrder::Single,
        };

        builder.add_bond(a1, a2, order).map_err(|e| MolParseError::InvalidBondLine {
            line: raw_lineno,
            detail: format!("bond {bond_i}: {e}"),
        })?;
    }

    // Consume lines until "M  END" (or end of iterator — the block may have
    // property lines before M  END which we skip).
    for (_i, l) in lines.by_ref() {
        if l.trim_start().starts_with("M  END") {
            break;
        }
    }

    Ok((builder.build(), metadata))
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write a `Molecule` to MOL V2000 format.
///
/// Coordinates are written as 0.0 because the core `Molecule` type does not
/// store 2D/3D coordinates.  All other atom and bond fields are derived from
/// the molecule graph.
pub fn write_mol(mol: &Molecule, metadata: &MolMetadata) -> String {
    let mut out = String::new();

    // Header lines 1–3
    out.push_str(&metadata.name);
    out.push('\n');
    out.push_str("  chematic\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    // Counts line (line 4)
    // Format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    // We write atom count, bond count, then pad the remaining fixed fields and
    // append the V2000 version tag.
    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        natoms, nbonds
    ));

    // Atom block
    for (_idx, atom) in mol.atoms() {
        let sym = atom.element.symbol();
        let charge_code = encode_charge(atom.charge);

        // x y z coords (10.4 format, 3 fields), then space, element symbol
        // (left-aligned, 3 chars), then mass diff (2 chars), charge (3 chars),
        // then trailing zeroed fields matching the standard atom line.
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0{:>3}  0  0  0  0  0  0  0  0  0\n",
            0.0_f64, 0.0_f64, 0.0_f64,
            sym,
            charge_code,
        ));
    }

    // Bond block
    for (_idx, bond) in mol.bonds() {
        let a1 = bond.atom1.0 + 1; // convert to 1-based
        let a2 = bond.atom2.0 + 1;
        let btype = match bond.order {
            BondOrder::Aromatic => 4,
            _ => bond.order.order_int(),
        };
        // Bond line: 111222tttsssxxxrrrccc (3-char fields)
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}  0\n",
            a1, a2, btype
        ));
    }

    // Terminator
    out.push_str("M  END\n");

    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal ethanol MOL V2000 block (CCO, 3 atoms, 2 bonds).
    const ETHANOL_MOL: &str = "\
ethanol
  chematic

  3  2  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
";

    #[test]
    fn test_parse_ethanol_counts() {
        let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse should succeed");
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        assert_eq!(meta.name, "ethanol");
    }

    #[test]
    fn test_parse_elements() {
        let (mol, _) = parse_mol(ETHANOL_MOL).expect("parse should succeed");
        let atoms: Vec<_> = mol.atoms().collect();
        assert_eq!(atoms[0].1.element, Element::C);
        assert_eq!(atoms[1].1.element, Element::C);
        assert_eq!(atoms[2].1.element, Element::O);
    }

    #[test]
    fn test_parse_bond_types() {
        // Two carbons: single, double, triple, aromatic bonds.
        let mol_str = "\
test
  chematic

  8  4  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  3  4  2  0
  5  6  3  0
  7  8  4  0
M  END
";
        let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
        let bonds: Vec<_> = mol.bonds().collect();
        assert_eq!(bonds[0].1.order, BondOrder::Single);
        assert_eq!(bonds[1].1.order, BondOrder::Double);
        assert_eq!(bonds[2].1.order, BondOrder::Triple);
        assert_eq!(bonds[3].1.order, BondOrder::Aromatic);
    }

    #[test]
    fn test_parse_charge() {
        // Nitrogen with charge code 3 (+1 formal charge).
        let mol_str = "\
charged
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
M  END
";
        let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
        assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
    }

    #[test]
    fn test_parse_negative_charge() {
        // Oxygen with charge code 5 (-1 formal charge).
        let mol_str = "\
negcharge
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  END
";
        let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
        assert_eq!(mol.atom(AtomIdx(0)).charge, -1);
    }

    #[test]
    fn test_round_trip() {
        // Parse → write → parse again; atom and bond counts must match.
        let (mol1, meta1) = parse_mol(ETHANOL_MOL).expect("first parse");
        let written = write_mol(&mol1, &meta1);
        let (mol2, _meta2) = parse_mol(&written).expect("second parse");
        assert_eq!(mol1.atom_count(), mol2.atom_count());
        assert_eq!(mol1.bond_count(), mol2.bond_count());
    }

    #[test]
    fn test_round_trip_elements_preserved() {
        let (mol1, meta1) = parse_mol(ETHANOL_MOL).expect("first parse");
        let written = write_mol(&mol1, &meta1);
        let (mol2, _) = parse_mol(&written).expect("second parse");
        for ((_, a1), (_, a2)) in mol1.atoms().zip(mol2.atoms()) {
            assert_eq!(a1.element, a2.element);
        }
    }

    #[test]
    fn test_error_missing_v2000() {
        let bad = "\
bad
  prog

  3  2  0  0  0  0  0  0  0  0  0 V3000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
";
        assert!(matches!(
            parse_mol(bad),
            Err(MolParseError::InvalidCountLine { .. })
        ));
    }

    #[test]
    fn test_error_truncated_input() {
        // Counts line says 3 atoms but only 1 is provided.
        let bad = "\
trunc
  prog

  3  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
";
        assert!(matches!(parse_mol(bad), Err(MolParseError::UnexpectedEnd)));
    }

    #[test]
    fn test_error_invalid_counts_line() {
        // Counts line too short (no V2000 tag at all).
        let bad = "\
mol
  prog

  X  Y
M  END
";
        assert!(matches!(
            parse_mol(bad),
            Err(MolParseError::InvalidCountLine { .. })
        ));
    }

    #[test]
    fn test_write_contains_m_end() {
        let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
        let written = write_mol(&mol, &meta);
        assert!(written.contains("M  END"));
    }

    #[test]
    fn test_write_contains_v2000() {
        let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
        let written = write_mol(&mol, &meta);
        assert!(written.contains("V2000"));
    }

    #[test]
    fn test_parse_stereo_up_bond() {
        // MOL V2000 with a stereo=1 (Up) bond
        let mol_str = "\n\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  1  0  0  0\nM  END\n";
        let (mol, _) = crate::parse_mol(mol_str).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(bond.order, chematic_core::BondOrder::Up);
    }

    #[test]
    fn test_parse_stereo_down_bond() {
        let mol_str = "\n\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  6  0  0  0\nM  END\n";
        let (mol, _) = crate::parse_mol(mol_str).unwrap();
        let bond = mol.bond(chematic_core::BondIdx(0));
        assert_eq!(bond.order, chematic_core::BondOrder::Down);
    }
}
