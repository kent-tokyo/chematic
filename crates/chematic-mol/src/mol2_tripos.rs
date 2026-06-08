//! Tripos MOL2 format reader and writer.
//!
//! Supports the three mandatory sections of the Tripos MOL2 format:
//! `@<TRIPOS>MOLECULE`, `@<TRIPOS>ATOM`, and `@<TRIPOS>BOND`.
//!
//! Reference: Tripos MOL2 format specification (SYBYL 7.x).

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Error returned when parsing a Tripos MOL2 file fails.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Mol2Error {
    /// Required section not found.
    MissingSection(String),
    /// An atom line is malformed.
    InvalidAtomLine { line: usize, detail: String },
    /// A bond line is malformed.
    InvalidBondLine { line: usize, detail: String },
    /// Unknown element symbol.
    UnknownElement { symbol: String, line: usize },
}

impl core::fmt::Display for Mol2Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingSection(s) => write!(f, "MOL2: missing section @<TRIPOS>{s}"),
            Self::InvalidAtomLine { line, detail } => {
                write!(f, "MOL2: invalid atom line {line}: {detail}")
            }
            Self::InvalidBondLine { line, detail } => {
                write!(f, "MOL2: invalid bond line {line}: {detail}")
            }
            Self::UnknownElement { symbol, line } => {
                write!(f, "MOL2: unknown element '{symbol}' at line {line}")
            }
        }
    }
}

impl std::error::Error for Mol2Error {}

// ---------------------------------------------------------------------------
// Parser helpers
// ---------------------------------------------------------------------------

/// Extract all lines belonging to a named `@<TRIPOS>SECTION` block.
///
/// Returns the (1-based line number, content) pairs starting immediately
/// after the section header, stopping at the next `@<TRIPOS>` header or EOF.
fn section_lines<'a>(lines: &'a [(usize, &'a str)], name: &str) -> Vec<(usize, &'a str)> {
    let header = format!("@<TRIPOS>{name}");
    let mut in_section = false;
    let mut result = Vec::new();
    for &(lineno, line) in lines {
        let trimmed = line.trim();
        if trimmed.eq_ignore_ascii_case(&header) {
            in_section = true;
            continue;
        }
        if in_section {
            if trimmed.starts_with("@<TRIPOS>") {
                break;
            }
            // Skip blank lines and comment lines (starting with #).
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                result.push((lineno, line));
            }
        }
    }
    result
}

/// Strip a Tripos atom type suffix (e.g. `C.ar` → `C`, `N.am` → `N`, `N.3` → `N`).
fn strip_atom_type(sym: &str) -> &str {
    sym.split('.').next().unwrap_or(sym)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse a Tripos MOL2 string into a `(Molecule, 3D-coordinates)` pair.
///
/// Only the `MOLECULE`, `ATOM`, and `BOND` sections are required; all other
/// sections (SUBSTRUCTURE, SET, …) are silently ignored.
///
/// 3D coordinates are returned as `Vec<(f64, f64, f64)>` aligned with the
/// molecule's atom indices.
#[allow(clippy::type_complexity)]
pub fn parse_mol2(s: &str) -> Result<(Molecule, Vec<(f64, f64, f64)>), Mol2Error> {
    let all_lines: Vec<(usize, &str)> = s.lines().enumerate().map(|(i, l)| (i + 1, l)).collect();

    // -- MOLECULE section: we read it but only need it for sanity; skip for now.

    // -- ATOM section ----------------------------------------------------------
    let atom_lines = section_lines(&all_lines, "ATOM");
    if atom_lines.is_empty() {
        return Err(Mol2Error::MissingSection("ATOM".into()));
    }

    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64, f64)> = Vec::new();
    // Map from MOL2 1-based atom_id → builder AtomIdx.
    let mut atom_id_map: Vec<(u32, AtomIdx)> = Vec::new();

    for (lineno, line) in &atom_lines {
        // MOL2 ATOM line format (space-separated):
        // atom_id  atom_name  x  y  z  atom_type  [subst_id  subst_name  charge]
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 6 {
            return Err(Mol2Error::InvalidAtomLine {
                line: *lineno,
                detail: format!("expected at least 6 fields, got {}", parts.len()),
            });
        }

        let atom_id: u32 = parts[0].parse().map_err(|_| Mol2Error::InvalidAtomLine {
            line: *lineno,
            detail: format!("cannot parse atom_id from '{}'", parts[0]),
        })?;

        let x: f64 = parts[2].parse().map_err(|_| Mol2Error::InvalidAtomLine {
            line: *lineno,
            detail: format!("cannot parse x from '{}'", parts[2]),
        })?;
        let y: f64 = parts[3].parse().map_err(|_| Mol2Error::InvalidAtomLine {
            line: *lineno,
            detail: format!("cannot parse y from '{}'", parts[3]),
        })?;
        let z: f64 = parts[4].parse().map_err(|_| Mol2Error::InvalidAtomLine {
            line: *lineno,
            detail: format!("cannot parse z from '{}'", parts[4]),
        })?;

        // atom_type may be "C.ar", "N.3", "O.2", etc. Strip the suffix.
        let atom_type_raw = parts[5];
        let sym_raw = strip_atom_type(atom_type_raw);
        // Capitalise first letter for element lookup.
        let sym = {
            let mut c = sym_raw.chars();
            match c.next() {
                None => String::new(),
                Some(first) => first.to_uppercase().collect::<String>() + c.as_str(),
            }
        };

        // Optional partial charge (field 9, 0-indexed as parts[8]).
        let charge_f: f64 = parts.get(8).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        // Round to nearest integer for formal charge (MOL2 uses partial charges).
        let formal_charge: i8 = charge_f.round() as i8;

        let element = Element::from_symbol(&sym).ok_or_else(|| Mol2Error::UnknownElement {
            symbol: sym.clone(),
            line: *lineno,
        })?;

        let mut atom = Atom::new(element);
        atom.charge = formal_charge;

        let builder_idx = builder.add_atom(atom);
        atom_id_map.push((atom_id, builder_idx));
        coords.push((x, y, z));
    }

    // -- BOND section ----------------------------------------------------------
    let bond_lines = section_lines(&all_lines, "BOND");
    // BOND section is optional (0-atom molecules are valid).

    for (lineno, line) in &bond_lines {
        // bond_id  origin_atom_id  target_atom_id  bond_type
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: format!("expected at least 4 fields, got {}", parts.len()),
            });
        }

        let a1_id: u32 = parts[1].parse().map_err(|_| Mol2Error::InvalidBondLine {
            line: *lineno,
            detail: format!("cannot parse origin atom_id from '{}'", parts[1]),
        })?;
        let a2_id: u32 = parts[2].parse().map_err(|_| Mol2Error::InvalidBondLine {
            line: *lineno,
            detail: format!("cannot parse target atom_id from '{}'", parts[2]),
        })?;

        let a1 = atom_id_map
            .iter()
            .find(|&&(k, _)| k == a1_id)
            .map(|&(_, v)| v)
            .ok_or_else(|| Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: format!("atom_id {a1_id} not found"),
            })?;
        let a2 = atom_id_map
            .iter()
            .find(|&&(k, _)| k == a2_id)
            .map(|&(_, v)| v)
            .ok_or_else(|| Mol2Error::InvalidBondLine {
                line: *lineno,
                detail: format!("atom_id {a2_id} not found"),
            })?;

        let bond_type = parts[3];
        let order = match bond_type {
            "1" | "1.5" => BondOrder::Single,
            "2" => BondOrder::Double,
            "3" => BondOrder::Triple,
            "ar" | "am" => BondOrder::Aromatic,
            "un" => BondOrder::QueryAny,
            "du" | "nc" => BondOrder::Zero,
            _ => BondOrder::Single,
        };

        // Ignore duplicate bond errors (some MOL2 files repeat bonds).
        let _ = builder.add_bond(a1, a2, order);
    }

    Ok((builder.build(), coords))
}

/// Write a molecule and its 3D coordinates to Tripos MOL2 format.
///
/// `coords` must have one entry per atom in `mol`.  Missing entries use
/// `(0.0, 0.0, 0.0)`.
pub fn write_mol2(mol: &Molecule, coords: &[(f64, f64, f64)]) -> String {
    let mut out = String::new();

    // MOLECULE section.
    out.push_str("@<TRIPOS>MOLECULE\n");
    out.push_str("chematic\n");
    out.push_str(&format!(
        "{} {} 0 0 0\n",
        mol.atom_count(),
        mol.bond_count()
    ));
    out.push_str("SMALL\n");
    out.push_str("GASTEIGER\n\n");

    // ATOM section.
    out.push_str("@<TRIPOS>ATOM\n");
    for (idx, atom) in mol.atoms() {
        let i = idx.0 + 1; // 1-based
        let sym = atom.element.symbol();
        let (x, y, z) = coords
            .get(idx.0 as usize)
            .copied()
            .unwrap_or((0.0, 0.0, 0.0));
        // Use symbol as atom type (simplified; no hybridisation-based suffix).
        out.push_str(&format!(
            "{i:>6} {sym:<4} {x:>10.4} {y:>10.4} {z:>10.4} {sym:<8} 1  LIG  0.0000\n"
        ));
    }

    // BOND section.
    out.push_str("@<TRIPOS>BOND\n");
    for (bidx, bond) in mol.bonds() {
        let bi = bidx.0 + 1;
        let a1 = bond.atom1.0 + 1;
        let a2 = bond.atom2.0 + 1;
        let btype = match bond.order {
            BondOrder::Zero => "nc",
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => "1",
            BondOrder::Double => "2",
            BondOrder::Triple => "3",
            BondOrder::Aromatic => "ar",
            BondOrder::Quadruple => "4",
            BondOrder::QueryAny
            | BondOrder::QuerySingleOrDouble
            | BondOrder::QuerySingleOrAromatic
            | BondOrder::QueryDoubleOrAromatic => "un",
        };
        out.push_str(&format!("{bi:>6} {a1:>6} {a2:>6} {btype}\n"));
    }

    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const ETHANOL_MOL2: &str = "\
@<TRIPOS>MOLECULE
ethanol
 3 2 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
      2 C2          1.5000    0.0000    0.0000 C.3     1  LIG1        0.0000
      3 O1          3.0000    0.0000    0.0000 O.3     1  LIG1       -0.3940
@<TRIPOS>BOND
     1     1     2    1
     2     2     3    1
";

    #[test]
    fn test_parse_mol2_atom_count() {
        let (mol, coords) = parse_mol2(ETHANOL_MOL2).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        assert_eq!(coords.len(), 3);
    }

    #[test]
    fn test_parse_mol2_coords() {
        let (_, coords) = parse_mol2(ETHANOL_MOL2).unwrap();
        assert!((coords[0].0 - 0.0).abs() < 1e-4);
        assert!((coords[1].0 - 1.5).abs() < 1e-4);
        assert!((coords[2].0 - 3.0).abs() < 1e-4);
    }

    #[test]
    fn test_parse_mol2_elements() {
        let (mol, _) = parse_mol2(ETHANOL_MOL2).unwrap();
        assert_eq!(mol.atom(AtomIdx(0)).element.symbol(), "C");
        assert_eq!(mol.atom(AtomIdx(1)).element.symbol(), "C");
        assert_eq!(mol.atom(AtomIdx(2)).element.symbol(), "O");
    }

    #[test]
    fn test_write_mol2_roundtrip() {
        let (mol, coords) = parse_mol2(ETHANOL_MOL2).unwrap();
        let written = write_mol2(&mol, &coords);
        let (mol2, coords2) = parse_mol2(&written).unwrap();
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert_eq!(mol.bond_count(), mol2.bond_count());
        // Coordinates should survive round-trip.
        for (a, b) in coords.iter().zip(coords2.iter()) {
            assert!((a.0 - b.0).abs() < 0.01);
            assert!((a.1 - b.1).abs() < 0.01);
            assert!((a.2 - b.2).abs() < 0.01);
        }
    }

    #[test]
    fn test_aromatic_bond_type() {
        let benzene_mol2 = "\
@<TRIPOS>MOLECULE
benzene
 6 6 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C1  0.0000  1.4000  0.0000 C.ar  1  BNZ  0.0000
      2 C2  1.2124  0.7000  0.0000 C.ar  1  BNZ  0.0000
      3 C3  1.2124 -0.7000  0.0000 C.ar  1  BNZ  0.0000
      4 C4  0.0000 -1.4000  0.0000 C.ar  1  BNZ  0.0000
      5 C5 -1.2124 -0.7000  0.0000 C.ar  1  BNZ  0.0000
      6 C6 -1.2124  0.7000  0.0000 C.ar  1  BNZ  0.0000
@<TRIPOS>BOND
     1     1     2   ar
     2     2     3   ar
     3     3     4   ar
     4     4     5   ar
     5     5     6   ar
     6     6     1   ar
";
        let (mol, _) = parse_mol2(benzene_mol2).unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        // All bonds should be aromatic.
        for (_, bond) in mol.bonds() {
            assert_eq!(bond.order, BondOrder::Aromatic);
        }
    }

    #[test]
    fn test_missing_atom_section() {
        let bad = "@<TRIPOS>MOLECULE\nbad\n";
        assert!(parse_mol2(bad).is_err());
    }
}
