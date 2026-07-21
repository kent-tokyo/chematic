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

/// Maximum number of atoms allowed in a MOL file (prevents memory exhaustion).
const MAX_ATOMS: usize = 100_000;

/// Maximum number of bonds allowed in a MOL file (prevents memory exhaustion).
const MAX_BONDS: usize = 200_000;

/// Metadata extracted from the three-line MOL header.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MolMetadata {
    /// Molecule name from header line 1.
    pub name: String,
    /// Comment string from header line 3.
    pub comment: String,
}

impl MolMetadata {
    /// Set the molecule name and return `self` (builder style).
    pub fn with_name(mut self, name: &str) -> Self {
        self.name = name.to_owned();
        self
    }

    /// Set the comment and return `self` (builder style).
    pub fn with_comment(mut self, comment: &str) -> Self {
        self.comment = comment.to_owned();
        self
    }
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
    let field = line
        .get(start..start + 3)
        .ok_or_else(|| make_err(line_num, format!("line too short at column {start}")))?;
    field
        .trim()
        .parse::<usize>()
        .map_err(|_| make_err(line_num, format!("cannot parse integer from '{field}'")))
}

/// Parse a MOL V2000 string into a `(Molecule, MolMetadata, coords)` triple.
///
/// The parser follows the MDL/CTfile fixed-width column layout.
/// `coords[i]` is the `(x, y)` position for atom `i` extracted from the
/// atom block.  Z-coordinates are discarded.
#[allow(clippy::type_complexity)]
pub fn parse_mol_with_coords(
    input: &str,
) -> Result<(Molecule, MolMetadata, Vec<(f64, f64)>), MolParseError> {
    // Yields (1-based line number, line text); short-circuits on EOF.
    let mut lines = input.lines().enumerate().map(|(i, l)| (i + 1, l));
    let mut next_line = || lines.next().ok_or(MolParseError::UnexpectedEnd);

    // -- Header block: lines 1–3 -------------------------------------------

    let name = next_line()?.1.to_string();
    next_line()?; // line 2: program/date — discarded
    let comment = next_line()?.1.to_string();

    let metadata = MolMetadata { name, comment };

    // -- Counts line (line 4) -----------------------------------------------

    let (counts_lineno, counts_line) = next_line()?;

    // Be lenient with shorter lines — just check the V2000 tag exists.
    if !counts_line.contains("V2000") {
        return Err(MolParseError::InvalidCountLine {
            line: counts_lineno,
            detail: "missing V2000 version tag".to_string(),
        });
    }

    let make_count_err = |ln: usize, d: String| MolParseError::InvalidCountLine {
        line: ln,
        detail: d,
    };

    let natoms = parse_field3(counts_line, 0, counts_lineno, make_count_err)?;
    let nbonds = parse_field3(counts_line, 3, counts_lineno, make_count_err)?;

    if natoms > MAX_ATOMS {
        return Err(MolParseError::InvalidCountLine {
            line: counts_lineno,
            detail: format!(
                "atom count {} exceeds maximum allowed {}",
                natoms, MAX_ATOMS
            ),
        });
    }

    if nbonds > MAX_BONDS {
        return Err(MolParseError::InvalidCountLine {
            line: counts_lineno,
            detail: format!(
                "bond count {} exceeds maximum allowed {}",
                nbonds, MAX_BONDS
            ),
        });
    }

    // -- Atom block ---------------------------------------------------------

    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64)> = Vec::with_capacity(natoms);
    let make_atom_err = |ln: usize, d: String| MolParseError::InvalidAtomLine {
        line: ln,
        detail: d,
    };

    for atom_i in 0..natoms {
        let (raw_lineno, atom_line) = next_line()?;

        // Coordinates: bytes 0–9 (x), 10–19 (y), 20–29 (z) — each 10 chars.
        let x: f64 = atom_line
            .get(0..10)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0.0);
        let y: f64 = atom_line
            .get(10..20)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0.0);
        coords.push((x, y));

        // Element symbol: bytes 31–33 (3 chars, left-padded with a space in
        // the spec, but writers vary; trim both ends).
        let sym = atom_line
            .get(31..34)
            .ok_or_else(|| {
                make_atom_err(
                    raw_lineno,
                    format!("atom line {atom_i} too short for element field"),
                )
            })?
            .trim();

        let element = Element::from_symbol(sym).ok_or_else(|| MolParseError::UnknownElement {
            symbol: sym.to_string(),
            line: raw_lineno,
        })?;

        // Charge code: bytes 36–38 (3 chars).
        let charge = atom_line
            .get(36..39)
            .map(|ccc| decode_charge(ccc.trim().parse().unwrap_or(0)))
            .unwrap_or(0);

        let mut atom = Atom::new(element);
        atom.charge = charge;
        builder.add_atom(atom);
    }

    // -- Bond block ---------------------------------------------------------

    let make_bond_err = |ln: usize, d: String| MolParseError::InvalidBondLine {
        line: ln,
        detail: d,
    };

    for bond_i in 0..nbonds {
        let (raw_lineno, bond_line) = next_line()?;

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
            5 => BondOrder::QuerySingleOrDouble,
            6 => BondOrder::QuerySingleOrAromatic,
            7 => BondOrder::QueryDoubleOrAromatic,
            8 => BondOrder::QueryAny,
            _ => BondOrder::Single,
        };

        builder
            .add_bond(a1, a2, order)
            .map_err(|e| MolParseError::InvalidBondLine {
                line: raw_lineno,
                detail: format!("bond {bond_i}: {e}"),
            })?;
    }

    // Skip property lines until "M  END" (or EOF if absent).
    for (_, l) in lines.by_ref() {
        if l.trim_start().starts_with("M  END") {
            break;
        }
    }

    Ok((builder.build(), metadata, coords))
}

/// Parse a MOL V2000 string into a `(Molecule, MolMetadata)` pair.
///
/// This is a convenience wrapper around [`parse_mol_with_coords`] that discards
/// the 2D coordinate data.
pub fn parse_mol(input: &str) -> Result<(Molecule, MolMetadata), MolParseError> {
    parse_mol_with_coords(input).map(|(mol, meta, _coords)| (mol, meta))
}

/// Parse all molecules from an SDF string, returning 2D coordinates.
///
/// Each entry contains the molecule, its metadata, and a `Vec<(x, y)>` of
/// 2D coordinates in atom-insertion order (the same order as `.atoms()`).
///
/// Stops and returns an error on the first parse failure.
#[allow(clippy::type_complexity)]
pub fn parse_sdf_with_coords(
    input: &str,
) -> Result<Vec<(Molecule, MolMetadata, Vec<(f64, f64)>)>, MolParseError> {
    // Re-use the SDF record splitter by borrowing its block-splitting logic,
    // but call parse_mol_with_coords on each block instead of parse_mol.
    let mut result = Vec::new();
    let mut remaining = input;
    loop {
        // Skip leading blank lines.
        while let Some(rest) = remaining
            .strip_prefix("\r\n")
            .or_else(|| remaining.strip_prefix('\n'))
        {
            remaining = rest;
        }
        if remaining.is_empty() {
            break;
        }

        // Find the $$$$ delimiter (line-by-line to avoid false matches inside data).
        let mut byte_offset = 0usize;
        let (end_byte, after_delim) = loop {
            let rest = &remaining[byte_offset..];
            match rest.find('\n') {
                Some(nl) => {
                    let line = rest[..nl].trim_end_matches('\r');
                    if line == "$$$$" {
                        break (byte_offset, &remaining[byte_offset + nl + 1..]);
                    }
                    byte_offset += nl + 1;
                }
                None => {
                    if rest.trim_end_matches('\r') == "$$$$" {
                        break (byte_offset, "");
                    }
                    break (remaining.len(), "");
                }
            }
        };

        let block = &remaining[..end_byte];
        remaining = after_delim;
        if block.trim().is_empty() {
            continue;
        }

        let (mol, meta, coords) = parse_mol_with_coords(block)?;
        result.push((mol, meta, coords));
    }
    Ok(result)
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
    write_mol_with_coords(mol, metadata, &[])
}

/// Serialize `mol` to a V2000 MOL block, using `coords` for atom positions.
///
/// `coords[i]` is the `(x, y)` position in Ångström for atom index `i`.
/// Atoms beyond `coords.len()` receive `(0.0, 0.0, 0.0)`.
pub fn write_mol_with_coords(
    mol: &Molecule,
    metadata: &MolMetadata,
    coords: &[(f64, f64)],
) -> String {
    let mut out = String::new();

    // Header lines 1–3
    out.push_str(&metadata.name);
    out.push('\n');
    out.push_str("  chematic\n");
    out.push_str(&metadata.comment);
    out.push('\n');

    // Counts line (line 4)
    let natoms = mol.atom_count();
    let nbonds = mol.bond_count();
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        natoms, nbonds
    ));

    // Atom block
    for (idx, atom) in mol.atoms() {
        let sym = atom.element.symbol();
        let charge_code = encode_charge(atom.charge);
        let (x, y) = coords.get(idx.0 as usize).copied().unwrap_or((0.0, 0.0));
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0{:>3}  0  0  0  0  0  0  0  0  0\n",
            x, y, 0.0_f64, sym, charge_code,
        ));
    }

    // Bond block
    for (_idx, bond) in mol.bonds() {
        let a1 = bond.atom1.0 + 1; // convert to 1-based
        let a2 = bond.atom2.0 + 1;
        let btype = match bond.order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::QuerySingleOrDouble => 5,
            BondOrder::QuerySingleOrAromatic => 6,
            BondOrder::QueryDoubleOrAromatic => 7,
            BondOrder::QueryAny | BondOrder::Zero => 8,
            BondOrder::Quadruple => 4,
        };
        out.push_str(&format!("{:>3}{:>3}{:>3}  0\n", a1, a2, btype));
    }

    // Terminator
    out.push_str("M  END\n");

    out
}

// ---------------------------------------------------------------------------
// SDF writer
// ---------------------------------------------------------------------------

/// Serialise one or more molecules to SDF format.
///
/// `records` — slice of `(molecule, metadata, coords)` tuples.
/// `coords` is optional; pass an empty slice to write zero coordinates.
/// Each molecule block is terminated with `$$$$`.
#[allow(clippy::type_complexity)]
pub fn write_sdf(records: &[(&Molecule, &MolMetadata, &[(f64, f64)])]) -> String {
    let mut out = String::new();
    for (mol, meta, coords) in records {
        out.push_str(&write_mol_with_coords(mol, meta, coords));
        out.push_str("$$$$\n");
    }
    out
}

/// Serialise one or more molecules to SDF format, appending per-atom partial
/// charges as an SD property `<PARTIAL_CHARGES>`.
///
/// `records` — slice of `(molecule, metadata, 2D-coords, charges)` tuples.
/// `charges[i]` is the partial charge for atom `i` (heavy atoms only).
/// Pass an empty charges slice to omit the property block.
///
/// Example SD block appended after `M  END`:
/// ```text
/// > <PARTIAL_CHARGES>
/// -0.2359 0.1076 -0.4500 0.1806
///
/// $$$$
/// ```
#[allow(clippy::type_complexity)]
pub fn write_sdf_with_charges(
    records: &[(&Molecule, &MolMetadata, &[(f64, f64)], &[f64])],
) -> String {
    let mut out = String::new();
    for (mol, meta, coords, charges) in records {
        out.push_str(&write_mol_with_coords(mol, meta, coords));
        if !charges.is_empty() {
            out.push_str("> <PARTIAL_CHARGES>\n");
            let vals: Vec<String> = charges.iter().map(|q| format!("{q:.4}")).collect();
            out.push_str(&vals.join(" "));
            out.push_str("\n\n");
        }
        out.push_str("$$$$\n");
    }
    out
}

/// Serialise a single molecule to one SDF record with arbitrary SD data fields.
///
/// Keys starting with `_` are treated as internal/computed properties and are
/// omitted from the SD block (e.g. `_Name` is written into the MOL header, not
/// as an SD field).  The record is terminated with `$$$$`.
pub fn write_sdf_record(
    mol: &Molecule,
    meta: &MolMetadata,
    coords: &[(f64, f64)],
    props: &std::collections::HashMap<String, String>,
) -> String {
    let mut out = write_mol_with_coords(mol, meta, coords);
    for (k, v) in props {
        if !k.starts_with('_') {
            out.push_str(&format!("> <{k}>\n{v}\n\n"));
        }
    }
    out.push_str("$$$$\n");
    out
}

/// Like [`write_sdf_record`] but emits a MOL V3000 (Extended Ctab) block
/// instead of V2000 — required for molecules with more than 999 atoms or
/// bonds, which don't fit V2000's fixed-width count fields.
pub fn write_sdf_record_v3000(
    mol: &Molecule,
    meta: &MolMetadata,
    coords: &[(f64, f64)],
    props: &std::collections::HashMap<String, String>,
) -> String {
    let mut out = crate::mol3000::write_mol_v3000(mol, meta, coords);
    for (k, v) in props {
        if !k.starts_with('_') {
            out.push_str(&format!("> <{k}>\n{v}\n\n"));
        }
    }
    out.push_str("$$$$\n");
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
    fn test_parse_query_bond_types_preserved() {
        let mol_str = "\
query_bonds
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
  1  2  5  0
  3  4  6  0
  5  6  7  0
  7  8  8  0
M  END
";
        let (mol, meta) = parse_mol(mol_str).expect("parse query bonds");
        let bonds: Vec<_> = mol.bonds().collect();
        assert_eq!(bonds[0].1.order, BondOrder::QuerySingleOrDouble);
        assert_eq!(bonds[1].1.order, BondOrder::QuerySingleOrAromatic);
        assert_eq!(bonds[2].1.order, BondOrder::QueryDoubleOrAromatic);
        assert_eq!(bonds[3].1.order, BondOrder::QueryAny);

        let written = write_mol(&mol, &meta);
        assert!(written.contains("  1  2  5  0"), "{written}");
        assert!(written.contains("  7  8  8  0"), "{written}");
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
    fn test_charged_atom_writes_implicit_h_in_smiles() {
        // MOL V2000 has no per-atom H-count field, so a charged atom read from
        // this format always has `hydrogen_count: None` (the format leaves
        // implicit H to be inferred, unlike bracket SMILES). Regression test
        // for a bug where the SMILES writer silently dropped these inferred
        // hydrogens on bracket atoms forced by charge/isotope/atom-map
        // (found via MRV oracle validation, see validation/mrv_io_parity_summary.json).
        let mol_str = "\
charged
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
M  END
";
        let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
        assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, None);
        assert_eq!(chematic_smiles::write(&mol), "[NH4+]");
        assert_eq!(chematic_smiles::canonical_smiles(&mol), "[NH4+]");
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

    #[test]
    fn test_molmetadata_builder() {
        let meta = MolMetadata::default()
            .with_name("aspirin")
            .with_comment("test molecule");
        assert_eq!(meta.name, "aspirin");
        assert_eq!(meta.comment, "test molecule");
    }

    #[test]
    fn test_molmetadata_with_name_roundtrip() {
        // Build a two-atom molecule and write it → name appears on MOL header line 1.
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();

        let meta = MolMetadata::default().with_name("acetic acid");
        let molblock = crate::write_mol(&mol, &meta);
        assert!(
            molblock.starts_with("acetic acid"),
            "MOL block must start with the molecule name"
        );
    }

    #[test]
    fn test_declared_max_atom_count_truncated_input_errors() {
        let bad = "\
max_atoms
  chematic

999  0  0  0  0  0  0  0  0  0  0 V2000
";
        assert!(matches!(parse_mol(bad), Err(MolParseError::UnexpectedEnd)));
    }

    #[test]
    fn test_declared_large_bond_count_truncated_input_errors() {
        let bad = "\
many_bonds
  chematic

  1 999  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
";
        assert!(matches!(parse_mol(bad), Err(MolParseError::UnexpectedEnd)));
    }

    #[test]
    fn test_write_sdf_record_props_roundtrip() {
        use crate::sdf::SdfRecordReader;
        use std::collections::HashMap;

        let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
        let coords = vec![(0.0, 0.0), (1.5, 0.0), (3.0, 0.0)];
        let mut props = HashMap::new();
        props.insert("Activity".to_string(), "7.2".to_string());
        props.insert("Source".to_string(), "test".to_string());
        props.insert("_Name".to_string(), "ethanol".to_string()); // internal, should be omitted

        let sdf = write_sdf_record(&mol, &meta, &coords, &props);

        // _Name must NOT appear as an SD field
        assert!(
            !sdf.contains("> <_Name>"),
            "internal prop leaked to SD block"
        );
        // Other props must appear
        assert!(sdf.contains("> <Activity>\n7.2"), "Activity missing");
        assert!(sdf.contains("> <Source>\ntest"), "Source missing");
        assert!(sdf.ends_with("$$$$\n"), "missing delimiter");

        // Parse back and verify
        let rec = SdfRecordReader::new(&sdf)
            .next()
            .expect("should have record")
            .expect("should parse");
        assert_eq!(rec.mol.atom_count(), mol.atom_count());
        assert_eq!(
            rec.properties.get("Activity").map(|s| s.as_str()),
            Some("7.2")
        );
        assert_eq!(
            rec.properties.get("Source").map(|s| s.as_str()),
            Some("test")
        );
        assert!(!rec.properties.contains_key("_Name"));
    }
}
