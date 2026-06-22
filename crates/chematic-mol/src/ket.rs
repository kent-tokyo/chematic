//! KET (Ketcher JSON) format parser and writer.
//!
//! KET is the native JSON format used by [Ketcher](https://github.com/epam/ketcher),
//! an open-source web-based chemical structure editor.  It is also accepted by EPAM's
//! Indigo toolkit.
//!
//! ## Format overview
//!
//! ```json
//! {
//!   "root": { "nodes": [{ "$ref": "mol0" }] },
//!   "mol0": {
//!     "type": "molecule",
//!     "atoms": [
//!       { "label": "C", "location": [0.0, 0.0, 0.0], "charge": 0 }
//!     ],
//!     "bonds": [
//!       { "type": 1, "atoms": [0, 1] }
//!     ]
//!   }
//! }
//! ```
//!
//! A simpler flat variant (no `root`/`molN` nesting) is also accepted:
//! ```json
//! { "atoms": [...], "bonds": [...] }
//! ```
//!
//! ## Bond type codes
//! 1=single, 2=double, 3=triple, 4=aromatic, 5=single-or-double,
//! 6=single-or-aromatic, 7=double-or-aromatic, 8=any
//!
//! ## Stereo bond stereo codes
//! 1=up (wedge), 6=down (hash/dash), 4=either (wavy)

#![forbid(unsafe_code)]

use chematic_core::{Atom, BondOrder, Element, Molecule, MoleculeBuilder};

// ─── Error ───────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KetError {
    InvalidJson(String),
    UnknownElement(String),
    InvalidAtomIndex {
        bond_idx: usize,
        atom_idx: usize,
        natoms: usize,
    },
    MissingField(&'static str),
}

impl std::fmt::Display for KetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KetError::InvalidJson(s) => write!(f, "KET: invalid JSON: {s}"),
            KetError::UnknownElement(sym) => write!(f, "KET: unknown element '{sym}'"),
            KetError::InvalidAtomIndex {
                bond_idx,
                atom_idx,
                natoms,
            } => write!(
                f,
                "KET: bond {bond_idx} references atom {atom_idx} (only {natoms} atoms)"
            ),
            KetError::MissingField(fld) => write!(f, "KET: missing required field '{fld}'"),
        }
    }
}
impl std::error::Error for KetError {}

// ─── Parser ──────────────────────────────────────────────────────────────────

/// Parse a KET (Ketcher JSON) string.
///
/// Returns `(molecule, coords_2d)` where `coords_2d` is the XY plane projection
/// of the KET `location` field (Z is discarded for 2D use).
/// For 3D use, call [`parse_ket_3d`] instead.
pub fn parse_ket(input: &str) -> Result<(Molecule, Vec<(f64, f64)>), KetError> {
    let (mol, coords3) = parse_ket_3d(input)?;
    let coords2 = coords3.iter().map(|&(x, y, _z)| (x, y)).collect();
    Ok((mol, coords2))
}

/// Parse a KET string and return 3D coordinates `(x, y, z)`.
#[allow(clippy::type_complexity)]
pub fn parse_ket_3d(input: &str) -> Result<(Molecule, Vec<(f64, f64, f64)>), KetError> {
    let v: serde_json::Value =
        serde_json::from_str(input).map_err(|e| KetError::InvalidJson(e.to_string()))?;

    // Locate the molecule node.  Two supported layouts:
    // 1. Flat:   { "atoms": [...], "bonds": [...] }
    // 2. Nested: { "root": { "nodes": [{"$ref":"mol0"}] }, "mol0": { "atoms":[...], "bonds":[...] } }
    let mol_node: &serde_json::Value = if v.get("atoms").is_some() {
        &v
    } else if let Some(root) = v.get("root") {
        // Find first $ref under root.nodes
        let mol_ref = root
            .get("nodes")
            .and_then(|n| n.as_array())
            .and_then(|arr| arr.first())
            .and_then(|n| n.get("$ref"))
            .and_then(|r| r.as_str())
            .ok_or(KetError::MissingField("root.nodes[0].$ref"))?;
        v.get(mol_ref)
            .ok_or(KetError::MissingField("molecule node referenced by $ref"))?
    } else {
        return Err(KetError::MissingField("atoms"));
    };

    let atoms_json = mol_node
        .get("atoms")
        .and_then(|a| a.as_array())
        .ok_or(KetError::MissingField("atoms"))?;

    // Guard against memory-DoS from enormous KET payloads.
    const MAX_ATOMS: usize = 10_000;
    const MAX_BONDS: usize = 20_000;
    if atoms_json.len() > MAX_ATOMS {
        return Err(KetError::InvalidJson(format!(
            "KET file exceeds atom limit ({} > {MAX_ATOMS})",
            atoms_json.len()
        )));
    }

    let empty_bonds = vec![];
    let bonds_json = mol_node
        .get("bonds")
        .and_then(|b| b.as_array())
        .unwrap_or(&empty_bonds); // bonds is optional (single-atom molecules)

    if bonds_json.len() > MAX_BONDS {
        return Err(KetError::InvalidJson(format!(
            "KET file exceeds bond limit ({} > {MAX_BONDS})",
            bonds_json.len()
        )));
    }

    // ── Atoms ─────────────────────────────────────────────────────────────────
    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64, f64)> = Vec::with_capacity(atoms_json.len());

    for atom_v in atoms_json {
        let label = atom_v
            .get("label")
            .and_then(|l| l.as_str())
            .ok_or(KetError::MissingField("atom.label"))?;

        // Handle special labels:
        // "*" = any-atom wildcard
        // "R", "R#", "R1"–"R9"… = R-group placeholders (Ketcher uses "R#" or "R1" etc.)
        // NOTE: must NOT catch real elements that start with 'R' (Ru, Rh, Re, Rn, Ra, Rb, Rf, Rg).
        let is_rgroup = label == "*"
            || label == "R"
            || (label.starts_with('R')
                && label
                    .chars()
                    .nth(1)
                    .is_some_and(|c| c.is_ascii_digit() || c == '#'));
        let element = if is_rgroup {
            Element::from_symbol("C").unwrap() // fallback wildcard → carbon placeholder
        } else {
            Element::from_symbol(label)
                .ok_or_else(|| KetError::UnknownElement(label.to_string()))?
        };

        let charge = atom_v.get("charge").and_then(|c| c.as_i64()).unwrap_or(0) as i8;

        let isotope = atom_v
            .get("isotope")
            .and_then(|i| i.as_u64())
            .filter(|&i| i != 0 && i <= u16::MAX as u64) // guard against silent u64→u16 truncation
            .map(|i| i as u16);

        let explicit_h = atom_v
            .get("explicitHydrogens")
            .and_then(|h| h.as_i64())
            .filter(|&h| h >= 0)
            .map(|h| h as u8);

        let wildcard = label == "*";

        let mut atom = Atom::new(element);
        atom.charge = charge;
        atom.isotope = isotope;
        atom.hydrogen_count = explicit_h;
        atom.wildcard = wildcard;
        builder.add_atom(atom);

        let loc = atom_v.get("location").and_then(|l| l.as_array());
        let (x, y, z) = if let Some(loc) = loc {
            let x = loc.first().and_then(|v| v.as_f64()).unwrap_or(0.0);
            let y = loc.get(1).and_then(|v| v.as_f64()).unwrap_or(0.0);
            let z = loc.get(2).and_then(|v| v.as_f64()).unwrap_or(0.0);
            (x, y, z)
        } else {
            (0.0, 0.0, 0.0)
        };
        coords.push((x, y, z));
    }

    let natoms = atoms_json.len();

    // ── Bonds ─────────────────────────────────────────────────────────────────
    for (bidx, bond_v) in bonds_json.iter().enumerate() {
        let bond_type = bond_v.get("type").and_then(|t| t.as_u64()).unwrap_or(1) as u8;

        let atom_refs = bond_v
            .get("atoms")
            .and_then(|a| a.as_array())
            .ok_or(KetError::MissingField("bond.atoms"))?;

        let a = atom_refs
            .first()
            .and_then(|v| v.as_u64())
            .ok_or(KetError::MissingField("bond.atoms[0]"))? as usize;
        let b = atom_refs
            .get(1)
            .and_then(|v| v.as_u64())
            .ok_or(KetError::MissingField("bond.atoms[1]"))? as usize;

        if a >= natoms {
            return Err(KetError::InvalidAtomIndex {
                bond_idx: bidx,
                atom_idx: a,
                natoms,
            });
        }
        if b >= natoms {
            return Err(KetError::InvalidAtomIndex {
                bond_idx: bidx,
                atom_idx: b,
                natoms,
            });
        }

        // Stereo from bond.stereo field (optional)
        let stereo = bond_v.get("stereo").and_then(|s| s.as_u64()).unwrap_or(0);

        let order = match bond_type {
            1 => match stereo {
                1 => BondOrder::Up,
                6 => BondOrder::Down,
                _ => BondOrder::Single,
            },
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Aromatic,
            _ => BondOrder::Single, // query / any bond types → single
        };

        use chematic_core::AtomIdx;
        let _ = builder.add_bond(AtomIdx(a as u32), AtomIdx(b as u32), order);
    }

    Ok((builder.build(), coords))
}

// ─── Writer ───────────────────────────────────────────────────────────────────

/// Write a molecule to KET (Ketcher JSON) format with 2D coordinates.
///
/// Produces a compact KET string using the nested `root`/`mol0` layout.
pub fn write_ket(mol: &Molecule, coords: &[(f64, f64)]) -> String {
    let coords3: Vec<(f64, f64, f64)> = coords.iter().map(|&(x, y)| (x, y, 0.0)).collect();
    write_ket_3d(mol, &coords3)
}

/// Write a molecule to KET format with 3D coordinates.
pub fn write_ket_3d(mol: &Molecule, coords: &[(f64, f64, f64)]) -> String {
    use serde_json::{Value, json};

    let atoms: Vec<Value> = (0..mol.atom_count())
        .map(|i| {
            use chematic_core::AtomIdx;
            let atom = mol.atom(AtomIdx(i as u32));
            let (x, y, z) = coords.get(i).copied().unwrap_or((0.0, 0.0, 0.0));
            let mut a = json!({
                "label": atom.element.symbol(),
                "location": [x, y, z],
            });
            if atom.charge != 0 {
                a["charge"] = json!(atom.charge);
            }
            if let Some(iso) = atom.isotope {
                a["isotope"] = json!(iso);
            }
            a
        })
        .collect();

    let bonds: Vec<Value> = (0..mol.bond_count())
        .map(|i| {
            use chematic_core::BondIdx;
            let bond = mol.bond(BondIdx(i as u32));
            let (bond_type, stereo) = match bond.order {
                BondOrder::Single => (1u8, 0u8),
                BondOrder::Up => (1, 1),
                BondOrder::Down => (1, 6),
                BondOrder::Double => (2, 0),
                BondOrder::Triple => (3, 0),
                BondOrder::Aromatic => (4, 0),
                _ => (1, 0),
            };
            let mut b = json!({
                "type": bond_type,
                "atoms": [bond.atom1.0, bond.atom2.0],
            });
            if stereo != 0 {
                b["stereo"] = json!(stereo);
            }
            b
        })
        .collect();

    let mol_obj = json!({
        "type": "molecule",
        "atoms": atoms,
        "bonds": bonds,
    });

    let root = json!({
        "root": {
            "nodes": [{ "$ref": "mol0" }]
        },
        "mol0": mol_obj
    });

    serde_json::to_string(&root).unwrap_or_default()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::BondOrder;

    const ETHANOL_KET: &str = r#"
    {
      "root": { "nodes": [{ "$ref": "mol0" }] },
      "mol0": {
        "type": "molecule",
        "atoms": [
          { "label": "C", "location": [0.0, 0.0, 0.0] },
          { "label": "C", "location": [1.5, 0.0, 0.0] },
          { "label": "O", "location": [3.0, 0.0, 0.0] }
        ],
        "bonds": [
          { "type": 1, "atoms": [0, 1] },
          { "type": 1, "atoms": [1, 2] }
        ]
      }
    }"#;

    const FLAT_BENZENE: &str = r#"
    {
      "atoms": [
        {"label": "C", "location": [0.0, 1.0, 0.0]},
        {"label": "C", "location": [0.866, 0.5, 0.0]},
        {"label": "C", "location": [0.866, -0.5, 0.0]},
        {"label": "C", "location": [0.0, -1.0, 0.0]},
        {"label": "C", "location": [-0.866, -0.5, 0.0]},
        {"label": "C", "location": [-0.866, 0.5, 0.0]}
      ],
      "bonds": [
        {"type": 4, "atoms": [0, 1]},
        {"type": 4, "atoms": [1, 2]},
        {"type": 4, "atoms": [2, 3]},
        {"type": 4, "atoms": [3, 4]},
        {"type": 4, "atoms": [4, 5]},
        {"type": 4, "atoms": [5, 0]}
      ]
    }"#;

    #[test]
    fn parse_ket_ethanol_nested() {
        let (mol, coords) = parse_ket(ETHANOL_KET).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        assert_eq!(coords.len(), 3);
        assert_eq!(mol.atom(chematic_core::AtomIdx(2)).element.symbol(), "O");
    }

    #[test]
    fn parse_ket_benzene_flat() {
        let (mol, _) = parse_ket(FLAT_BENZENE).unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        // All bonds should be aromatic
        for i in 0..mol.bond_count() {
            assert_eq!(
                mol.bond(chematic_core::BondIdx(i as u32)).order,
                BondOrder::Aromatic
            );
        }
    }

    #[test]
    fn ket_round_trip_ethanol() {
        let (mol, coords) = parse_ket(ETHANOL_KET).unwrap();
        let out = write_ket(&mol, &coords);
        let (mol2, coords2) = parse_ket(&out).unwrap();
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert_eq!(mol.bond_count(), mol2.bond_count());
        for (i, ((x1, y1), (x2, y2))) in coords.iter().zip(coords2.iter()).enumerate() {
            assert!((x1 - x2).abs() < 1e-9, "coord x mismatch at {i}");
            assert!((y1 - y2).abs() < 1e-9, "coord y mismatch at {i}");
        }
    }

    #[test]
    fn parse_ket_charged_atom() {
        let ket = r#"{"atoms":[{"label":"N","location":[0,0,0],"charge":1}],"bonds":[]}"#;
        let (mol, _) = parse_ket(ket).unwrap();
        assert_eq!(mol.atom(chematic_core::AtomIdx(0)).charge, 1);
    }

    #[test]
    fn parse_ket_stereo_bonds() {
        let ket = r#"{"atoms":[
            {"label":"C","location":[0,0,0]},
            {"label":"F","location":[1,0,0]},
            {"label":"Cl","location":[0,1,0]}
        ],"bonds":[
            {"type":1,"atoms":[0,1],"stereo":1},
            {"type":1,"atoms":[0,2],"stereo":6}
        ]}"#;
        let (mol, _) = parse_ket(ket).unwrap();
        assert_eq!(mol.bond(chematic_core::BondIdx(0)).order, BondOrder::Up);
        assert_eq!(mol.bond(chematic_core::BondIdx(1)).order, BondOrder::Down);
    }

    #[test]
    fn write_ket_produces_valid_json() {
        use chematic_smiles::parse;
        let mol = parse("CCO").unwrap();
        let coords = vec![(0.0, 0.0), (1.5, 0.0), (3.0, 0.0)];
        let out = write_ket(&mol, &coords);
        // Must parse as JSON
        let _v: serde_json::Value = serde_json::from_str(&out).unwrap();
        // Must round-trip
        let (mol2, _) = parse_ket(&out).unwrap();
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert_eq!(mol.bond_count(), mol2.bond_count());
    }

    #[test]
    fn parse_ket_unknown_element_errors() {
        let ket = r#"{"atoms":[{"label":"Xx","location":[0,0,0]}],"bonds":[]}"#;
        assert!(matches!(parse_ket(ket), Err(KetError::UnknownElement(_))));
    }

    #[test]
    fn parse_ket_bad_bond_index_errors() {
        let ket =
            r#"{"atoms":[{"label":"C","location":[0,0,0]}],"bonds":[{"type":1,"atoms":[0,99]}]}"#;
        assert!(matches!(
            parse_ket(ket),
            Err(KetError::InvalidAtomIndex { .. })
        ));
    }
}
