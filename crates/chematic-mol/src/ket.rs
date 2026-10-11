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
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
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
            KetError::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(f, "KET {resource} limit exceeded: {actual} > {limit}"),
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

/// Resource limits applied before and during KET interpretation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KetParseLimits {
    pub max_input_bytes: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
}

impl Default for KetParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_atoms: 10_000,
            max_bonds: 20_000,
        }
    }
}

// ─── Parser ──────────────────────────────────────────────────────────────────

/// Parse a KET (Ketcher JSON) string.
///
/// Returns `(molecule, coords_2d)` where `coords_2d` is the XY plane projection
/// of the KET `location` field (Z is discarded for 2D use).
/// For 3D use, call [`parse_ket_3d`] instead.
pub fn parse_ket(input: &str) -> Result<(Molecule, Vec<(f64, f64)>), KetError> {
    parse_ket_with_limits(input, &KetParseLimits::default())
}

/// Parse KET with explicit resource limits.
pub fn parse_ket_with_limits(
    input: &str,
    limits: &KetParseLimits,
) -> Result<(Molecule, Vec<(f64, f64)>), KetError> {
    let (mol, coords3) = parse_ket_3d_with_limits(input, limits)?;
    let coords2 = coords3.iter().map(|&(x, y, _z)| (x, y)).collect();
    Ok((mol, coords2))
}

/// Parse a KET string and return 3D coordinates `(x, y, z)`.
#[allow(clippy::type_complexity)]
pub fn parse_ket_3d(input: &str) -> Result<(Molecule, Vec<(f64, f64, f64)>), KetError> {
    parse_ket_3d_with_limits(input, &KetParseLimits::default())
}

fn ket_molecule_node(root: &serde_json::Value) -> Result<&serde_json::Value, KetError> {
    if root.get("root").is_none() {
        return Ok(root);
    }
    let reference = root
        .get("root")
        .and_then(|value| value.get("nodes"))
        .and_then(serde_json::Value::as_array)
        .and_then(|nodes| nodes.first())
        .and_then(|node| node.get("$ref"))
        .and_then(serde_json::Value::as_str)
        .ok_or(KetError::MissingField("root.nodes[0].$ref"))?;
    root.get(reference)
        .ok_or(KetError::MissingField("molecule node referenced by $ref"))
}

fn ket_records<'a>(
    node: &'a serde_json::Value,
    limits: &KetParseLimits,
) -> Result<(&'a [serde_json::Value], &'a [serde_json::Value]), KetError> {
    let atoms = node
        .get("atoms")
        .and_then(serde_json::Value::as_array)
        .ok_or(KetError::MissingField("atoms"))?;
    let bonds = node
        .get("bonds")
        .and_then(serde_json::Value::as_array)
        .map(Vec::as_slice)
        .unwrap_or(&[]);
    for (resource, actual, limit) in [
        ("atoms", atoms.len(), limits.max_atoms),
        ("bonds", bonds.len(), limits.max_bonds),
    ] {
        if actual > limit {
            return Err(KetError::ResourceLimit {
                resource,
                actual,
                limit,
            });
        }
    }
    Ok((atoms, bonds))
}

struct KetAtom {
    atom: Atom,
    coord: (f64, f64, f64),
}

fn parse_ket_atom(value: &serde_json::Value) -> Result<KetAtom, KetError> {
    let label = value
        .get("label")
        .and_then(serde_json::Value::as_str)
        .ok_or(KetError::MissingField("atom.label"))?;
    let placeholder = label == "*"
        || label == "R"
        || (label.starts_with('R')
            && label
                .chars()
                .nth(1)
                .is_some_and(|character| character.is_ascii_digit() || character == '#'));
    let element = if placeholder {
        Element::from_symbol("C").expect("carbon is a supported element")
    } else {
        Element::from_symbol(label).ok_or_else(|| KetError::UnknownElement(label.to_owned()))?
    };
    let mut atom = Atom::new(element);
    atom.charge = value
        .get("charge")
        .and_then(serde_json::Value::as_i64)
        .unwrap_or(0) as i8;
    atom.isotope = value
        .get("isotope")
        .and_then(serde_json::Value::as_u64)
        .filter(|&isotope| isotope != 0 && isotope <= u16::MAX as u64)
        .map(|isotope| isotope as u16);
    atom.hydrogen_count = value
        .get("explicitHydrogens")
        .and_then(serde_json::Value::as_i64)
        .filter(|&count| count >= 0)
        .map(|count| count as u8);
    atom.wildcard = label == "*";
    let location = value.get("location").and_then(serde_json::Value::as_array);
    let coordinate = |index| {
        location
            .and_then(|values| values.get(index))
            .and_then(serde_json::Value::as_f64)
            .unwrap_or(0.0)
    };
    Ok(KetAtom {
        atom,
        coord: (coordinate(0), coordinate(1), coordinate(2)),
    })
}

fn ket_bond_order(value: &serde_json::Value) -> BondOrder {
    let bond_type = value
        .get("type")
        .and_then(serde_json::Value::as_u64)
        .unwrap_or(1);
    let stereo = value
        .get("stereo")
        .and_then(serde_json::Value::as_u64)
        .unwrap_or(0);
    match (bond_type, stereo) {
        (1, 1) => BondOrder::Up,
        (1, 6) => BondOrder::Down,
        (2, _) => BondOrder::Double,
        (3, _) => BondOrder::Triple,
        (4, _) => BondOrder::Aromatic,
        _ => BondOrder::Single,
    }
}

fn add_ket_bonds(
    bonds: &[serde_json::Value],
    count: usize,
    builder: &mut MoleculeBuilder,
) -> Result<(), KetError> {
    use chematic_core::AtomIdx;
    for (index, value) in bonds.iter().enumerate() {
        let atoms = value
            .get("atoms")
            .and_then(serde_json::Value::as_array)
            .ok_or(KetError::MissingField("bond.atoms"))?;
        let endpoint = |position: usize, field: &'static str| {
            atoms
                .get(position)
                .and_then(serde_json::Value::as_u64)
                .map(|atom| atom as usize)
                .ok_or(KetError::MissingField(field))
        };
        let source = endpoint(0, "bond.atoms[0]")?;
        let target = endpoint(1, "bond.atoms[1]")?;
        for atom in [source, target] {
            if atom >= count {
                return Err(KetError::InvalidAtomIndex {
                    bond_idx: index,
                    atom_idx: atom,
                    natoms: count,
                });
            }
        }
        let _ = builder.add_bond(
            AtomIdx(source as u32),
            AtomIdx(target as u32),
            ket_bond_order(value),
        );
    }
    Ok(())
}

/// Parse KET with explicit resource limits and return 3D coordinates.
#[allow(clippy::type_complexity)]
pub fn parse_ket_3d_with_limits(
    input: &str,
    limits: &KetParseLimits,
) -> Result<(Molecule, Vec<(f64, f64, f64)>), KetError> {
    if input.len() > limits.max_input_bytes {
        return Err(KetError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let root: serde_json::Value =
        serde_json::from_str(input).map_err(|e| KetError::InvalidJson(e.to_string()))?;
    let node = ket_molecule_node(&root)?;
    let (atoms, bonds) = ket_records(node, limits)?;
    let mut builder = MoleculeBuilder::new();
    let mut coords = Vec::with_capacity(atoms.len());
    for value in atoms {
        let parsed = parse_ket_atom(value)?;
        builder.add_atom(parsed.atom);
        coords.push(parsed.coord);
    }
    add_ket_bonds(bonds, atoms.len(), &mut builder)?;
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

    #[test]
    fn parse_ket_limits_are_typed_and_shared_by_2d_and_3d_paths() {
        let ket = r#"{"atoms":[{"label":"C","location":[0,0,0]}],"bonds":[]}"#;
        assert!(matches!(
            parse_ket_with_limits(
                ket,
                &KetParseLimits {
                    max_input_bytes: 4,
                    ..Default::default()
                }
            ),
            Err(KetError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));

        let many_atoms = r#"{"atoms":[{"label":"C"},{"label":"C"}],"bonds":[]}"#;
        let limits = KetParseLimits {
            max_atoms: 1,
            ..Default::default()
        };
        assert!(matches!(
            parse_ket_3d_with_limits(many_atoms, &limits),
            Err(KetError::ResourceLimit {
                resource: "atoms",
                actual: 2,
                limit: 1
            })
        ));

        let many_bonds =
            r#"{"atoms":[{"label":"C"},{"label":"C"}],"bonds":[{"atoms":[0,1]},{"atoms":[0,1]}]}"#;
        assert!(matches!(
            parse_ket_3d_with_limits(
                many_bonds,
                &KetParseLimits {
                    max_bonds: 1,
                    ..Default::default()
                }
            ),
            Err(KetError::ResourceLimit {
                resource: "bonds",
                actual: 2,
                limit: 1
            })
        ));
    }
}
