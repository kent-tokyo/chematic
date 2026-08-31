//! ChemicalJSON (.cjson) format parser and writer.
//!
//! ChemicalJSON is an open JSON format used by Avogadro, Avogadro2, and the
//! MolSSI Open Chemistry toolkit.  It stores atoms, bonds, and optionally
//! 3D coordinates, formal charges, and isotope labels.
//!
//! ## Format overview
//!
//! ```json
//! {
//!   "chemicalJson": 1,
//!   "atoms": {
//!     "coords": { "3d": [x0,y0,z0, x1,y1,z1, ...] },
//!     "elements": { "number": [6, 7, 8, ...] },
//!     "formalCharges": [0, 1, -1, ...],
//!     "isotopes": [0, 13, 0, ...]
//!   },
//!   "bonds": {
//!     "connections": { "index": [0,1, 1,2, ...] },
//!     "order": [1, 2, 1.5, ...]
//!   }
//! }
//! ```
//!
//! ## Bond order values
//! `1` = single, `2` = double, `3` = triple, `1.5` or `4` = aromatic.
//!
//! ## References
//! Open Chemistry wiki: <https://wiki.openchemistry.org/Chemical_JSON>

#![forbid(unsafe_code)]

use chematic_core::{Atom, BondOrder, Element, Molecule, MoleculeBuilder};

/// Resource limits for ChemicalJSON parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CjsonParseLimits {
    pub max_input_bytes: usize,
    pub max_json_depth: usize,
    pub max_array_items: usize,
    pub max_string_bytes: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
}

impl Default for CjsonParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_json_depth: 128,
            max_array_items: 10_000_000,
            max_string_bytes: 16 * 1024 * 1024,
            max_atoms: 1_000_000,
            max_bonds: 2_000_000,
        }
    }
}

// ─── Error ───────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CjsonError {
    InvalidJson(String),
    UnknownAtomicNumber(u64),
    InvalidBondIndex {
        pair: usize,
        atom_idx: u64,
        n_atoms: usize,
    },
    MissingField(&'static str),
    OddBondIndexList,
    InvalidCoordCount {
        expected: usize,
        got: usize,
    },
    InvalidValue {
        field: &'static str,
    },
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl std::fmt::Display for CjsonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CjsonError::InvalidJson(s) => write!(f, "CJSON: invalid JSON: {s}"),
            CjsonError::UnknownAtomicNumber(n) => {
                write!(f, "CJSON: unknown atomic number {n}")
            }
            CjsonError::InvalidBondIndex {
                pair,
                atom_idx,
                n_atoms,
            } => write!(
                f,
                "CJSON: bond pair {pair} references atom {atom_idx} (only {n_atoms} atoms)"
            ),
            CjsonError::MissingField(fld) => {
                write!(f, "CJSON: missing required field '{fld}'")
            }
            CjsonError::OddBondIndexList => {
                write!(
                    f,
                    "CJSON: bonds.connections.index must have an even length (pairs)"
                )
            }
            CjsonError::InvalidCoordCount { expected, got } => write!(
                f,
                "CJSON: coords.3d length {got} is not divisible by 3 \
                 or does not match atom count (expected {expected} × 3)"
            ),
            CjsonError::InvalidValue { field } => {
                write!(f, "CJSON: invalid value in '{field}'")
            }
            CjsonError::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "CJSON: {resource} has size {actual}, exceeding limit {limit}"
            ),
        }
    }
}

impl std::error::Error for CjsonError {}

// ─── Parser ──────────────────────────────────────────────────────────────────

/// Parse a ChemicalJSON string.
///
/// Returns `(molecule, coords)` where `coords` is a `Vec<(x, y, z)>` in Å.
/// `coords` is empty when the file has no `atoms.coords.3d` field.
#[allow(clippy::type_complexity)]
pub fn parse_cjson(input: &str) -> Result<(Molecule, Vec<(f64, f64, f64)>), CjsonError> {
    parse_cjson_with_limits(input, &CjsonParseLimits::default())
}

/// Parse ChemicalJSON with explicit resource limits.
#[allow(clippy::type_complexity)]
pub fn parse_cjson_with_limits(
    input: &str,
    limits: &CjsonParseLimits,
) -> Result<(Molecule, Vec<(f64, f64, f64)>), CjsonError> {
    if input.len() > limits.max_input_bytes {
        return Err(CjsonError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let v: serde_json::Value =
        serde_json::from_str(input).map_err(|e| CjsonError::InvalidJson(e.to_string()))?;
    check_json_limits(&v, 0, limits)?;

    let atoms_obj = v.get("atoms").ok_or(CjsonError::MissingField("atoms"))?;

    // ── Atomic numbers ───────────────────────────────────────────────────────
    let numbers_value = atoms_obj
        .get("elements")
        .and_then(|e| e.get("number"))
        .and_then(|n| n.as_array())
        .ok_or(CjsonError::MissingField("atoms.elements.number"))?;
    if numbers_value.len() > limits.max_atoms {
        return Err(CjsonError::ResourceLimit {
            resource: "atom records",
            actual: numbers_value.len(),
            limit: limits.max_atoms,
        });
    }
    let numbers: Vec<u64> = numbers_value
        .iter()
        .map(|v| {
            v.as_u64().ok_or(CjsonError::InvalidValue {
                field: "atoms.elements.number",
            })
        })
        .collect::<Result<_, _>>()?;
    let n_atoms = numbers.len();

    // ── Optional: formal charges ─────────────────────────────────────────────
    let formal_charges: Vec<i8> =
        if let Some(arr) = atoms_obj.get("formalCharges").and_then(|fc| fc.as_array()) {
            arr.iter()
                .map(|v| {
                    let value = v.as_i64().ok_or(CjsonError::InvalidValue {
                        field: "atoms.formalCharges",
                    })?;
                    i8::try_from(value).map_err(|_| CjsonError::InvalidValue {
                        field: "atoms.formalCharges",
                    })
                })
                .collect::<Result<Vec<_>, _>>()?
        } else {
            vec![0i8; n_atoms]
        };

    // ── Optional: isotopes (0 = natural abundance) ───────────────────────────
    let isotopes: Vec<Option<u16>> =
        if let Some(arr) = atoms_obj.get("isotopes").and_then(|iso| iso.as_array()) {
            arr.iter()
                .map(|v| {
                    let n = v.as_u64().ok_or(CjsonError::InvalidValue {
                        field: "atoms.isotopes",
                    })?;
                    if n == 0 {
                        Ok(None)
                    } else {
                        u16::try_from(n)
                            .map(Some)
                            .map_err(|_| CjsonError::InvalidValue {
                                field: "atoms.isotopes",
                            })
                    }
                })
                .collect::<Result<Vec<_>, _>>()?
        } else {
            vec![None; n_atoms]
        };

    // ── 3D coordinates (optional) ────────────────────────────────────────────
    let coords: Vec<(f64, f64, f64)> = {
        let flat: Option<Result<Vec<f64>, CjsonError>> = atoms_obj
            .get("coords")
            .and_then(|c| c.get("3d"))
            .and_then(|c| c.as_array())
            .map(|arr| {
                arr.iter()
                    .map(|v| {
                        v.as_f64().ok_or(CjsonError::InvalidValue {
                            field: "atoms.coords.3d",
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
            });

        match flat.transpose()? {
            None => Vec::new(),
            Some(flat) => {
                if flat.len() % 3 != 0 || flat.len() / 3 != n_atoms {
                    return Err(CjsonError::InvalidCoordCount {
                        expected: n_atoms,
                        got: flat.len(),
                    });
                }
                flat.as_chunks::<3>()
                    .0
                    .iter()
                    .map(|c| (c[0], c[1], c[2]))
                    .collect()
            }
        }
    };

    // ── Build atoms ──────────────────────────────────────────────────────────
    let mut builder = MoleculeBuilder::new();
    for (i, &an) in numbers.iter().enumerate() {
        let element = u8::try_from(an)
            .ok()
            .and_then(Element::from_atomic_number)
            .ok_or(CjsonError::UnknownAtomicNumber(an))?;
        let mut atom = Atom::new(element);
        atom.charge = *formal_charges.get(i).unwrap_or(&0);
        atom.isotope = isotopes.get(i).copied().flatten();
        builder.add_atom(atom);
    }

    // ── Bonds ────────────────────────────────────────────────────────────────
    if let Some(bonds_obj) = v.get("bonds") {
        let index_value = bonds_obj
            .get("connections")
            .and_then(|c| c.get("index"))
            .and_then(|i| i.as_array())
            .ok_or(CjsonError::MissingField("bonds.connections.index"))?;
        if index_value.len() / 2 > limits.max_bonds {
            return Err(CjsonError::ResourceLimit {
                resource: "bond records",
                actual: index_value.len() / 2,
                limit: limits.max_bonds,
            });
        }
        let index: Vec<u64> = index_value
            .iter()
            .map(|v| {
                v.as_u64().ok_or(CjsonError::InvalidValue {
                    field: "bonds.connections.index",
                })
            })
            .collect::<Result<_, _>>()?;

        if !index.len().is_multiple_of(2) {
            return Err(CjsonError::OddBondIndexList);
        }

        let orders: Vec<f64> = if let Some(arr) = bonds_obj.get("order").and_then(|o| o.as_array())
        {
            arr.iter()
                .map(|v| {
                    v.as_f64().ok_or(CjsonError::InvalidValue {
                        field: "bonds.order",
                    })
                })
                .collect::<Result<Vec<_>, _>>()?
        } else {
            vec![1.0; index.len() / 2]
        };

        for (pair_idx, chunk) in index.as_chunks::<2>().0.iter().enumerate() {
            let a = chunk[0];
            let b = chunk[1];
            let a_idx = usize::try_from(a).map_err(|_| CjsonError::InvalidBondIndex {
                pair: pair_idx,
                atom_idx: a,
                n_atoms,
            })?;
            let b_idx = usize::try_from(b).map_err(|_| CjsonError::InvalidBondIndex {
                pair: pair_idx,
                atom_idx: b,
                n_atoms,
            })?;
            if a_idx >= n_atoms {
                return Err(CjsonError::InvalidBondIndex {
                    pair: pair_idx,
                    atom_idx: a,
                    n_atoms,
                });
            }
            if b_idx >= n_atoms {
                return Err(CjsonError::InvalidBondIndex {
                    pair: pair_idx,
                    atom_idx: b,
                    n_atoms,
                });
            }
            let order = float_to_bond_order(*orders.get(pair_idx).unwrap_or(&1.0));
            let _ = builder.add_bond(
                chematic_core::AtomIdx(a_idx as u32),
                chematic_core::AtomIdx(b_idx as u32),
                order,
            );
        }
    }

    Ok((builder.build(), coords))
}

fn check_json_limits(
    value: &serde_json::Value,
    depth: usize,
    limits: &CjsonParseLimits,
) -> Result<(), CjsonError> {
    if depth > limits.max_json_depth {
        return Err(CjsonError::ResourceLimit {
            resource: "JSON depth",
            actual: depth,
            limit: limits.max_json_depth,
        });
    }
    match value {
        serde_json::Value::String(s) if s.len() > limits.max_string_bytes => {
            return Err(CjsonError::ResourceLimit {
                resource: "string bytes",
                actual: s.len(),
                limit: limits.max_string_bytes,
            });
        }
        serde_json::Value::Array(values) => {
            if values.len() > limits.max_array_items {
                return Err(CjsonError::ResourceLimit {
                    resource: "array items",
                    actual: values.len(),
                    limit: limits.max_array_items,
                });
            }
            for child in values {
                check_json_limits(child, depth + 1, limits)?;
            }
        }
        serde_json::Value::Object(object) => {
            for (key, child) in object {
                if key.len() > limits.max_string_bytes {
                    return Err(CjsonError::ResourceLimit {
                        resource: "object key bytes",
                        actual: key.len(),
                        limit: limits.max_string_bytes,
                    });
                }
                check_json_limits(child, depth + 1, limits)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Map a float bond-order value to [`BondOrder`].
fn float_to_bond_order(v: f64) -> BondOrder {
    // 1.5 or 4 → Aromatic (different implementations use different values)
    if (v - 1.5).abs() < 0.1 || (v - 4.0).abs() < 0.1 {
        return BondOrder::Aromatic;
    }
    match v.round() as i64 {
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        _ => BondOrder::Single,
    }
}

// ─── Writer ──────────────────────────────────────────────────────────────────

/// Write a [`Molecule`] and optional 3D coordinates as a ChemicalJSON string.
///
/// `coords` may be empty; if so the `atoms.coords` section is omitted and
/// Avogadro will open the file without 3D geometry.
pub fn write_cjson(mol: &Molecule, coords: &[(f64, f64, f64)]) -> String {
    use serde_json::{Value, json};

    let n = mol.atom_count();

    // Atomic numbers
    let numbers: Vec<Value> = (0..n)
        .map(|i| {
            json!(
                mol.atom(chematic_core::AtomIdx(i as u32))
                    .element
                    .atomic_number()
            )
        })
        .collect();

    // Formal charges (omit section if all zero)
    let charges: Vec<i8> = (0..n)
        .map(|i| mol.atom(chematic_core::AtomIdx(i as u32)).charge)
        .collect();
    let has_charges = charges.iter().any(|&c| c != 0);

    // Isotopes (omit section if none set)
    let isotopes: Vec<Option<u16>> = (0..n)
        .map(|i| mol.atom(chematic_core::AtomIdx(i as u32)).isotope)
        .collect();
    let has_isotopes = isotopes.iter().any(|iso| iso.is_some());

    // Build atoms object
    let mut atoms_obj = json!({
        "elements": { "number": numbers }
    });

    if !coords.is_empty() && coords.len() == n {
        let flat: Vec<f64> = coords.iter().flat_map(|&(x, y, z)| [x, y, z]).collect();
        atoms_obj["coords"] = json!({ "3d": flat });
    }

    if has_charges {
        atoms_obj["formalCharges"] = json!(charges.iter().map(|&c| c as i64).collect::<Vec<_>>());
    }

    if has_isotopes {
        let iso_vals: Vec<u64> = isotopes
            .iter()
            .map(|iso| iso.map(|n| n as u64).unwrap_or(0))
            .collect();
        atoms_obj["isotopes"] = json!(iso_vals);
    }

    // Bonds
    let mut bond_index: Vec<u64> = Vec::new();
    let mut bond_orders: Vec<Value> = Vec::new();
    for (_, bond) in mol.bonds() {
        bond_index.push(bond.atom1.0 as u64);
        bond_index.push(bond.atom2.0 as u64);
        bond_orders.push(json!(bond_order_to_float(bond.order)));
    }

    let bonds_obj = json!({
        "connections": { "index": bond_index },
        "order": bond_orders
    });

    let root = json!({
        "chemicalJson": 1,
        "atoms": atoms_obj,
        "bonds": bonds_obj,
    });

    serde_json::to_string(&root).unwrap_or_default()
}

/// Map a [`BondOrder`] to the float value used in ChemicalJSON.
fn bond_order_to_float(order: BondOrder) -> f64 {
    match order {
        BondOrder::Single
        | BondOrder::Up
        | BondOrder::Down
        | BondOrder::Zero
        | BondOrder::Dative => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Aromatic => 1.5,
        // Query bonds: fall back to single
        BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => 1.0,
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::BondOrder;

    fn make_ethanol_cjson() -> &'static str {
        r#"{
  "chemicalJson": 1,
  "atoms": {
    "coords": { "3d": [0.0,0.0,0.0, 1.54,0.0,0.0, 2.5,1.0,0.0] },
    "elements": { "number": [6, 6, 8] },
    "formalCharges": [0, 0, 0]
  },
  "bonds": {
    "connections": { "index": [0,1, 1,2] },
    "order": [1, 1]
  }
}"#
    }

    fn make_benzene_cjson() -> &'static str {
        r#"{
  "chemicalJson": 1,
  "atoms": {
    "elements": { "number": [6,6,6,6,6,6] }
  },
  "bonds": {
    "connections": { "index": [0,1, 1,2, 2,3, 3,4, 4,5, 5,0] },
    "order": [1.5, 1.5, 1.5, 1.5, 1.5, 1.5]
  }
}"#
    }

    #[test]
    fn parse_ethanol() {
        let (mol, coords) = parse_cjson(make_ethanol_cjson()).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        assert_eq!(coords.len(), 3);
        assert!((coords[0].0 - 0.0).abs() < 1e-9);
        assert!((coords[1].0 - 1.54).abs() < 1e-9);
        // O is atomic number 8
        let o_idx = chematic_core::AtomIdx(2);
        assert_eq!(mol.atom(o_idx).element.atomic_number(), 8);
    }

    #[test]
    fn parse_benzene_aromatic_bonds() {
        let (mol, coords) = parse_cjson(make_benzene_cjson()).unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        // coords absent → empty
        assert!(coords.is_empty());
        // All bonds should be aromatic
        for (_, bond) in mol.bonds() {
            assert_eq!(
                bond.order,
                BondOrder::Aromatic,
                "expected Aromatic, got {:?}",
                bond.order
            );
        }
    }

    #[test]
    fn parse_with_isotope() {
        let cjson = r#"{
  "chemicalJson": 1,
  "atoms": {
    "elements": { "number": [6, 1] },
    "isotopes": [13, 2]
  },
  "bonds": {
    "connections": { "index": [0,1] },
    "order": [1]
  }
}"#;
        let (mol, _) = parse_cjson(cjson).unwrap();
        assert_eq!(mol.atom(chematic_core::AtomIdx(0)).isotope, Some(13));
        assert_eq!(mol.atom(chematic_core::AtomIdx(1)).isotope, Some(2));
    }

    #[test]
    fn parse_with_charge() {
        let cjson = r#"{
  "chemicalJson": 1,
  "atoms": {
    "elements": { "number": [7] },
    "formalCharges": [1]
  },
  "bonds": { "connections": { "index": [] }, "order": [] }
}"#;
        let (mol, _) = parse_cjson(cjson).unwrap();
        assert_eq!(mol.atom(chematic_core::AtomIdx(0)).charge, 1);
    }

    #[test]
    fn roundtrip_with_coords() {
        let coords = vec![(0.0, 0.0, 0.0), (1.54, 0.0, 0.0), (2.5, 1.0, 0.0)];
        let (orig, _) = parse_cjson(make_ethanol_cjson()).unwrap();
        let written = write_cjson(&orig, &coords);
        let (parsed, coords2) = parse_cjson(&written).unwrap();
        assert_eq!(orig.atom_count(), parsed.atom_count());
        assert_eq!(orig.bond_count(), parsed.bond_count());
        assert_eq!(coords2.len(), 3);
        assert!((coords2[1].0 - 1.54).abs() < 1e-6);
    }

    #[test]
    fn write_no_coords() {
        let (mol, _) = parse_cjson(make_benzene_cjson()).unwrap();
        let written = write_cjson(&mol, &[]);
        let v: serde_json::Value = serde_json::from_str(&written).unwrap();
        // No coords.3d field when no coords given
        assert!(v["atoms"].get("coords").is_none());
        assert_eq!(v["chemicalJson"].as_u64(), Some(1));
    }

    #[test]
    fn error_unknown_atomic_number() {
        let cjson = r#"{"chemicalJson":1,"atoms":{"elements":{"number":[999]}},"bonds":{"connections":{"index":[]},"order":[]}}"#;
        assert!(matches!(
            parse_cjson(cjson),
            Err(CjsonError::UnknownAtomicNumber(999))
        ));
    }

    #[test]
    fn error_missing_atoms() {
        let cjson = r#"{"chemicalJson":1}"#;
        assert!(matches!(
            parse_cjson(cjson),
            Err(CjsonError::MissingField("atoms"))
        ));
    }

    #[test]
    fn bond_order_roundtrip() {
        for &(order, expected_float) in &[
            (BondOrder::Single, 1.0_f64),
            (BondOrder::Double, 2.0),
            (BondOrder::Triple, 3.0),
            (BondOrder::Aromatic, 1.5),
        ] {
            let f = bond_order_to_float(order);
            assert!((f - expected_float).abs() < 1e-9);
            assert_eq!(float_to_bond_order(f), order);
        }
    }

    #[test]
    fn bounded_parser_rejects_input_and_record_limits() {
        let input = make_benzene_cjson();
        assert!(matches!(
            parse_cjson_with_limits(
                input,
                &CjsonParseLimits {
                    max_input_bytes: 8,
                    ..Default::default()
                }
            ),
            Err(CjsonError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));
        assert!(matches!(
            parse_cjson_with_limits(
                input,
                &CjsonParseLimits {
                    max_atoms: 2,
                    ..Default::default()
                }
            ),
            Err(CjsonError::ResourceLimit {
                resource: "atom records",
                ..
            })
        ));
    }

    #[test]
    fn parser_rejects_numeric_truncation_and_wrong_types() {
        let oversized_atomic_number =
            r#"{"atoms":{"elements":{"number":[257]}},"bonds":{"connections":{"index":[]}}}"#;
        assert!(matches!(
            parse_cjson(oversized_atomic_number),
            Err(CjsonError::UnknownAtomicNumber(257))
        ));
        let oversized_isotope = r#"{"atoms":{"elements":{"number":[6]},"isotopes":[65536]},"bonds":{"connections":{"index":[]}}}"#;
        assert!(matches!(
            parse_cjson(oversized_isotope),
            Err(CjsonError::InvalidValue {
                field: "atoms.isotopes"
            })
        ));
        let invalid_charge = r#"{"atoms":{"elements":{"number":[6]},"formalCharges":[128]},"bonds":{"connections":{"index":[]}}}"#;
        assert!(matches!(
            parse_cjson(invalid_charge),
            Err(CjsonError::InvalidValue {
                field: "atoms.formalCharges"
            })
        ));
    }
}
