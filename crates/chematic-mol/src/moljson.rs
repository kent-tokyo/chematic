//! MolJSON format parser and writer.
//!
//! MolJSON is a JSON-based molecular representation designed for LLM
//! (large language model) compatibility.  Unlike SMILES, it makes atoms,
//! bonds, and connectivity explicit without requiring domain-specific parsing
//! rules, which makes it easier for LLMs to generate and interpret.
//!
//! The atom-record dialect below is chematic's existing representation; it
//! differs from the upstream oxpig schema (which uses `charges` and
//! `aromatic_n_h`). This API does not implement that upstream schema.
//!
//! ## Format overview
//!
//! ```json
//! {
//!   "atoms": [
//!     {"id": "a1", "element": "C", "charge": 0, "isotope": null,
//!      "hydrogens": 3, "aromatic": false},
//!     {"id": "a2", "element": "O", "charge": 0, "isotope": null,
//!      "hydrogens": 1, "aromatic": false}
//!   ],
//!   "bonds": [
//!     {"id": "b1", "source_id": "a1", "target_id": "a2",
//!      "order": 1.0, "aromatic": false}
//!   ]
//! }
//! ```
//!
//! ## References
//! - Nicholas T. Runcie, Fergus Imrie, Charlotte M. Deane, "Molecular
//!   Representations for Large Language Models," arXiv:2605.01822, 2026.
//! - Official schema and reference implementation: oxpig/MolJSON
//!   <https://github.com/oxpig/MolJSON> (reference implementation licensed
//!   BSD-3-Clause).

#![forbid(unsafe_code)]

use chematic_core::{
    Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder, implicit_hcount,
};
use std::collections::{HashMap, HashSet};

mod fields;

/// Resource limits for MolJSON parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolJsonParseLimits {
    pub max_input_bytes: usize,
    pub max_json_depth: usize,
    pub max_array_items: usize,
    pub max_string_bytes: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
}

impl Default for MolJsonParseLimits {
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
pub enum MolJsonError {
    InvalidJson(String),
    UnknownElement(String),
    InvalidBondRef {
        bond_id: String,
        ref_id: String,
    },
    MissingField(&'static str),
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl std::fmt::Display for MolJsonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MolJsonError::InvalidJson(s) => write!(f, "MolJSON: invalid JSON: {s}"),
            MolJsonError::UnknownElement(s) => write!(f, "MolJSON: unknown element symbol '{s}'"),
            MolJsonError::InvalidBondRef { bond_id, ref_id } => write!(
                f,
                "MolJSON: bond '{bond_id}' references unknown atom id '{ref_id}'"
            ),
            MolJsonError::MissingField(fld) => {
                write!(f, "MolJSON: missing required field '{fld}'")
            }
            MolJsonError::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "MolJSON: {resource} has size {actual}, exceeding limit {limit}"
            ),
        }
    }
}

impl std::error::Error for MolJsonError {}

// ─── Parser ──────────────────────────────────────────────────────────────────

/// Parse a MolJSON string into a [`Molecule`].
///
/// Atom IDs (`"a1"`, `"a2"`, …) are used only internally during parsing to
/// resolve bond references; they are not stored on the returned molecule.
/// A present `hydrogens` field is an integer in 0..=255 giving the atom's
/// non-node H count, including an explicit zero. H atoms present as graph
/// nodes are additional and are not counted here. When the field is absent,
/// the existing valence-based inference remains in effect. This is not a
/// valence or radical-state validator.
///
/// Known charge/isotope/boolean fields are checked before conversion. Numeric
/// bond orders are exactly 1, 1.5, 2, 3 or 4, with consistent aromatic flags;
/// 4 denotes a quadruple bond. Duplicate JSON keys, IDs, self/duplicate bonds
/// and malformed present fields are errors. Unknown fields remain outside the
/// graph projection; atom/bond IDs are used for references, not source retention.
pub fn parse_moljson(input: &str) -> Result<Molecule, MolJsonError> {
    parse_moljson_with_limits(input, &MolJsonParseLimits::default())
}

fn hydrogen_count(value: &serde_json::Value) -> Result<Option<u8>, MolJsonError> {
    value
        .get("hydrogens")
        .map(|count| {
            count
                .as_u64()
                .and_then(|count| u8::try_from(count).ok())
                .ok_or_else(|| {
                    MolJsonError::InvalidJson(
                        "atoms[].hydrogens must be an integer in 0..=255".into(),
                    )
                })
        })
        .transpose()
}

fn parse_moljson_atoms(
    atoms: &[serde_json::Value],
) -> Result<(MoleculeBuilder, HashMap<String, AtomIdx>), MolJsonError> {
    let mut builder = MoleculeBuilder::new();
    let mut ids = HashMap::new();
    for value in atoms {
        let id = fields::id(value, "id", "atoms[].id")?.to_owned();
        if ids.contains_key(&id) {
            return Err(fields::invalid(format!("duplicate atom id {id:?}")));
        }
        let symbol = value
            .get("element")
            .and_then(|value| value.as_str())
            .ok_or(MolJsonError::MissingField("atoms[].element"))?;
        let element = Element::from_symbol(symbol)
            .ok_or_else(|| MolJsonError::UnknownElement(symbol.to_owned()))?;
        let mut atom = Atom::new(element);
        atom.hydrogen_count = hydrogen_count(value)?;
        atom.charge = fields::charge(value)?;
        atom.isotope = fields::isotope(value)?;
        atom.aromatic = fields::flag(value, "atoms[].aromatic")?.unwrap_or(false);
        let index = builder.add_atom(atom);
        ids.insert(id, index);
    }
    Ok((builder, ids))
}

fn parse_moljson_bonds(
    bonds: &[serde_json::Value],
    ids: &HashMap<String, AtomIdx>,
    builder: &mut MoleculeBuilder,
) -> Result<(), MolJsonError> {
    let mut bond_ids = HashSet::new();
    for value in bonds {
        if !value.is_object() {
            return Err(fields::invalid("bonds[] must be an object"));
        }
        let bond_id = if value.get("id").is_some() {
            let id = fields::id(value, "id", "bonds[].id")?;
            if !bond_ids.insert(id) {
                return Err(fields::invalid(format!("duplicate bond id {id:?}")));
            }
            id
        } else {
            "?"
        }
        .to_owned();
        let endpoint = |key: &str, field: &'static str| {
            let id = value
                .get(key)
                .and_then(|value| value.as_str())
                .ok_or(MolJsonError::MissingField(field))?;
            ids.get(id)
                .copied()
                .ok_or_else(|| MolJsonError::InvalidBondRef {
                    bond_id: bond_id.clone(),
                    ref_id: id.to_owned(),
                })
        };
        let source = endpoint("source_id", "bonds[].source_id")?;
        let target = endpoint("target_id", "bonds[].target_id")?;
        if source == target {
            return Err(fields::invalid(format!("bond {bond_id:?} is a self bond")));
        }
        let order = fields::bond_order(value)?;
        builder
            .add_bond(source, target, order)
            .map_err(|error| fields::invalid(format!("bond {bond_id:?}: {error}")))?;
    }
    Ok(())
}

/// Parse MolJSON with explicit resource limits.
pub fn parse_moljson_with_limits(
    input: &str,
    limits: &MolJsonParseLimits,
) -> Result<Molecule, MolJsonError> {
    if input.len() > limits.max_input_bytes {
        return Err(MolJsonError::ResourceLimit {
            resource: "input bytes",
            actual: input.len(),
            limit: limits.max_input_bytes,
        });
    }
    let v: crate::cjson_record::UniqueValue =
        serde_json::from_str(input).map_err(|e| MolJsonError::InvalidJson(e.to_string()))?;
    let v = v.0;
    check_json_limits(&v, "$", 0, limits)?;

    let atoms_arr = v
        .get("atoms")
        .and_then(|a| a.as_array())
        .ok_or(MolJsonError::MissingField("atoms"))?;
    if atoms_arr.len() > limits.max_atoms {
        return Err(MolJsonError::ResourceLimit {
            resource: "atom records",
            actual: atoms_arr.len(),
            limit: limits.max_atoms,
        });
    }

    let (mut builder, id_to_idx) = parse_moljson_atoms(atoms_arr)?;

    // ── Bonds ────────────────────────────────────────────────────────────────
    if let Some(bonds) = v.get("bonds") {
        let bonds_arr = bonds
            .as_array()
            .ok_or_else(|| fields::invalid("bonds must be an array when present"))?;
        if bonds_arr.len() > limits.max_bonds {
            return Err(MolJsonError::ResourceLimit {
                resource: "bond records",
                actual: bonds_arr.len(),
                limit: limits.max_bonds,
            });
        }
        parse_moljson_bonds(bonds_arr, &id_to_idx, &mut builder)?;
    }

    Ok(builder.build())
}

fn check_json_limits(
    value: &serde_json::Value,
    path: &str,
    depth: usize,
    limits: &MolJsonParseLimits,
) -> Result<(), MolJsonError> {
    if depth > limits.max_json_depth {
        return Err(MolJsonError::ResourceLimit {
            resource: "JSON depth",
            actual: depth,
            limit: limits.max_json_depth,
        });
    }
    match value {
        serde_json::Value::String(s) if s.len() > limits.max_string_bytes => {
            return Err(MolJsonError::ResourceLimit {
                resource: "string bytes",
                actual: s.len(),
                limit: limits.max_string_bytes,
            });
        }
        serde_json::Value::Array(values) => {
            if values.len() > limits.max_array_items {
                return Err(MolJsonError::ResourceLimit {
                    resource: "array items",
                    actual: values.len(),
                    limit: limits.max_array_items,
                });
            }
            for (index, child) in values.iter().enumerate() {
                check_json_limits(child, &format!("{path}[{index}]"), depth + 1, limits)?;
            }
        }
        serde_json::Value::Object(object) => {
            for (key, child) in object {
                if key.len() > limits.max_string_bytes {
                    return Err(MolJsonError::ResourceLimit {
                        resource: "object key bytes",
                        actual: key.len(),
                        limit: limits.max_string_bytes,
                    });
                }
                check_json_limits(child, &format!("{path}.{key}"), depth + 1, limits)?;
            }
        }
        _ => {}
    }
    Ok(())
}

// ─── Writer ──────────────────────────────────────────────────────────────────

/// Serialize a [`Molecule`] to a MolJSON string (pretty-printed).
///
/// Atom IDs are assigned as `"a1"`, `"a2"`, … in molecule atom order.
/// Bond IDs are assigned as `"b1"`, `"b2"`, … in bond iteration order.
/// The `hydrogens` field reflects computed implicit H count.
/// This graph-only writer retains Single/Double/Triple/Quadruple/Aromatic orders.
/// Native zero/dative/stereo/query orders still project to single bonds; atom
/// stereo, atom maps, coordinates and source metadata are not encoded.
pub fn write_moljson(mol: &Molecule) -> String {
    use serde_json::{Value, json};

    let n = mol.atom_count();

    let atoms: Vec<Value> = (0..n)
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let atom = mol.atom(idx);
            let h = implicit_hcount(mol, idx);
            let isotope: Value = match atom.isotope {
                Some(iso) => json!(iso as u64),
                None => Value::Null,
            };
            json!({
                "id": format!("a{}", i + 1),
                "element": atom.element.symbol(),
                "charge": atom.charge as i64,
                "isotope": isotope,
                "hydrogens": h as u64,
                "aromatic": atom.aromatic
            })
        })
        .collect();

    let bonds: Vec<Value> = mol
        .bonds()
        .enumerate()
        .map(|(i, (_, bond))| {
            let aromatic = bond.order == BondOrder::Aromatic;
            let order = bond_order_to_float(bond.order);
            json!({
                "id": format!("b{}", i + 1),
                "source_id": format!("a{}", bond.atom1.0 + 1),
                "target_id": format!("a{}", bond.atom2.0 + 1),
                "order": order,
                "aromatic": aromatic
            })
        })
        .collect();

    let root = json!({
        "atoms": atoms,
        "bonds": bonds
    });

    serde_json::to_string_pretty(&root).unwrap_or_default()
}

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
    use chematic_smiles::{canonical_smiles, parse};

    fn smiles_roundtrip(smi: &str) -> String {
        let mol = parse(smi).expect("parse SMILES");
        let json = write_moljson(&mol);
        let mol2 = parse_moljson(&json).expect("parse MolJSON");
        canonical_smiles(&mol2)
    }

    #[test]
    fn roundtrip_ethanol() {
        let cs = smiles_roundtrip("CCO");
        // Canonical SMILES for ethanol
        assert_eq!(cs, canonical_smiles(&parse("CCO").unwrap()));
    }

    #[test]
    fn shared_moljson_roundtrip_contract_matches() {
        let document: serde_json::Value = serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cross_binding_contract.json"
        )))
        .expect("shared contract JSON");
        let contract = &document["moljson_contract"];
        let mol = parse_moljson(contract["input"].as_str().unwrap()).unwrap();
        let roundtripped = parse_moljson(&write_moljson(&mol)).unwrap();
        let expected = &contract["expected"];
        assert_eq!(
            mol.atom_count(),
            expected["atom_count"].as_u64().unwrap() as usize
        );
        assert_eq!(
            canonical_smiles(&roundtripped),
            expected["canonical_smiles"].as_str().unwrap()
        );
    }

    #[test]
    fn roundtrip_benzene() {
        let cs = smiles_roundtrip("c1ccccc1");
        assert_eq!(cs, canonical_smiles(&parse("c1ccccc1").unwrap()));
    }

    #[test]
    fn roundtrip_charged() {
        // [Cl-]: unambiguous H count (0), roundtrips cleanly
        let cs = smiles_roundtrip("[Cl-]");
        assert_eq!(cs, canonical_smiles(&parse("[Cl-]").unwrap()));
    }

    #[test]
    fn roundtrip_charged_carbon() {
        // Simple negatively charged carbon
        let mol = parse("CC[CH2-]").unwrap();
        let json_str = write_moljson(&mol);
        let mol2 = parse_moljson(&json_str).unwrap();
        assert_eq!(mol2.atom_count(), 3);
        assert_eq!(mol2.atom(AtomIdx(2)).charge, -1);
    }

    #[test]
    fn roundtrip_isotope() {
        let mol = parse("[13C]").unwrap();
        let json = write_moljson(&mol);
        let mol2 = parse_moljson(&json).unwrap();
        assert_eq!(mol2.atom(AtomIdx(0)).isotope, Some(13));
    }

    #[test]
    fn roundtrip_aromatic_nh() {
        // The pyrrolic N-H count is semantic data, not a valence hint.
        let mol = parse("c1cc[nH]c1").unwrap();
        let json_str = write_moljson(&mol);
        let mol2 = parse_moljson(&json_str).unwrap();
        assert_eq!(mol2.atom_count(), 5);
        assert_eq!(canonical_smiles(&mol2), canonical_smiles(&mol));
        // All atoms should retain their aromatic flags
        for (idx, atom) in mol.atoms() {
            assert_eq!(
                atom.aromatic,
                mol2.atom(idx).aromatic,
                "atom {idx:?} aromatic mismatch"
            );
        }
    }

    #[test]
    fn hydrogens_count_in_output() {
        let mol = parse("CCO").unwrap();
        let json_str = write_moljson(&mol);
        let v: serde_json::Value = serde_json::from_str(&json_str).unwrap();
        // First atom (C in CH3-): 3 implicit H
        assert_eq!(v["atoms"][0]["hydrogens"].as_u64(), Some(3));
        // Third atom (O in -OH): 1 implicit H
        assert_eq!(v["atoms"][2]["hydrogens"].as_u64(), Some(1));
    }

    #[test]
    fn atom_ids_sequential() {
        let mol = parse("CC").unwrap();
        let json_str = write_moljson(&mol);
        let v: serde_json::Value = serde_json::from_str(&json_str).unwrap();
        assert_eq!(v["atoms"][0]["id"].as_str(), Some("a1"));
        assert_eq!(v["atoms"][1]["id"].as_str(), Some("a2"));
        assert_eq!(v["bonds"][0]["id"].as_str(), Some("b1"));
        assert_eq!(v["bonds"][0]["source_id"].as_str(), Some("a1"));
        assert_eq!(v["bonds"][0]["target_id"].as_str(), Some("a2"));
    }

    #[test]
    fn parse_rejects_invalid_json() {
        assert!(matches!(
            parse_moljson("not json"),
            Err(MolJsonError::InvalidJson(_))
        ));
    }

    #[test]
    fn parse_rejects_bad_bond_ref() {
        let bad = r#"{"atoms":[{"id":"a1","element":"C","charge":0,"isotope":null,"hydrogens":4,"aromatic":false}],"bonds":[{"id":"b1","source_id":"a1","target_id":"a99","order":1.0,"aromatic":false}]}"#;
        assert!(matches!(
            parse_moljson(bad),
            Err(MolJsonError::InvalidBondRef { .. })
        ));
    }

    #[test]
    fn parse_rejects_unknown_element() {
        let bad = r#"{"atoms":[{"id":"a1","element":"Xx","charge":0,"isotope":null,"hydrogens":0,"aromatic":false}],"bonds":[]}"#;
        assert!(matches!(
            parse_moljson(bad),
            Err(MolJsonError::UnknownElement(_))
        ));
    }

    #[test]
    fn roundtrip_aspirin() {
        let cs = smiles_roundtrip("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(
            cs,
            canonical_smiles(&parse("CC(=O)Oc1ccccc1C(=O)O").unwrap())
        );
    }

    #[test]
    fn bounded_parser_rejects_input_and_json_limits() {
        let input = r#"{"atoms":[],"bonds":[]}"#;
        let array_input = r#"{"atoms":[{}],"bonds":[]}"#;
        assert!(matches!(
            parse_moljson_with_limits(
                input,
                &MolJsonParseLimits {
                    max_input_bytes: 8,
                    ..Default::default()
                }
            ),
            Err(MolJsonError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));
        assert!(matches!(
            parse_moljson_with_limits(
                array_input,
                &MolJsonParseLimits {
                    max_array_items: 0,
                    ..Default::default()
                }
            ),
            Err(MolJsonError::ResourceLimit {
                resource: "array items",
                ..
            })
        ));
    }

    #[test]
    fn bounded_parser_rejects_atom_and_bond_limits() {
        let input = r#"{"atoms":[{"id":"a1","element":"C"},{"id":"a2","element":"O"}],"bonds":[{"id":"b1","source_id":"a1","target_id":"a2","order":1.0}]}"#;
        assert!(matches!(
            parse_moljson_with_limits(
                input,
                &MolJsonParseLimits {
                    max_atoms: 1,
                    ..Default::default()
                }
            ),
            Err(MolJsonError::ResourceLimit {
                resource: "atom records",
                ..
            })
        ));
        assert!(matches!(
            parse_moljson_with_limits(
                input,
                &MolJsonParseLimits {
                    max_bonds: 0,
                    ..Default::default()
                }
            ),
            Err(MolJsonError::ResourceLimit {
                resource: "bond records",
                ..
            })
        ));
    }
}
