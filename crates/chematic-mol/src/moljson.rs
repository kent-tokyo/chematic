//! MolJSON format parser and writer.
//!
//! MolJSON is a JSON-based molecular representation designed for LLM
//! (large language model) compatibility.  Unlike SMILES, it makes atoms,
//! bonds, and connectivity explicit without requiring domain-specific parsing
//! rules, which makes it easier for LLMs to generate and interpret.
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
use std::collections::HashMap;

// ─── Error ───────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MolJsonError {
    InvalidJson(String),
    UnknownElement(String),
    InvalidBondRef { bond_id: String, ref_id: String },
    MissingField(&'static str),
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
        }
    }
}

impl std::error::Error for MolJsonError {}

// ─── Parser ──────────────────────────────────────────────────────────────────

/// Parse a MolJSON string into a [`Molecule`].
///
/// Atom IDs (`"a1"`, `"a2"`, …) are used only internally during parsing to
/// resolve bond references; they are not stored on the returned molecule.
/// The `hydrogens` field is informational and is ignored during parsing —
/// implicit H counts are recomputed from valence rules by chematic-core.
pub fn parse_moljson(input: &str) -> Result<Molecule, MolJsonError> {
    let v: serde_json::Value =
        serde_json::from_str(input).map_err(|e| MolJsonError::InvalidJson(e.to_string()))?;

    let atoms_arr = v
        .get("atoms")
        .and_then(|a| a.as_array())
        .ok_or(MolJsonError::MissingField("atoms"))?;

    // ── Build atoms and id→AtomIdx map ──────────────────────────────────────
    let mut builder = MoleculeBuilder::new();
    let mut id_to_idx: HashMap<String, AtomIdx> = HashMap::new();

    for atom_val in atoms_arr {
        let id = atom_val
            .get("id")
            .and_then(|v| v.as_str())
            .ok_or(MolJsonError::MissingField("atoms[].id"))?
            .to_owned();

        let element_sym = atom_val
            .get("element")
            .and_then(|v| v.as_str())
            .ok_or(MolJsonError::MissingField("atoms[].element"))?;

        let element = Element::from_symbol(element_sym)
            .ok_or_else(|| MolJsonError::UnknownElement(element_sym.to_owned()))?;

        let charge = atom_val.get("charge").and_then(|v| v.as_i64()).unwrap_or(0) as i8;

        let isotope = atom_val.get("isotope").and_then(|v| {
            if v.is_null() {
                None
            } else {
                v.as_u64().map(|n| n as u16)
            }
        });

        let aromatic = atom_val
            .get("aromatic")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let mut atom = Atom::new(element);
        atom.charge = charge;
        atom.isotope = isotope;
        atom.aromatic = aromatic;

        let idx = builder.add_atom(atom);
        id_to_idx.insert(id, idx);
    }

    // ── Bonds ────────────────────────────────────────────────────────────────
    if let Some(bonds_arr) = v.get("bonds").and_then(|b| b.as_array()) {
        for bond_val in bonds_arr {
            let bond_id = bond_val
                .get("id")
                .and_then(|v| v.as_str())
                .unwrap_or("?")
                .to_owned();

            let src = bond_val
                .get("source_id")
                .and_then(|v| v.as_str())
                .ok_or(MolJsonError::MissingField("bonds[].source_id"))?;
            let tgt = bond_val
                .get("target_id")
                .and_then(|v| v.as_str())
                .ok_or(MolJsonError::MissingField("bonds[].target_id"))?;

            let src_idx =
                id_to_idx
                    .get(src)
                    .copied()
                    .ok_or_else(|| MolJsonError::InvalidBondRef {
                        bond_id: bond_id.clone(),
                        ref_id: src.to_owned(),
                    })?;
            let tgt_idx =
                id_to_idx
                    .get(tgt)
                    .copied()
                    .ok_or_else(|| MolJsonError::InvalidBondRef {
                        bond_id: bond_id.clone(),
                        ref_id: tgt.to_owned(),
                    })?;

            let aromatic = bond_val
                .get("aromatic")
                .and_then(|v| v.as_bool())
                .unwrap_or(false);
            let order_f = bond_val
                .get("order")
                .and_then(|v| v.as_f64())
                .unwrap_or(1.0);
            let order = if aromatic {
                BondOrder::Aromatic
            } else {
                float_to_bond_order(order_f)
            };

            let _ = builder.add_bond(src_idx, tgt_idx, order);
        }
    }

    Ok(builder.build())
}

fn float_to_bond_order(v: f64) -> BondOrder {
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

/// Serialize a [`Molecule`] to a MolJSON string (pretty-printed).
///
/// Atom IDs are assigned as `"a1"`, `"a2"`, … in molecule atom order.
/// Bond IDs are assigned as `"b1"`, `"b2"`, … in bond iteration order.
/// The `hydrogens` field reflects computed implicit H count.
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
        // Verify aromaticity flag and atom count survive roundtrip.
        // Note: bracket H notation ([nH]) is not stored explicitly and is
        // recomputed from valence rules; SMILES string may differ.
        let mol = parse("c1cc[nH]c1").unwrap();
        let json_str = write_moljson(&mol);
        let mol2 = parse_moljson(&json_str).unwrap();
        assert_eq!(mol2.atom_count(), 5);
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
}
