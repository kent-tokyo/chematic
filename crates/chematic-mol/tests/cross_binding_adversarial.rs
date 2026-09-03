//! Shared adversarial parser contract.
//!
//! This is intentionally driven by the same fixture consumed by the Python
//! and Node-hosted WASM tests. Every common topology parser must reject the
//! malformed input instead of panicking, accepting a partial molecule, or
//! silently falling back to an empty result.

use serde_json::Value;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cross_binding_contract.json"
));

#[test]
fn every_common_topology_parser_rejects_shared_adversarial_inputs() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let cases = document["adversarial"]
        .as_array()
        .expect("adversarial array");
    assert_eq!(cases.len(), 8);

    for case in cases {
        let id = case["id"].as_str().unwrap();
        let format = case["format"].as_str().unwrap();
        let input = case["input"].as_str().unwrap();
        let accepted = match format {
            "smiles" => chematic_smiles::parse(input).is_ok(),
            "mol" => chematic_mol::parse_mol(input).is_ok(),
            "mol_v3000" => chematic_mol::parse_mol_v3000(input).is_ok(),
            "mol2" => chematic_mol::parse_mol2(input).is_ok(),
            "cml" => chematic_mol::parse_cml(input).is_ok(),
            "cjson" => chematic_mol::parse_cjson(input).is_ok(),
            "moljson" => chematic_mol::parse_moljson(input).is_ok(),
            "cdxml" => chematic_mol::parse_cdxml(input).is_ok(),
            other => panic!("unregistered adversarial parser format: {other}"),
        };
        assert!(!accepted, "adversarial fixture was accepted: {id}");
    }
}
