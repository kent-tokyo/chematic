//! Shared binding contract fixture.
//!
//! The same JSON file is consumed by this Rust test, the Python binding test,
//! and the Node-hosted WASM test. Keeping the expectations in one checked-in
//! fixture prevents the binding suites from silently drifting apart.

use serde_json::Value;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cross_binding_contract.json"
));

#[test]
fn shared_parse_and_canonical_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    assert_eq!(document["schema_version"], 1);
    let fixtures = document["fixtures"].as_array().expect("fixtures array");
    assert_eq!(fixtures.len(), 4);

    for fixture in fixtures {
        let id = fixture["id"].as_str().unwrap();
        let smiles = fixture["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles)
            .unwrap_or_else(|error| panic!("fixture {id} ({smiles}) must parse: {error}"));
        assert_eq!(
            chematic_smiles::canonical_smiles(&molecule),
            fixture["canonical_smiles"].as_str().unwrap(),
            "canonical output mismatch for {id}"
        );
        assert_eq!(
            molecule.atom_count(),
            fixture["heavy_atoms"].as_u64().unwrap() as usize,
            "atom count mismatch for {id}"
        );
    }
}
