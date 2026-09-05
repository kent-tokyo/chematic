//! Shared descriptor contract consumed by the Rust, Python, and WASM suites.

use serde_json::Value;

use chematic_chem::{hba_count, hbd_count, molecular_weight, tpsa};

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cross_binding_contract.json"
));

#[test]
fn shared_descriptor_fixture_matches_rust_source_of_truth() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["descriptor_contract"];
    assert_eq!(contract["schema_version"], 1);
    assert_eq!(contract["fields"]["molecular_weight"]["unit"], "Da");
    assert_eq!(contract["fields"]["tpsa"]["unit"], "A2");

    for fixture in contract["fixtures"]
        .as_array()
        .expect("descriptor fixtures")
    {
        let id = fixture["id"].as_str().unwrap();
        let smiles = fixture["smiles"].as_str().unwrap();
        let mol = chematic_smiles::parse(smiles)
            .unwrap_or_else(|error| panic!("descriptor fixture {id} must parse: {error}"));
        let close = |actual: f64, key: &str| {
            let expected = fixture[key].as_f64().unwrap();
            assert!(
                (actual - expected).abs() <= 1e-6,
                "{id} {key}: actual={actual} expected={expected}"
            );
        };
        close(molecular_weight(&mol), "molecular_weight");
        close(tpsa(&mol), "tpsa");
        assert_eq!(
            hbd_count(&mol),
            fixture["hbd"].as_u64().unwrap() as usize,
            "{id} hbd"
        );
        assert_eq!(
            hba_count(&mol),
            fixture["hba"].as_u64().unwrap() as usize,
            "{id} hba"
        );
        assert_eq!(
            chematic_chem::heavy_atom_count(&mol),
            fixture["heavy_atoms"].as_u64().unwrap() as usize,
            "{id} heavy_atoms"
        );
    }
}
