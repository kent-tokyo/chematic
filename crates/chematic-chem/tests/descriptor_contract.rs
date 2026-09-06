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

#[test]
fn shared_fingerprint_fixture_freezes_shape_and_configuration() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["fingerprint_contract"];
    assert_eq!(contract["schema_version"], 1);
    assert_eq!(contract["operations"]["ecfp4"]["bits"], 2048);
    assert_eq!(contract["operations"]["ecfp4"]["bytes"], 256);
    assert_eq!(
        contract["operations"]["ecfp4"]["configuration"]["radius"],
        2
    );
    assert_eq!(contract["operations"]["maccs"]["bits"], 166);
    assert_eq!(contract["operations"]["maccs"]["bytes"], 21);

    for fixture in contract["fixtures"]
        .as_array()
        .expect("fingerprint fixtures")
    {
        let id = fixture["id"].as_str().unwrap();
        let smiles = fixture["smiles"].as_str().unwrap();
        let mol = chematic_smiles::parse(smiles)
            .unwrap_or_else(|error| panic!("fingerprint fixture {id} must parse: {error}"));
        let ecfp4 = chematic_fp::ecfp4(&mol);
        let maccs = chematic_fp::maccs(&mol);
        assert_eq!(ecfp4.to_bitvecn().bit_width(), 2048, "{id} ECFP4 shape");
        assert_eq!(
            maccs.to_bitvecn().bit_width(),
            2048,
            "{id} MACCS backing shape"
        );
        assert!(ecfp4.popcount() > 0, "{id} ECFP4 must not be empty");
        assert!(maccs.popcount() > 0, "{id} MACCS must not be empty");
        assert!(
            (166..2048).all(|bit| !maccs.get(bit)),
            "{id} MACCS upper bits"
        );
    }
}
