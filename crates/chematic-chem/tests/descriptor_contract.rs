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
fn shared_standardization_profile_matches_rust_source_of_truth() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["standardization_contract"];
    assert_eq!(contract["schema_version"], 1);
    assert_eq!(contract["profile"]["largest_fragment_only"], true);
    for fixture in contract["fixtures"]
        .as_array()
        .expect("standardization fixtures")
    {
        let id = fixture["id"].as_str().unwrap();
        let smiles = fixture["smiles"].as_str().unwrap();
        let mol = chematic_smiles::parse(smiles)
            .unwrap_or_else(|error| panic!("standardization fixture {id} must parse: {error}"));
        let options = chematic_chem::StandardizeOptions {
            largest_fragment_only: true,
            ..Default::default()
        };
        let output = chematic_chem::standardize(&mol, &options);
        assert_eq!(
            chematic_smiles::canonical_smiles(&output),
            fixture["output_smiles"].as_str().unwrap(),
            "{id}"
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

#[test]
fn shared_fingerprint_detail_contract_matches_rust_source_of_truth() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["fingerprint_detail_contract"];
    assert_eq!(contract["schema_version"], 1);
    let operation = &contract["operations"]["rdkit_ecfp4_detail"];
    assert_eq!(operation["bits"], 2048);
    assert_eq!(operation["bytes"], 256);
    assert_eq!(operation["configuration"]["radius"], 2);
    assert_eq!(
        operation["configuration"]["include_redundant_environments"],
        false
    );
    assert_eq!(
        operation["explanation"]["radius_range"],
        serde_json::json!([0, 2])
    );

    for fixture in contract["fixtures"].as_array().expect("detail fixtures") {
        let id = fixture["id"].as_str().unwrap();
        let smiles = fixture["smiles"].as_str().unwrap();
        let mol = chematic_smiles::parse(smiles)
            .unwrap_or_else(|error| panic!("detail fixture {id} must parse: {error}"));
        let detail = chematic_fp::rdkit_morgan_ecfp4_experimental(&mol)
            .unwrap_or_else(|error| panic!("detail fixture {id} failed: {error}"));
        assert_eq!(
            detail.fingerprint.to_bitvecn().bit_width(),
            2048,
            "{id} bit shape"
        );
        assert!(!detail.sparse_counts.is_empty(), "{id} sparse counts");
        assert!(!detail.raw_bit_info.is_empty(), "{id} raw explanation");
        assert!(
            !detail.folded_bit_info.is_empty(),
            "{id} folded explanation"
        );
        for pairs in detail.raw_bit_info.values() {
            for &(atom, radius) in pairs {
                assert!((radius as u64) <= 2, "{id} radius");
                assert!((atom as usize) < mol.atom_count(), "{id} atom provenance");
            }
        }
        for pairs in detail.folded_bit_info.values() {
            for &(atom, radius) in pairs {
                assert!((radius as u64) <= 2, "{id} radius");
                assert!((atom as usize) < mol.atom_count(), "{id} atom provenance");
            }
        }
    }
}
