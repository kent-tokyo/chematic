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
    let operations = document["operation_manifest"]["operations"]
        .as_array()
        .expect("operation manifest operations");
    assert_eq!(operations.len(), 56);
    let mut operation_ids = std::collections::HashSet::new();
    for operation in operations {
        let id = operation["id"].as_str().expect("operation id");
        assert!(operation_ids.insert(id), "duplicate operation id: {id}");
        assert_eq!(operation["bindings"].as_array().unwrap().len(), 4);
        assert!(!operation["test_anchors"].as_array().unwrap().is_empty());
    }

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

#[test]
fn shared_smiles_validity_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["smiles_validity_contract"];
    assert_eq!(contract["schema_version"], 1);
    for value in contract["accepted"].as_array().unwrap() {
        let smiles = value.as_str().unwrap();
        assert!(chematic_smiles::parse(smiles).is_ok(), "accepted: {smiles}");
    }
    for value in contract["rejected"].as_array().unwrap() {
        let smiles = value.as_str().unwrap();
        assert!(
            chematic_smiles::parse(smiles).is_err(),
            "rejected: {smiles}"
        );
    }
}

#[test]
fn shared_batch_canonicalization_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["batch_canonicalization_contract"];
    assert_eq!(contract["schema_version"], 1);
    let inputs: Vec<&str> = contract["inputs"]
        .as_array()
        .expect("batch inputs")
        .iter()
        .map(|value| value.as_str().unwrap())
        .collect();
    let actual = chematic_smiles::SmilesBatchCanonicalizer::default().canonicalize(inputs);
    for (record, expected) in actual.iter().zip(contract["expected"].as_array().unwrap()) {
        assert_eq!(
            record.input_index,
            expected["input_index"].as_u64().unwrap() as usize
        );
        match (&record.result, expected["status"].as_str().unwrap()) {
            (chematic_smiles::BatchCanonicalization::Accepted { canonical_smiles }, "accepted") => {
                assert_eq!(
                    canonical_smiles,
                    expected["canonical_smiles"].as_str().unwrap()
                )
            }
            (chematic_smiles::BatchCanonicalization::Rejected { .. }, "rejected") => {}
            (result, status) => panic!("unexpected batch result {result:?} for {status}"),
        }
    }
}

#[test]
fn shared_rdkit_rdk_fingerprint_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let operation = &document["fingerprint_contract"]["operations"]["rdkit_rdk"];
    assert_eq!(operation["bits"], 2048);
    assert_eq!(operation["bytes"], 256);
    for fixture in document["fingerprint_contract"]["rdkit_rdk_fixtures"]
        .as_array()
        .expect("RDK fixtures")
    {
        let molecule = chematic_smiles::parse(fixture["smiles"].as_str().unwrap()).unwrap();
        let fp = chematic_fp::rdkit_rdk_fp(&molecule);
        let expected: Vec<usize> = fixture["rdkit_rdk_bits"]
            .as_array()
            .unwrap()
            .iter()
            .map(|bit| bit.as_u64().unwrap() as usize)
            .collect();
        let actual: Vec<usize> = (0..2048).filter(|&bit| fp.get(bit)).collect();
        assert_eq!(
            actual, expected,
            "RDK fingerprint mismatch for {}",
            fixture["id"]
        );
    }
}

#[test]
fn shared_rdkit_path_fingerprint_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    assert_eq!(
        document["fingerprint_contract"]["operations"]["rdkit_path"]["bytes"],
        256
    );
    for fixture in document["fingerprint_contract"]["rdkit_path_fixtures"]
        .as_array()
        .expect("path fixtures")
    {
        let molecule = chematic_smiles::parse(fixture["smiles"].as_str().unwrap()).unwrap();
        let fp = chematic_fp::rdkit_path_fp(&molecule);
        let expected: Vec<usize> = fixture["rdkit_path_bits"]
            .as_array()
            .unwrap()
            .iter()
            .map(|bit| bit.as_u64().unwrap() as usize)
            .collect();
        let actual: Vec<usize> = (0..2048).filter(|&bit| fp.get(bit)).collect();
        assert_eq!(
            actual, expected,
            "path fingerprint mismatch for {}",
            fixture["id"]
        );
    }
}
