use chematic_chem::rdkit_num_rings;
use chematic_smiles::parse;
use serde_json::Value;

#[test]
fn pinned_rdkit_ring_oracle_versions_are_not_conflated() {
    let manifest: Value = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/rdkit_ring_count_profiles_v1.json"
    )))
    .unwrap();
    assert_eq!(manifest["schema_version"], 1);
    assert_eq!(manifest["operation"], "CalcNumRings");

    let fixtures = manifest["fixtures"].as_array().unwrap();
    assert_eq!(fixtures.len(), 4);
    for fixture in fixtures {
        assert_eq!(fixture["classification"], "oracle_version_difference");
        let expected = &fixture["expected"];
        assert_ne!(
            expected["rdkit_2026_03_6"], expected["rdkit_2026_09_1"],
            "fixture {} no longer distinguishes the pinned oracles",
            fixture["id"]
        );
        let molecule = parse(fixture["smiles"].as_str().unwrap()).unwrap();
        assert_eq!(
            rdkit_num_rings(&molecule) as u64,
            expected["chematic_default"].as_u64().unwrap(),
            "{}",
            fixture["id"]
        );
    }
}
