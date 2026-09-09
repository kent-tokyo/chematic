use serde_json::Value;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cross_binding_contract.json"
));

#[test]
fn shared_inchi_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["inchi_contract"];
    assert_eq!(contract["schema_version"], 1);
    for case in contract["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles).expect("fixture SMILES must parse");
        let inchi = chematic_inchi::inchi(&molecule);
        assert_eq!(inchi, case["inchi"], "InChI mismatch for {smiles}");
        assert_eq!(
            chematic_inchi::inchi_key(&inchi),
            case["inchikey"],
            "InChIKey mismatch for {smiles}"
        );
    }
}
