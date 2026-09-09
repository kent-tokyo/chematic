use serde_json::Value;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cross_binding_contract.json"
));

#[test]
fn shared_smarts_validity_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["smarts_validity_contract"];
    assert_eq!(contract["schema_version"], 1);
    for value in contract["accepted"].as_array().unwrap() {
        let smarts = value.as_str().unwrap();
        assert!(
            chematic_smarts::parse_smarts(smarts).is_ok(),
            "accepted: {smarts}"
        );
    }
    for value in contract["rejected"].as_array().unwrap() {
        let smarts = value.as_str().unwrap();
        assert!(
            chematic_smarts::parse_smarts(smarts).is_err(),
            "rejected: {smarts}"
        );
    }
}
