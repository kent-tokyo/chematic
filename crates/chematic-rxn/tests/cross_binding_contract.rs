use serde_json::Value;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cross_binding_contract.json"
));

#[test]
fn reaction_document_contract_round_trips_through_typed_model() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["rxn_document_contract"];
    assert_eq!(contract["schema_version"], 1);
    let authored = serde_json::to_string(&contract["document"]).unwrap();
    let parsed = chematic_rxn::ReactionDocument::from_json_str(&authored).unwrap();
    let observed: Vec<Value> = parsed.steps[0]
        .components
        .iter()
        .map(|component| serde_json::json!({"role": component.role, "smiles": component.smiles}))
        .collect();
    let expected: Vec<Value> =
        serde_json::from_value(contract["expected_components"].clone()).unwrap();
    assert_eq!(observed, expected);
}
