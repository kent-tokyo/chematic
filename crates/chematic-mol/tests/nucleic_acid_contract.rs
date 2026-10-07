use chematic_mol::{
    NucleicAcidDocument, NucleicAcidLimits, apply_nucleic_acid_json_command,
    validate_nucleic_acid_json,
};
use serde_json::Value;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/nucleic_acid_document_contract.json"
));

#[test]
fn rust_contract_roundtrips_and_edits_valid_documents() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    for case in fixture["valid_cases"].as_array().unwrap() {
        let document = NucleicAcidDocument::from_json(&case["document"]).unwrap();
        assert_eq!(document.to_json(), case["document"], "{}", case["id"]);

        let envelope = apply_nucleic_acid_json_command(
            &case["document"].to_string(),
            &case["command"].to_string(),
            &NucleicAcidLimits::default(),
        );
        assert_eq!(envelope["ok"], true, "{}: {envelope}", case["id"]);
        let mut value = &envelope["document"];
        for segment in case["edited_path"].as_array().unwrap() {
            value = if let Some(key) = segment.as_str() {
                &value[key]
            } else {
                &value[segment.as_u64().unwrap() as usize]
            };
        }
        assert_eq!(*value, case["edited_value"], "{}", case["id"]);
    }
}

#[test]
fn rust_contract_reports_the_expected_typed_rejections() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    for case in fixture["invalid_cases"].as_array().unwrap() {
        let envelope = validate_nucleic_acid_json(
            &case["document"].to_string(),
            &NucleicAcidLimits::default(),
        );
        assert_eq!(envelope["ok"], false, "{}: {envelope}", case["id"]);
        assert_eq!(
            envelope["error"]["code"], case["expected_code"],
            "{}: {envelope}",
            case["id"]
        );
    }
}

#[test]
fn collection_limit_is_enforced_before_acceptance() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let limits = NucleicAcidLimits {
        max_strands: 0,
        ..NucleicAcidLimits::default()
    };
    let envelope =
        validate_nucleic_acid_json(&fixture["valid_cases"][0]["document"].to_string(), &limits);
    assert_eq!(envelope["ok"], false);
    assert_eq!(envelope["error"]["code"], "resource_limit");
}
