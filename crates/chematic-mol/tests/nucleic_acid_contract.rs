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

const SEQUENCES: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/nucleic_acid_edit_sequences.json"
));

/// Atom id -> owning residue id; panics when an atom is owned twice or a
/// declared atom is unowned.
fn ownership(document: &Value) -> std::collections::BTreeMap<String, String> {
    let mut owner = std::collections::BTreeMap::new();
    for strand in document["strands"].as_array().unwrap() {
        for residue in strand["residues"].as_array().unwrap() {
            for atom in residue["atom_refs"].as_array().unwrap() {
                let previous = owner.insert(
                    atom.as_str().unwrap().to_string(),
                    residue["id"].as_str().unwrap().to_string(),
                );
                assert!(previous.is_none(), "atom {atom} owned twice");
            }
        }
    }
    for atom in document["atom_ids"].as_array().unwrap() {
        assert!(owner.contains_key(atom.as_str().unwrap()), "atom {atom} unowned");
    }
    owner
}

#[test]
fn edit_sequences_reload_and_keep_atom_ownership() {
    let contract: Value = serde_json::from_str(FIXTURE).unwrap();
    let fixture: Value = serde_json::from_str(SEQUENCES).unwrap();
    let limits = NucleicAcidLimits::default();
    for sequence in fixture["sequences"].as_array().unwrap() {
        let id = sequence["id"].as_str().unwrap();
        let mut document = contract["valid_cases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|case| case["id"] == sequence["start"])
            .unwrap()["document"]
            .clone();
        let owners = ownership(&document);
        for (k, step) in sequence["steps"].as_array().unwrap().iter().enumerate() {
            let envelope = apply_nucleic_acid_json_command(
                &document.to_string(),
                &step["command"].to_string(),
                &limits,
            );
            assert_eq!(envelope["ok"], step["ok"], "{id} step {k}: {envelope}");
            if step["ok"] == true {
                assert_eq!(envelope["document"], step["document"], "{id} step {k}");
                document = envelope["document"].clone();
                // Serialized and read again, the edited document is the same.
                let reread = validate_nucleic_acid_json(&document.to_string(), &limits);
                assert_eq!(reread["ok"], true, "{id} step {k}");
                assert_eq!(reread["document"], document, "{id} step {k}");
                assert_eq!(ownership(&document), owners, "{id} step {k}");
            } else {
                assert_eq!(envelope["error"]["code"], step["code"], "{id} step {k}");
                assert_eq!(envelope["error"]["path"], step["path"], "{id} step {k}");
            }
        }
    }
}
