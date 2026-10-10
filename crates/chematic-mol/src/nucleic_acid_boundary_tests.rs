use crate::nucleic_acid::*;
use serde_json::{Value, json};

fn document() -> Value {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/nucleic_acid_document_contract.json"
    ))
    .unwrap();
    fixture["valid_cases"][0]["document"].clone()
}

#[test]
fn invalid_nucleic_acid_topology_reports_the_source_path() {
    let original = document();
    for (path, bad, diagnostic) in [
        ("/schema", json!("bad"), "schema"),
        ("/atom_ids/1", json!("d1-5"), "duplicate"),
        ("/atom_ids/0", json!(""), "id"),
        ("/strands", json!([]), "strand"),
        ("/strands/0/id", json!(""), "id"),
        ("/strands/0/residues", json!([]), "residue"),
        ("/strands/0/residues/1/id", json!("d1"), "duplicate"),
        ("/strands/0/residues/0/atom_refs", json!([]), "atom_refs"),
        ("/strands/0/residues/0/atom_refs/0", json!("gone"), "gone"),
        (
            "/strands/0/residues/0/atom_refs/1",
            json!("d1-5"),
            "ambiguous",
        ),
        ("/strands/0/residues/0/sugar", json!("ribose"), "sugar"),
        ("/strands/0/residues/0/base", json!("U"), "base"),
        ("/strands/0/residues/0/base", json!("other"), "modification"),
        (
            "/strands/0/residues/0/modification",
            json!("5mC"),
            "modification",
        ),
        (
            "/strands/0/residues/0/modification",
            json!("unknown"),
            "unknown",
        ),
        ("/linkages/0/id", json!("dna-1"), "duplicate"),
        ("/linkages/0/kind", json!("other"), "topology"),
        ("/linkages/0/three_prime_residue_id", json!("gone"), "gone"),
        ("/linkages/0/five_prime_residue_id", json!("gone"), "gone"),
        ("/linkages/0/three_prime_atom_ref", json!("d2-5"), "atom"),
        ("/linkages/0/five_prime_atom_ref", json!("d1-5"), "atom"),
        ("/linkages", json!([]), "linkages"),
        (
            "/strands/0/residues/0/five_prime_linkage",
            json!("dl1"),
            "terminus",
        ),
        (
            "/strands/0/residues/1/three_prime_linkage",
            json!("dl1"),
            "terminus",
        ),
        (
            "/strands/0/residues/0/three_prime_linkage",
            Value::Null,
            "linkage",
        ),
        (
            "/strands/0/residues/1/five_prime_linkage",
            json!("gone"),
            "linkage",
        ),
    ] {
        let mut bad_document = original.clone();
        if let Some(slot) = bad_document.pointer_mut(path) {
            *slot = bad;
        } else {
            let (parent, key) = path.rsplit_once('/').unwrap();
            bad_document
                .pointer_mut(parent)
                .unwrap()
                .as_object_mut()
                .unwrap()
                .insert(key.into(), bad);
        }
        let error = NucleicAcidDocument::from_json(&bad_document).unwrap_err();
        assert!(error.to_string().contains(diagnostic), "{path}: {error}");
        assert_eq!(document(), original);
    }
}

#[test]
fn each_nucleic_acid_resource_limit_is_enforced() {
    let original = NucleicAcidDocument::from_json(&document()).unwrap();
    let defaults = NucleicAcidLimits::default();
    for limits in [
        NucleicAcidLimits {
            max_strands: 0,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_residues: 1,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_atom_ids: 3,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_atom_refs_per_residue: 1,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_linkages: 0,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_annotations: 0,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_annotation_bytes: 0,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_id_bytes: 1,
            ..defaults.clone()
        },
    ] {
        let error = original.validate_with_limits(&limits).unwrap_err();
        assert!(matches!(
            error,
            NucleicAcidError::ResourceLimit { .. } | NucleicAcidError::InvalidId { .. }
        ));
        assert!(!error.to_string().is_empty());
    }
}

#[test]
fn nucleic_acid_diagnostics_keep_machine_codes_and_source_paths() {
    let limits = NucleicAcidLimits::default();
    let original = document();
    let mut cases = Vec::new();
    for (path, value, code, error_path) in [
        ("/schema", json!("bad"), "invalid_schema", "/"),
        ("/atom_ids/0", json!(""), "invalid_id", "/atom_ids/0"),
        ("/atom_ids/1", json!("d1-5"), "duplicate_id", "/atom_ids/1"),
        (
            "/strands/0/residues/0/atom_refs/0",
            json!("gone"),
            "missing_reference",
            "/strands/0/residues/0/atom_refs/0",
        ),
        (
            "/strands/0/residues/0/sugar",
            json!("ribose"),
            "invalid_residue_identity",
            "/strands/0/residues/0/sugar",
        ),
    ] {
        let mut bad = original.clone();
        *bad.pointer_mut(path).unwrap() = value;
        cases.push((bad, code, error_path));
    }
    let mut duplicate_strand = original.clone();
    duplicate_strand["strands"]
        .as_array_mut()
        .unwrap()
        .push(original["strands"][0].clone());
    cases.push((duplicate_strand, "duplicate_id", "/strands/1/id"));
    let mut unassigned = original.clone();
    unassigned["atom_ids"]
        .as_array_mut()
        .unwrap()
        .push(json!("unassigned"));
    cases.push((unassigned, "ambiguous_atom_mapping", "/atom_ids"));
    let mut missing_linkage = original.clone();
    missing_linkage["strands"][0]["residues"][0]["three_prime_linkage"] = json!("gone");
    missing_linkage["strands"][0]["residues"][1]["five_prime_linkage"] = json!("gone");
    cases.push((
        missing_linkage,
        "missing_reference",
        "/strands/0/residues/0/three_prime_linkage",
    ));
    let mut reversed = original.clone();
    reversed["linkages"][0]["three_prime_residue_id"] = json!("d2");
    reversed["linkages"][0]["five_prime_residue_id"] = json!("d1");
    reversed["linkages"][0]["three_prime_atom_ref"] = json!("d2-3");
    reversed["linkages"][0]["five_prime_atom_ref"] = json!("d1-5");
    cases.push((reversed, "unsupported_linkage_topology", "/linkages/dl1"));
    for (value, code, path) in cases {
        let input = value.to_string();
        let error = NucleicAcidDocument::from_json_str(&input).unwrap_err();
        assert_eq!(error.code(), code);
        assert_eq!(error.path(), path);
        assert_eq!(
            validate_nucleic_acid_json(&input, &limits)["error"],
            error.to_json()
        );
        assert!(error.to_json()["message"].as_str().unwrap().len() > 5);
    }
    let invalid = NucleicAcidDocument::from_json_str("{").unwrap_err();
    assert_eq!(invalid.code(), "invalid_json");
    assert_eq!(invalid.path(), "/");
    let parsed = NucleicAcidDocument::from_json_str(&original.to_string()).unwrap();
    parsed.validate().unwrap();
    assert_eq!(parsed.to_json(), original);
}

#[test]
fn failed_nucleic_acid_edits_preserve_the_document_and_report_missing_targets() {
    let limits = NucleicAcidLimits::default();
    let original = NucleicAcidDocument::from_json(&document()).unwrap();
    let before = original.to_json();
    for kind in ["strand", "residue", "linkage"] {
        for action in ["set_annotation", "remove_annotation"] {
            let mut command =
                json!({"kind": action, "target": {"kind":kind,"id":"gone"}, "key":"label"});
            if action == "set_annotation" {
                command["value"] = json!("test");
            }
            let error = original
                .apply_json_command_with_limits(&command, &limits)
                .unwrap_err();
            assert_eq!(error.code(), "missing_target");
            assert_eq!(error.path(), "/command/target/id");
            assert_eq!(original.to_json(), before);
        }
    }
    let command =
        json!({"kind":"set_residue_identity","residue_id":"gone","base":"A","sugar":"deoxyribose"});
    let error = original
        .apply_json_command_with_limits(&command, &limits)
        .unwrap_err();
    assert_eq!(error.code(), "missing_target");
    assert_eq!(error.path(), "/command/residue_id");
    for command in ["{", "{\"kind\":\"unknown\"}"] {
        let result = apply_nucleic_acid_json_command(&before.to_string(), command, &limits);
        assert_eq!(result["ok"], false);
        assert_eq!(result["error"]["code"], "invalid_command");
        assert_eq!(result["error"]["path"], "/");
    }
    assert_eq!(original.to_json(), before);
}
