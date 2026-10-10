use crate::semantic::*;
use serde_json::{Value, json};

fn model() -> SemanticModel {
    SemanticModel::from_json(&json!({
        "schema":"chematic.semantic.v1", "atom_ids":["a1", "a2"], "bond_ids":[],
        "r_groups":[{"id":"r1", "attachment_atoms":["a1"], "alternatives":["[*]C"], "selected_alternative":0}],
        "polymer_units":[{"id":"p1", "attachment_atoms":["a1", "a2"], "end_groups":[], "repeat_count":2, "repeat_smiles":"[*]CC[*]"}],
        "s_groups":[{"id":"s1", "kind":"SRU", "member_atoms":["a1"], "repeat_unit_id":"p1",
            "linkage":{"id":"l1", "left_atom":"a1", "right_atom":"a2", "bond_order":1}}]
    })).unwrap()
}

#[test]
fn malformed_semantic_json_has_field_specific_diagnostics() {
    let original = model().to_json();
    for (pointer, bad, diagnostic) in [
        ("/schema", json!("other"), "schema"),
        ("/atom_ids", json!(false), "atom_ids"),
        ("/bond_ids", json!([5]), "bond_ids"),
        ("/r_groups", json!(false), "r_groups"),
        (
            "/r_groups/0/selected_alternative",
            json!(-1),
            "selected_alternative",
        ),
        ("/r_groups/0/nested_groups", json!(true), "nested_groups"),
        ("/r_groups/0/alternatives", json!([true]), "alternatives"),
        ("/polymer_units", json!(true), "polymer_units"),
        (
            "/polymer_units/0/repeat_count",
            json!(4294967296u64),
            "repeat_count",
        ),
        (
            "/polymer_units/0/end_group_definitions",
            json!(true),
            "end_group_definitions",
        ),
        (
            "/polymer_units/0/repeat_endpoint_atoms",
            json!(false),
            "repeat_endpoint_atoms",
        ),
        (
            "/polymer_units/0/repeat_endpoint_atoms",
            json!([0]),
            "two indices",
        ),
        (
            "/polymer_units/0/repeat_endpoint_atoms",
            json!([-1, 2]),
            "u32",
        ),
        ("/s_groups", json!(false), "s_groups"),
    ] {
        let mut value = original.clone();
        if let Some(slot) = value.pointer_mut(pointer) {
            *slot = bad;
        } else {
            let (parent, key) = pointer.rsplit_once('/').unwrap();
            value
                .pointer_mut(parent)
                .unwrap()
                .as_object_mut()
                .unwrap()
                .insert(key.into(), bad);
        }
        let error = SemanticModel::from_json(&value).unwrap_err();
        assert!(error.to_string().contains(diagnostic), "{pointer}: {error}");
    }
}

#[test]
fn semantic_references_and_polymer_topology_are_validated() {
    let original = model();
    let mut cases: Vec<(SemanticModel, &str)> = Vec::new();
    let mut bad = original.clone();
    bad.bond_ids.push("a1".into());
    cases.push((bad, "duplicate"));
    let mut bad = original.clone();
    bad.polymer_units[0].id = "r1".into();
    cases.push((bad, "duplicate"));
    let mut bad = original.clone();
    bad.polymer_units[0].attachment_atoms.clear();
    cases.push((bad, "attachment"));
    let mut bad = original.clone();
    bad.polymer_units[0].repeat_count = Some(0);
    cases.push((bad, "repeat count"));
    let mut bad = original.clone();
    bad.polymer_units[0].end_groups = vec!["[*]O".into()];
    cases.push((bad, "end_groups"));
    let mut bad = original.clone();
    bad.polymer_units[0].end_groups = vec!["".into(), "[*]O".into()];
    cases.push((bad, "empty"));
    let mut bad = original.clone();
    bad.polymer_units[0].repeat_smiles = Some("CC".into());
    cases.push((bad, "endpoint"));
    let mut bad = original.clone();
    bad.polymer_units[0].repeat_endpoint_atoms = Some([0, 0]);
    cases.push((bad, "distinct"));
    let mut bad = original.clone();
    bad.polymer_units[0].attachment_atoms[0].atom_id = "gone".into();
    cases.push((bad, "gone"));
    let mut bad = original.clone();
    bad.s_groups[0].id = "r1".into();
    cases.push((bad, "duplicate"));
    let mut bad = original.clone();
    bad.s_groups[0].kind.clear();
    cases.push((bad, "empty kind"));
    let mut bad = original.clone();
    bad.s_groups[0].member_atoms.clear();
    cases.push((bad, "member atom"));
    let mut bad = original.clone();
    bad.s_groups[0].member_atoms[0].atom_id = "gone".into();
    cases.push((bad, "gone"));
    let mut bad = original.clone();
    bad.s_groups[0].repeat_unit_id = Some("gone".into());
    cases.push((bad, "gone"));
    let mut bad = original.clone();
    bad.s_groups[0].linkage.as_mut().unwrap().id.clear();
    cases.push((bad, "duplicate"));
    let mut bad = original.clone();
    bad.s_groups[0].linkage.as_mut().unwrap().left_atom.atom_id = "gone".into();
    cases.push((bad, "gone"));
    let mut bad = original.clone();
    bad.s_groups[0].linkage.as_mut().unwrap().right_atom.atom_id = "a1".into();
    cases.push((bad, "attachment"));
    let mut bad = original.clone();
    bad.s_groups[0].linkage.as_mut().unwrap().bond_order = 8;
    cases.push((bad, "bond order"));
    for (bad, diagnostic) in cases {
        let error = bad.validate().unwrap_err();
        assert!(error.to_string().contains(diagnostic), "{error}");
    }
    assert_eq!(model(), original);
}

#[test]
fn edits_of_missing_objects_are_atomic() {
    let original = model();
    for command in [
        SemanticCommand::SelectRGroupAlternative {
            group_id: "gone".into(),
            alternative: 0,
        },
        SemanticCommand::ReplaceRGroupAlternatives {
            group_id: "gone".into(),
            alternatives: vec![],
        },
        SemanticCommand::ClearRGroupAlternative {
            group_id: "gone".into(),
        },
        SemanticCommand::SetSGroupKind {
            group_id: "gone".into(),
            kind: "SRU".into(),
        },
        SemanticCommand::ClearPolymerRepeatCount {
            unit_id: "gone".into(),
        },
        SemanticCommand::SetPolymerRepeatCount {
            unit_id: "gone".into(),
            repeat_count: 2,
        },
    ] {
        assert!(
            original
                .apply(&command)
                .unwrap_err()
                .to_string()
                .contains("gone")
        );
        assert_eq!(original, model());
    }
    let error = original
        .apply(&SemanticCommand::SetPolymerRepeatCount {
            unit_id: "p1".into(),
            repeat_count: 0,
        })
        .unwrap_err();
    assert!(error.to_string().contains("greater than zero"));
    for command in [
        json!({"replace_group_id":"r1"}),
        json!({"replace_group_id":"r1", "alternatives":[7]}),
        json!({"s_group_id":"s1"}),
        Value::Null,
    ] {
        assert!(matches!(
            original.apply_json_command(&command),
            Err(SemanticError::InvalidJson(_))
        ));
    }
}

#[test]
fn expansion_checks_resource_limits_before_modifying_the_base() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        ..Default::default()
    };
    for limits in [
        SemanticExpansionLimits {
            max_atoms: 0,
            max_repeat_count: 1,
        },
        SemanticExpansionLimits {
            max_atoms: 2,
            max_repeat_count: 0,
        },
    ] {
        assert!(
            model
                .expand_with_limits(&base, &limits)
                .err()
                .unwrap()
                .to_string()
                .contains("greater than zero")
        );
    }
    let error = model
        .expand_with_limits(
            &base,
            &SemanticExpansionLimits {
                max_atoms: 1,
                max_repeat_count: 1,
            },
        )
        .err()
        .unwrap();
    assert!(matches!(
        error,
        SemanticError::ExpansionLimit {
            resource: "atoms",
            requested: 2,
            limit: 1,
            ..
        }
    ));
    assert!(error.to_string().contains("limit is 1"));
    let wrong = SemanticModel::default();
    assert!(
        wrong
            .expand(&base)
            .err()
            .unwrap()
            .to_string()
            .contains("atom_ids")
    );
    assert_eq!(base.atom_count(), 2);
    assert_eq!(base.bond_count(), 1);
}
