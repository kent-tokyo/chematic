use super::*;

#[test]
fn validates_r_group_references_and_serializes_schema() {
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        r_groups: vec![RGroupDefinition {
            id: "r1".into(),
            attachment_atoms: vec![AtomRef {
                atom_id: "a1".into(),
            }],
            alternatives: vec!["[*]C".into()],
            selected_alternative: Some(0),
            nested_groups: vec![],
        }],
        ..Default::default()
    };
    model.validate().unwrap();
    assert_eq!(model.to_json()["schema"], "chematic.semantic.v1");
}

#[test]
fn typed_s_group_linkage_round_trips_and_validates_references() {
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        s_groups: vec![SemanticSGroup {
            id: "sg1".into(),
            kind: "SRU".into(),
            member_atoms: vec![AtomRef {
                atom_id: "a1".into(),
            }],
            repeat_unit_id: None,
            linkage: Some(PolymerLinkage {
                id: "link1".into(),
                left_atom: AtomRef {
                    atom_id: "a1".into(),
                },
                right_atom: AtomRef {
                    atom_id: "a2".into(),
                },
                bond_order: 1,
            }),
        }],
        ..Default::default()
    };
    model.validate().unwrap();
    let decoded = SemanticModel::from_json(&model.to_json()).unwrap();
    assert_eq!(decoded, model);

    let invalid = SemanticModel {
        s_groups: vec![SemanticSGroup {
            id: "sg1".into(),
            kind: "SRU".into(),
            member_atoms: vec![AtomRef {
                atom_id: "missing".into(),
            }],
            repeat_unit_id: None,
            linkage: None,
        }],
        ..model.clone()
    };
    assert!(matches!(
        invalid.validate(),
        Err(SemanticError::MissingAtom(id)) if id == "missing"
    ));

    let renamed = model
        .apply_json_command(&serde_json::json!({
            "s_group_id": "sg1",
            "kind": "CROSSLINK"
        }))
        .unwrap();
    assert_eq!(renamed.s_groups[0].kind, "CROSSLINK");
}

#[test]
fn rejects_unimplemented_biomolecule_json_explicitly() {
    let input = serde_json::json!({
        "schema": "chematic.semantic.v1",
        "atom_ids": [],
        "bond_ids": [],
        "r_groups": [],
        "polymer_units": [],
        "biomolecules": [],
        "extensions": {}
    });
    assert!(matches!(
        SemanticModel::from_json(&input),
        Err(SemanticError::Unsupported { construct, .. }) if construct == "biomolecule"
    ));
}

#[test]
fn nested_markush_round_trips_selects_and_expands_deterministically() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        r_groups: vec![RGroupDefinition {
            id: "r1".into(),
            attachment_atoms: vec![AtomRef {
                atom_id: "a1".into(),
            }],
            alternatives: vec!["[*]O".into()],
            selected_alternative: Some(0),
            nested_groups: vec![RGroupDefinition {
                id: "r1a".into(),
                attachment_atoms: vec![AtomRef {
                    atom_id: "a2".into(),
                }],
                alternatives: vec!["[*]N".into()],
                selected_alternative: Some(0),
                nested_groups: vec![],
            }],
        }],
        ..Default::default()
    };
    let roundtrip = SemanticModel::from_json(&model.to_json()).unwrap();
    let expanded = roundtrip.expand(&base).unwrap();
    assert_eq!(expanded.molecule.atom_count(), 4);
    assert_eq!(expanded.source_to_expanded["r1"].len(), 1);
    assert_eq!(expanded.source_to_expanded["r1a"].len(), 1);
    let contracted = expanded.contract();
    assert_eq!(contracted.atom_count(), base.atom_count());
    assert_eq!(contracted.bond_count(), base.bond_count());
    assert_eq!(
        chematic_smiles::canonical_smiles(&contracted),
        chematic_smiles::canonical_smiles(&base)
    );
}

#[test]
fn rejects_dangling_and_ambiguous_constructs() {
    let model = SemanticModel {
        atom_ids: vec!["a1".into()],
        polymer_units: vec![PolymerRepeatUnit {
            id: "p1".into(),
            attachment_atoms: vec![AtomRef {
                atom_id: "missing".into(),
            }],
            end_groups: vec![],
            end_group_definitions: vec![],
            repeat_count: None,
            repeat_smiles: None,
            repeat_endpoint_atoms: None,
        }],
        ..Default::default()
    };
    assert!(matches!(
        model.validate(),
        Err(SemanticError::AmbiguousAttachment(_))
    ));
}

#[test]
fn command_selects_explicit_r_group_and_expands_with_mapping() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        r_groups: vec![RGroupDefinition {
            id: "r1".into(),
            attachment_atoms: vec![AtomRef {
                atom_id: "a2".into(),
            }],
            alternatives: vec!["[*]O".into()],
            selected_alternative: None,
            nested_groups: vec![],
        }],
        ..Default::default()
    };
    let selected = model
        .apply(&SemanticCommand::SelectRGroupAlternative {
            group_id: "r1".into(),
            alternative: 0,
        })
        .unwrap();
    let expanded = selected.expand(&base).unwrap();
    assert_eq!(expanded.molecule.atom_count(), 3);
    assert_eq!(expanded.source_to_expanded["r1"].len(), 1);
}

#[test]
fn command_replaces_allowed_r_group_set_and_clears_selection() {
    let model = SemanticModel {
        atom_ids: vec!["a1".into()],
        r_groups: vec![RGroupDefinition {
            id: "r1".into(),
            attachment_atoms: vec![AtomRef {
                atom_id: "a1".into(),
            }],
            alternatives: vec!["[*]O".into()],
            selected_alternative: Some(0),
            nested_groups: vec![],
        }],
        ..Default::default()
    };
    let changed = model
        .apply_json_command(&serde_json::json!({
            "replace_group_id": "r1",
            "alternatives": ["[*]N", "[*]Cl"]
        }))
        .unwrap();
    assert_eq!(changed.r_groups[0].alternatives, ["[*]N", "[*]Cl"]);
    assert_eq!(changed.r_groups[0].selected_alternative, None);
    assert!(matches!(
        model.apply_json_command(&serde_json::json!({
            "replace_group_id": "r1",
            "alternatives": []
        })),
        Err(SemanticError::InvalidAlternative { .. })
    ));
}

#[test]
fn command_clears_r_group_selection_for_lossless_contraction() {
    let model = SemanticModel {
        atom_ids: vec!["a1".into()],
        r_groups: vec![RGroupDefinition {
            id: "r1".into(),
            attachment_atoms: vec![AtomRef {
                atom_id: "a1".into(),
            }],
            alternatives: vec!["[*]O".into()],
            selected_alternative: Some(0),
            nested_groups: vec![],
        }],
        ..Default::default()
    };
    let contracted = model
        .apply_json_command(&serde_json::json!({"clear_group_id": "r1"}))
        .unwrap();
    assert_eq!(contracted.r_groups[0].selected_alternative, None);
    assert!(matches!(
        contracted.expand(&chematic_smiles::parse("C").unwrap()),
        Err(SemanticError::Unsupported { .. })
    ));
}

#[test]
fn command_selects_polymer_repeat_count_from_unselected_json_state() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        polymer_units: vec![PolymerRepeatUnit {
            id: "p1".into(),
            attachment_atoms: vec![
                AtomRef {
                    atom_id: "a1".into(),
                },
                AtomRef {
                    atom_id: "a2".into(),
                },
            ],
            end_groups: vec![],
            end_group_definitions: vec![],
            repeat_count: None,
            repeat_smiles: Some("[*]CC[*]".into()),
            repeat_endpoint_atoms: None,
        }],
        ..Default::default()
    };
    let model = SemanticModel::from_json(&model.to_json()).unwrap();
    model.validate().unwrap();
    assert!(matches!(
        model.expand(&base),
        Err(SemanticError::Unsupported { .. })
    ));

    let selected = model
        .apply_json_command(&serde_json::json!({
            "unit_id": "p1",
            "repeat_count": 3
        }))
        .unwrap();
    assert_eq!(selected.polymer_units[0].repeat_count, Some(3));
    let expanded = selected.expand(&base).unwrap();
    assert_eq!(expanded.molecule.atom_count(), 8);
    assert_eq!(expanded.source_to_expanded["p1"].len(), 6);
}

#[test]
fn command_clears_polymer_repeat_count_for_lossless_contraction() {
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        polymer_units: vec![PolymerRepeatUnit {
            id: "p1".into(),
            attachment_atoms: vec![
                AtomRef {
                    atom_id: "a1".into(),
                },
                AtomRef {
                    atom_id: "a2".into(),
                },
            ],
            end_groups: vec![],
            end_group_definitions: vec![],
            repeat_count: Some(2),
            repeat_smiles: Some("[*]CC[*]".into()),
            repeat_endpoint_atoms: None,
        }],
        ..Default::default()
    };
    let contracted = model
        .apply_json_command(&serde_json::json!({"clear_unit_id": "p1"}))
        .unwrap();
    assert_eq!(contracted.polymer_units[0].repeat_count, None);
    assert!(matches!(
        contracted.expand(&chematic_smiles::parse("CC").unwrap()),
        Err(SemanticError::Unsupported { .. })
    ));
}

#[test]
fn rejects_json_repeat_count_that_exceeds_u32() {
    let value = serde_json::json!({
        "schema": "chematic.semantic.v1",
        "atom_ids": ["a1", "a2"],
        "bond_ids": [],
        "r_groups": [],
        "polymer_units": [{
            "id": "p1",
            "attachment_atoms": ["a1", "a2"],
            "end_groups": [],
            "repeat_count": 4294967296u64,
            "repeat_smiles": "[*]CC[*]"
        }]
    });
    assert!(matches!(
        SemanticModel::from_json(&value),
        Err(SemanticError::InvalidJson(message))
            if message.contains("repeat_count")
    ));
}

#[test]
fn rejects_json_repeat_endpoint_that_exceeds_u32() {
    let value = serde_json::json!({
        "schema": "chematic.semantic.v1",
        "atom_ids": ["a1", "a2"],
        "bond_ids": [],
        "r_groups": [],
        "polymer_units": [{
            "id": "p1",
            "attachment_atoms": ["a1", "a2"],
            "end_groups": [],
            "repeat_smiles": "CC",
            "repeat_endpoint_atoms": [4294967296u64, 1]
        }]
    });
    assert!(matches!(
        SemanticModel::from_json(&value),
        Err(SemanticError::InvalidJson(message))
            if message.contains("endpoint indices")
    ));
}

#[test]
fn expands_explicit_two_ended_polymer_repeat() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        polymer_units: vec![PolymerRepeatUnit {
            id: "p1".into(),
            attachment_atoms: vec![
                AtomRef {
                    atom_id: "a1".into(),
                },
                AtomRef {
                    atom_id: "a2".into(),
                },
            ],
            end_groups: vec![],
            end_group_definitions: vec![],
            repeat_count: Some(2),
            repeat_smiles: Some("[*]CC[*]".into()),
            repeat_endpoint_atoms: None,
        }],
        ..Default::default()
    };
    let expanded = model.expand(&base).unwrap();
    assert_eq!(expanded.molecule.atom_count(), 6);
    assert_eq!(expanded.source_to_expanded["p1"].len(), 4);
}

#[test]
fn expands_explicit_polymer_endpoints_without_placeholders() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        polymer_units: vec![PolymerRepeatUnit {
            id: "p1".into(),
            attachment_atoms: vec![
                AtomRef {
                    atom_id: "a1".into(),
                },
                AtomRef {
                    atom_id: "a2".into(),
                },
            ],
            end_groups: vec![],
            end_group_definitions: vec![],
            repeat_count: Some(2),
            repeat_smiles: Some("CCO".into()),
            repeat_endpoint_atoms: Some([0, 2]),
        }],
        ..Default::default()
    };
    let expanded = model.expand(&base).unwrap();
    assert_eq!(expanded.molecule.atom_count(), 8);
    assert_eq!(expanded.source_to_expanded["p1"].len(), 6);
}

#[test]
fn json_contract_round_trips_and_expands_selected_markush() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        r_groups: vec![RGroupDefinition {
            id: "r1".into(),
            attachment_atoms: vec![AtomRef {
                atom_id: "a2".into(),
            }],
            alternatives: vec!["[*]O".into()],
            selected_alternative: Some(0),
            nested_groups: vec![],
        }],
        ..Default::default()
    };
    let decoded = SemanticModel::from_json(&model.to_json()).unwrap();
    let expanded = decoded.expand(&base).unwrap();
    assert_eq!(expanded.molecule.atom_count(), 3);
    assert_eq!(expanded.to_json()["source_to_expanded"]["r1"][0], 2);
    assert_eq!(expanded.to_json()["contracted_smiles"], "CC");
}

#[test]
fn expands_markush_with_multiple_attachment_points_in_marker_order() {
    let base = chematic_smiles::parse("CCCC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into(), "a3".into(), "a4".into()],
        r_groups: vec![RGroupDefinition {
            id: "r1".into(),
            attachment_atoms: vec![
                AtomRef {
                    atom_id: "a1".into(),
                },
                AtomRef {
                    atom_id: "a4".into(),
                },
            ],
            alternatives: vec!["[*]O[*]".into()],
            selected_alternative: Some(0),
            nested_groups: vec![],
        }],
        ..Default::default()
    };
    let expanded = model.expand(&base).unwrap();
    assert_eq!(expanded.molecule.atom_count(), 5);
    assert_eq!(expanded.source_to_expanded["r1"].len(), 1);
}

#[test]
fn expands_polymer_end_groups_with_stable_mapping() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        polymer_units: vec![PolymerRepeatUnit {
            id: "p1".into(),
            attachment_atoms: vec![
                AtomRef {
                    atom_id: "a1".into(),
                },
                AtomRef {
                    atom_id: "a2".into(),
                },
            ],
            end_groups: vec!["[*]O".into(), "[*]N".into()],
            end_group_definitions: vec![],
            repeat_count: Some(1),
            repeat_smiles: Some("[*]CC[*]".into()),
            repeat_endpoint_atoms: None,
        }],
        ..Default::default()
    };
    let expanded = model.expand(&base).unwrap();
    assert_eq!(expanded.molecule.atom_count(), 6);
    assert_eq!(expanded.source_to_expanded["p1"].len(), 4);
    assert_eq!(expanded.source_to_expanded["p1.end_group_left"].len(), 1);
    assert_eq!(expanded.source_to_expanded["p1.end_group_right"].len(), 1);
}

#[test]
fn expands_typed_polymer_end_groups_with_named_provenance() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        polymer_units: vec![PolymerRepeatUnit {
            id: "p1".into(),
            attachment_atoms: vec![
                AtomRef {
                    atom_id: "a1".into(),
                },
                AtomRef {
                    atom_id: "a2".into(),
                },
            ],
            end_groups: vec![],
            end_group_definitions: vec![
                PolymerEndGroup {
                    id: "cap-left".into(),
                    smiles: "[*]O".into(),
                },
                PolymerEndGroup {
                    id: "cap-right".into(),
                    smiles: "[*]N".into(),
                },
            ],
            repeat_count: Some(1),
            repeat_smiles: Some("[*]CC[*]".into()),
            repeat_endpoint_atoms: None,
        }],
        ..Default::default()
    };
    let decoded = SemanticModel::from_json(&model.to_json()).unwrap();
    let expanded = decoded.expand(&base).unwrap();
    assert_eq!(expanded.source_to_expanded["cap-left"].len(), 1);
    assert_eq!(expanded.source_to_expanded["cap-right"].len(), 1);
}

#[test]
fn rejects_zero_repeat_count() {
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        polymer_units: vec![PolymerRepeatUnit {
            id: "p1".into(),
            attachment_atoms: vec![
                AtomRef {
                    atom_id: "a1".into(),
                },
                AtomRef {
                    atom_id: "a2".into(),
                },
            ],
            end_groups: vec![],
            end_group_definitions: vec![],
            repeat_count: Some(0),
            repeat_smiles: Some("[*]CC[*]".into()),
            repeat_endpoint_atoms: None,
        }],
        ..Default::default()
    };
    assert!(matches!(
        model.validate(),
        Err(SemanticError::InvalidExpansion { .. })
    ));
}

#[test]
fn rejects_repeat_expansion_over_explicit_budget_before_allocating() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        polymer_units: vec![PolymerRepeatUnit {
            id: "p1".into(),
            attachment_atoms: vec![
                AtomRef {
                    atom_id: "a1".into(),
                },
                AtomRef {
                    atom_id: "a2".into(),
                },
            ],
            end_groups: vec![],
            end_group_definitions: vec![],
            repeat_count: Some(11),
            repeat_smiles: Some("[*]CC[*]".into()),
            repeat_endpoint_atoms: None,
        }],
        ..Default::default()
    };
    let limits = SemanticExpansionLimits {
        max_atoms: 100,
        max_repeat_count: 10,
    };
    assert!(matches!(
        model.expand_with_limits(&base, &limits),
        Err(SemanticError::ExpansionLimit {
            resource: "repeats",
            requested: 11,
            limit: 10,
            ..
        })
    ));
}

#[test]
fn rejects_atom_expansion_over_budget_before_mutation() {
    let base = chematic_smiles::parse("CC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into()],
        r_groups: vec![RGroupDefinition {
            id: "r1".into(),
            attachment_atoms: vec![AtomRef {
                atom_id: "a2".into(),
            }],
            alternatives: vec!["[*]CCCC".into()],
            selected_alternative: Some(0),
            nested_groups: vec![],
        }],
        ..Default::default()
    };
    let limits = SemanticExpansionLimits {
        max_atoms: 4,
        max_repeat_count: 10,
    };
    assert!(matches!(
        model.expand_with_limits(&base, &limits),
        Err(SemanticError::ExpansionLimit {
            resource: "atoms",
            requested: 6,
            limit: 4,
            ..
        })
    ));
}

#[test]
fn rejects_base_molecule_over_budget_without_semantic_additions() {
    let base = chematic_smiles::parse("CCCC").unwrap();
    let model = SemanticModel {
        atom_ids: vec!["a1".into(), "a2".into(), "a3".into(), "a4".into()],
        ..Default::default()
    };
    let limits = SemanticExpansionLimits {
        max_atoms: 3,
        max_repeat_count: 10,
    };
    assert!(matches!(
        model.expand_with_limits(&base, &limits),
        Err(SemanticError::ExpansionLimit {
            resource: "atoms",
            requested: 4,
            limit: 3,
            ..
        })
    ));
}
