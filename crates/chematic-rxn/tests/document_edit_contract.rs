use chematic_rxn::ReactionDocument;

#[test]
fn edits_update_by_stable_id_and_leave_original_unchanged() {
    let original = ReactionDocument::from_reaction_smiles("CC>>O").unwrap();
    let renamed = original
        .apply_json_edit(r#"{"kind":"set_document_id","id":"doc-42"}"#)
        .unwrap();
    assert_eq!(renamed.id, "doc-42");
    assert_ne!(original.id, renamed.id);
    let with_condition = renamed.apply_json_edit(r#"{"kind":"set_step_condition","step_id":"step-1","key":"temperature","value":"20 C"}"#).unwrap();
    let updated = with_condition.apply_json_edit(r#"{"kind":"set_step_condition","step_id":"step-1","key":"temperature","value":"25 C"}"#).unwrap();
    assert_eq!(with_condition.steps[0].conditions[0].value, "20 C");
    assert_eq!(updated.steps[0].conditions.len(), 1);
    assert_eq!(updated.steps[0].conditions[0].value, "25 C");
    let weighted = updated
        .apply_json_edit(
            r#"{"kind":"set_component_coefficient","component_id":"reactant-1","coefficient":2}"#,
        )
        .unwrap();
    assert_eq!(weighted.steps[0].components[0].coefficient, 2);
    assert_eq!(updated.steps[0].components[0].coefficient, 1);
    let error = weighted
        .apply_json_edit(
            r#"{"kind":"set_step_condition","step_id":"missing","key":"pH","value":"7"}"#,
        )
        .unwrap_err();
    assert!(error.to_string().contains("unknown step id 'missing'"));
    assert_eq!(weighted.steps[0].conditions[0].value, "25 C");
}

#[test]
fn document_validation_rejects_missing_and_duplicate_identifiers() {
    let original = ReactionDocument::from_reaction_smiles("CC>>O").unwrap();
    let mut missing_id = original.clone();
    missing_id.id.clear();
    assert!(
        missing_id
            .validate()
            .unwrap_err()
            .to_string()
            .contains("document id and at least one step")
    );

    let mut no_steps = original.clone();
    no_steps.steps.clear();
    assert!(
        no_steps
            .validate()
            .unwrap_err()
            .to_string()
            .contains("document id and at least one step")
    );

    let mut no_component_id = original.clone();
    no_component_id.steps[0].components[0].id.clear();
    assert!(
        no_component_id
            .validate()
            .unwrap_err()
            .to_string()
            .contains("component IDs must be non-empty")
    );

    let mut duplicate = original;
    duplicate.steps[0].components[1].id = duplicate.steps[0].components[0].id.clone();
    assert!(
        duplicate
            .validate()
            .unwrap_err()
            .to_string()
            .contains("globally unique")
    );
}
