use chematic_smiles::{
    CanonicalizationError, CanonicalizationLimits, canonical_smiles, canonical_smiles_with_limits,
    parse,
};

#[test]
fn canonical_budget_failure_reports_search_effort_without_partial_result() {
    let molecule = parse("c1ccccc1").unwrap();
    let limit = CanonicalizationLimits {
        max_search_nodes: Some(0),
        max_automorphism_tests: None,
    };
    let error = canonical_smiles_with_limits(&molecule, &limit).unwrap_err();
    let CanonicalizationError::SearchBudgetExceeded {
        nodes_visited,
        automorphism_tests,
    } = error
    else {
        panic!("unexpected canonicalization failure");
    };
    let message = error.to_string();
    assert!(message.contains(&format!("{nodes_visited} nodes")));
    assert!(message.contains(&format!("{automorphism_tests} automorphism tests")));
    assert_eq!(
        canonical_smiles_with_limits(&molecule, &CanonicalizationLimits::unbounded()).unwrap(),
        canonical_smiles(&molecule)
    );
}
