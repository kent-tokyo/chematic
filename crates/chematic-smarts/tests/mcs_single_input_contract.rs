use chematic_smarts::{
    MatchConfig, McsConfig, McsOutcome, find_matches_with_config, find_mcs_with_config_checked,
};

#[test]
fn single_input_mcs_respects_zero_search_budget() {
    let molecule = chematic_smiles::parse("CCO").unwrap();
    let config = McsConfig {
        max_search_nodes: Some(0),
        ..Default::default()
    };
    let result = find_mcs_with_config_checked(&[&molecule], &config);
    assert!(matches!(&result, McsOutcome::ResourceLimited(_)));
    assert_eq!(result.into_query().atom_count(), 0);
}

#[test]
fn single_input_mcs_preserves_tetrahedral_stereochemistry_when_requested() {
    let template = chematic_smiles::parse("F[C@](Cl)(Br)I").unwrap();
    let opposite = chematic_smiles::parse("F[C@@](Cl)(Br)I").unwrap();
    let config = McsConfig {
        match_chiral_tag: true,
        ..Default::default()
    };
    let result = find_mcs_with_config_checked(&[&template], &config);
    assert!(matches!(&result, McsOutcome::Exhaustive(_)));
    let query = result.into_query();
    assert_eq!(query.atom_count(), template.atom_count());
    let match_config = MatchConfig {
        use_chirality: true,
        ..Default::default()
    };
    assert!(!find_matches_with_config(&query, &template, &match_config).is_empty());
    assert!(find_matches_with_config(&query, &opposite, &match_config).is_empty());
}
