use chematic_rxn::query::parse_reaction_smarts;

#[test]
fn single_arrow_places_atom_maps_on_the_correct_side() {
    let pattern = parse_reaction_smarts("[C:1][O:2]>[C:1]=[O:2]").unwrap();
    assert_eq!(pattern.reactant_patterns.len(), 1);
    assert!(pattern.agent_patterns.is_empty());
    assert_eq!(pattern.product_patterns.len(), 1);
    let maps = pattern.map_number_info;
    assert_eq!(maps.reactant_maps, [1, 2].into_iter().collect());
    assert!(maps.agent_maps.is_empty());
    assert_eq!(maps.product_maps, [1, 2].into_iter().collect());
    assert_eq!(maps.all_map_numbers, [1, 2].into_iter().collect());
}
