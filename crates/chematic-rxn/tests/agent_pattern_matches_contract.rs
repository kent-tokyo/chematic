use chematic_rxn::query::get_reaction_smarts_matches;
use chematic_rxn::{parse_reaction, parse_reaction_query};

#[test]
fn agent_pattern_results_can_be_retrieved_by_alternative() {
    let reaction = parse_reaction("CC>O>CC").unwrap();
    let query = parse_reaction_query("[#6]>[#8]|[#7]>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&reaction, &query);
    let oxygen = matches.agent_matches.get_pattern_matches(0).unwrap();
    assert_eq!(oxygen.len(), 1);
    assert_eq!(oxygen[0].molecule_index, 0);
    assert_eq!(oxygen[0].pattern_index, 0);
    assert!(
        matches
            .agent_matches
            .get_pattern_matches(1)
            .unwrap()
            .is_empty()
    );
    assert!(matches.agent_matches.get_pattern_matches(2).is_none());
}
