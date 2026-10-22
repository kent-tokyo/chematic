use super::*;
use crate::parse_reaction;

fn rxn(s: &str) -> Reaction {
    crate::reaction::parse_reaction(s).unwrap()
}

#[test]
fn test_parse_reaction_query_basic() {
    let query = parse_reaction_query("[#6]>>[#6]").unwrap();
    assert_eq!(query.reactant_patterns.len(), 1);
    assert_eq!(query.product_patterns.len(), 1);
}

#[test]
fn test_parse_reaction_query_multiple_patterns() {
    let query = parse_reaction_query("[#6]|[#7]>>[#8]|[#9]").unwrap();
    assert_eq!(query.reactant_patterns.len(), 2);
    assert_eq!(query.product_patterns.len(), 2);
}

#[test]
fn test_has_reaction_substructure_match_simple() {
    // Reaction: ethane to ethane (trivial)
    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#6]>>[#6]").unwrap();
    assert!(has_reaction_substructure_match(&rxn, &query));
}

#[test]
fn test_has_reaction_substructure_match_no_match() {
    // Reaction: ethane to ethane
    let rxn = rxn("CC>>CC");
    // Query: looking for nitrogen (not present)
    let query = parse_reaction_query("[#7]>>[#7]").unwrap();
    assert!(!has_reaction_substructure_match(&rxn, &query));
}

#[test]
fn test_has_reaction_substructure_match_product_mismatch() {
    // Reaction: ethane to methane
    let rxn = rxn("CC>>C");
    // Query: looking for ethane in products (not present)
    let query = parse_reaction_query("[#6]>>CC").unwrap();
    assert!(!has_reaction_substructure_match(&rxn, &query));
}

#[test]
fn test_reaction_query_matches_agent_section() {
    let rxn = rxn("CC>[Pd]>CC");
    let query = parse_reaction_query("[#6]>[Pd]>[#6]").unwrap();
    assert!(has_reaction_substructure_match(&rxn, &query));

    let non_matching = parse_reaction_query("[#6]>[Ni]>[#6]").unwrap();
    assert!(!has_reaction_substructure_match(&rxn, &non_matching));
}

#[test]
fn test_reaction_query_agent_alternatives_use_or_semantics() {
    let reaction = rxn("CC>[Ni]>CC");
    let query = parse_reaction_query("[#6]>[Pd]|[Ni]>[#6]").unwrap();
    assert!(has_reaction_substructure_match(&reaction, &query));

    let absent = rxn("CC>[Cu]>CC");
    assert!(!has_reaction_substructure_match(&absent, &query));
}

#[test]
fn test_reaction_query_respects_target_atom_maps_when_present() {
    let mapped = rxn("[C:1][C:2]>>[C:1][C:2]");
    let query = parse_reaction_query("[C:1][C:2]>>[C:1][C:2]").unwrap();
    assert!(has_reaction_substructure_match(&mapped, &query));

    let differently_mapped = rxn("[C:9][C:10]>>[C:9][C:10]");
    assert!(!has_reaction_substructure_match(
        &differently_mapped,
        &query
    ));
}

#[test]
fn test_has_reaction_substructure_match_empty_query() {
    // Empty query should match any reaction (trivial)
    let rxn = rxn("CC>>C");
    let query = parse_reaction_query(">>").unwrap();
    assert!(has_reaction_substructure_match(&rxn, &query));
}

#[test]
fn test_reaction_query_matches_disconnected_components_without_reuse() {
    let query = parse_reaction_query("[#6].[#7]>>[#6].[#7]").unwrap();
    assert!(has_reaction_substructure_match(&rxn("C.N>>C.N"), &query));
    assert!(!has_reaction_substructure_match(&rxn("C>>C"), &query));
}

#[test]
fn test_reaction_query_pipe_alternatives_are_or_on_reactant_and_product_sides() {
    let query = parse_reaction_query("[#6]|[#7]>>[#8]|[#9]").unwrap();
    assert!(has_reaction_substructure_match(&rxn("C>>O"), &query));
    assert!(has_reaction_substructure_match(&rxn("N>>F"), &query));
    assert!(!has_reaction_substructure_match(&rxn("C>>N"), &query));
}

// ===== Phase 1: Agents Section Support =====

#[test]
fn test_parse_reaction_smarts_legacy_format() {
    // Legacy 2-part format should still work
    let pattern = parse_reaction_smarts("[C:1][C:2]>>[C:2][C:1]").unwrap();
    assert_eq!(pattern.reactant_patterns.len(), 1);
    assert_eq!(pattern.agent_patterns.len(), 0);
    assert_eq!(pattern.product_patterns.len(), 1);
}

#[test]
fn test_parse_reaction_smarts_with_agents() {
    // New 3-part format with agents
    let pattern = parse_reaction_smarts("[C:1][C:2]>[Pd]>[C:2][C:1]").unwrap();
    assert_eq!(pattern.reactant_patterns.len(), 1);
    assert_eq!(pattern.agent_patterns.len(), 1);
    assert_eq!(pattern.product_patterns.len(), 1);
}

#[test]
fn test_parse_reaction_smarts_empty_agents() {
    // Format: reactants > (empty) > products
    let pattern = parse_reaction_smarts("[C:1][C:2]>>[C:2][C:1]").unwrap();
    assert_eq!(pattern.agent_patterns.len(), 0);
}

#[test]
fn test_parse_reaction_smarts_multiple_agent_patterns() {
    // Multiple agent patterns with pipe-separated OR logic
    let pattern = parse_reaction_smarts("[C:1]>[Pd]|[Ni]>[C:1]").unwrap();
    assert_eq!(pattern.agent_patterns.len(), 2);
}

#[test]
fn test_extract_map_numbers_basic() {
    let pattern = parse_reaction_smarts("[C:1][C:2]>>[C:2][C:1]").unwrap();
    assert!(pattern.map_number_info.reactant_maps.contains(&1));
    assert!(pattern.map_number_info.reactant_maps.contains(&2));
    assert!(pattern.map_number_info.product_maps.contains(&1));
    assert!(pattern.map_number_info.product_maps.contains(&2));
}

#[test]
fn test_extract_map_numbers_with_agents() {
    let pattern = parse_reaction_smarts("[C:1]>[Pd:3]>[C:1]").unwrap();
    assert!(pattern.map_number_info.reactant_maps.contains(&1));
    assert!(pattern.map_number_info.agent_maps.contains(&3));
    assert!(pattern.map_number_info.product_maps.contains(&1));
}

#[test]
fn test_map_number_validation_missing_in_products() {
    // Map number :1 in reactants but not in products should fail
    let result = parse_reaction_smarts("[C:1][C:2]>>[C:2]");
    assert!(result.is_err());
    if let Err(ReactionQueryError::MapNumberMismatch { map_num, message }) = result {
        assert_eq!(map_num, 1);
        assert!(message.contains("missing from products"));
    } else {
        panic!("expected MapNumberMismatch error");
    }
}

#[test]
fn test_map_number_validation_undefined_in_reactants() {
    // Map number :1 in products but not in reactants should fail
    let result = parse_reaction_smarts("[C:2]>>[C:1][C:2]");
    assert!(result.is_err());
    if let Err(ReactionQueryError::MapNumberMismatch { map_num, message }) = result {
        assert_eq!(map_num, 1);
        assert!(message.contains("not found in reactants"));
    } else {
        panic!("expected MapNumberMismatch error");
    }
}

#[test]
fn test_map_number_validation_consistent() {
    // Valid: all map numbers present in both reactants and products
    let pattern = parse_reaction_smarts("[C:1][C:2][C:3]>>[C:3][C:1][C:2]").unwrap();
    assert_eq!(pattern.map_number_info.reactant_maps.len(), 3);
    assert_eq!(pattern.map_number_info.product_maps.len(), 3);
}

#[test]
fn test_map_number_mismatch_diagnostic_uses_smallest_offending_map() {
    // The source sets are hash-based, so this also protects the error's
    // machine-readable map_num field from process-dependent ordering.
    let missing = parse_reaction_smarts("[C:7][C:2][C:9]>>[C:9]");
    assert!(matches!(
        missing,
        Err(ReactionQueryError::MapNumberMismatch { map_num: 2, .. })
    ));

    let undefined = parse_reaction_smarts("[C:7]>>[C:7][C:9][C:2]");
    assert!(matches!(
        undefined,
        Err(ReactionQueryError::MapNumberMismatch { map_num: 2, .. })
    ));
}

#[test]
fn test_map_numbers_no_maps() {
    // Pattern with no map numbers is valid
    let pattern = parse_reaction_smarts("[#6]>>[#6]").unwrap();
    assert_eq!(pattern.map_number_info.all_map_numbers.len(), 0);
}

#[test]
fn test_map_numbers_large_numbers() {
    // Handle large map numbers (e.g., :100, :999)
    let pattern = parse_reaction_smarts("[C:100][N:999]>>[N:999][C:100]").unwrap();
    assert!(pattern.map_number_info.reactant_maps.contains(&100));
    assert!(pattern.map_number_info.reactant_maps.contains(&999));
}

#[test]
fn test_parse_reaction_smarts_missing_delimiter() {
    // No '>' or '>>' delimiter should error
    let result = parse_reaction_smarts("[C:1][C:2][C:1][C:2]");
    assert!(result.is_err());
}

#[test]
fn test_agents_section_in_map_numbers() {
    // Agents can have their own map numbers (independent of reactant/product validation)
    let pattern = parse_reaction_smarts("[C:1]>[Pd:99]>[C:1]").unwrap();
    assert!(pattern.map_number_info.agent_maps.contains(&99));
    // Agent map :99 is not required in reactants/products
    assert!(!pattern.map_number_info.reactant_maps.contains(&99));
}

// ===== Phase 2: Detailed Match Information =====

#[test]
fn test_get_reaction_smarts_matches_basic() {
    // Simple match: ethane to ethane with carbon pattern
    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#6]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Should have one reactant pattern and one product pattern
    assert_eq!(matches.reactant_matches.pattern_matches.len(), 1);
    assert_eq!(matches.product_matches.pattern_matches.len(), 1);

    // Both patterns should match
    assert!(matches.all_reactants_matched());
    assert!(matches.all_products_matched());
    assert!(matches.is_complete_match);
}

#[test]
fn test_get_reaction_smarts_matches_reports_agent_alternatives() {
    let reaction = rxn("CC>[Ni]>CC");
    let query = parse_reaction_query("[#6]>[Pd]|[Ni]>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&reaction, &query);

    assert!(matches.is_complete_match);
    assert_eq!(matches.agent_matches.pattern_matches.len(), 2);
    assert!(matches.agent_matches.pattern_matches[0].is_empty());
    assert_eq!(matches.agent_matches.pattern_matches[1].len(), 1);
    assert!(matches.agent_matches.any_pattern_matched());
}

#[test]
fn detailed_match_uses_alternative_group_semantics() {
    // The flattened detail vectors contain one unmatched alternative,
    // but the reaction query is satisfied because C|N means either one.
    // The summary helpers must agree with is_complete_match rather than
    // treating every flattened alternative as independently required.
    let reaction = rxn("CC>>CC");
    let query = parse_reaction_query("[#7]|[#6]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&reaction, &query);

    assert!(matches.is_complete_match);
    assert!(matches.all_reactants_matched());
    assert!(matches.all_products_matched());
    assert!(matches.reactants_satisfied);
    assert!(matches.products_satisfied);
    assert!(matches.agents_satisfied);
    assert!(matches.reactant_matches.pattern_matches[0].is_empty());
    assert!(!matches.reactant_matches.pattern_matches[1].is_empty());
}

#[test]
fn test_get_reaction_smarts_matches_multiple_reactants() {
    // Multiple reactants: C + C >> CC
    let rxn = rxn("C.C>>CC");
    let query = parse_reaction_query("[#6]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Reactant pattern should match both molecules
    assert_eq!(matches.reactant_matches.pattern_matches[0].len(), 2);
    // Product pattern should match the combined molecule
    assert_eq!(matches.product_matches.pattern_matches[0].len(), 1);
}

#[test]
fn test_get_reaction_smarts_matches_incomplete() {
    // Query for nitrogen (not present) should have incomplete match
    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#7]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Reactant pattern should not match
    assert!(!matches.all_reactants_matched());
    assert!(!matches.is_complete_match);
}

#[test]
fn test_get_reaction_smarts_matches_empty_patterns() {
    // Empty query should match anything
    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query(">>").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Empty patterns should be considered complete matches
    assert!(matches.is_complete_match);
}

#[test]
fn test_molecule_match_has_correct_indices() {
    // Verify that matched atoms have correct indices
    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#6]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Reactant carbon pattern matches (first match = first atom)
    let reactant_matches = &matches.reactant_matches.pattern_matches[0];
    assert_eq!(reactant_matches.len(), 1);
    assert_eq!(reactant_matches[0].molecule_index, 0);
    assert!(!reactant_matches[0].atom_indices.is_empty()); // At least one carbon matched

    // Product carbon pattern matches
    let product_matches = &matches.product_matches.pattern_matches[0];
    assert_eq!(product_matches.len(), 1);
    assert_eq!(product_matches[0].molecule_index, 0);
    assert!(!product_matches[0].atom_indices.is_empty()); // At least one carbon matched
}

#[test]
fn molecule_match_atom_indices_follow_query_order() {
    let rxn = rxn("CO>>CO");
    let query = parse_reaction_query("[O:2][C:1]>>[O:2][C:1]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Query atom 0 is oxygen (target index 1), followed by carbon
    // (target index 0).  This must not depend on hash-map iteration.
    assert_eq!(
        matches.reactant_matches.pattern_matches[0][0].atom_indices,
        vec![1, 0]
    );
    assert_eq!(
        matches.product_matches.pattern_matches[0][0].atom_indices,
        vec![1, 0]
    );
}

#[test]
fn test_reactant_matches_get_pattern_matches() {
    // Test the get_pattern_matches helper method
    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#6]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Should be able to retrieve matches for a specific pattern
    let reactant_pattern_0 = matches.reactant_matches.get_pattern_matches(0);
    assert!(reactant_pattern_0.is_some());
    assert_eq!(reactant_pattern_0.unwrap().len(), 1);

    // Out of bounds pattern should return None
    let reactant_pattern_1 = matches.reactant_matches.get_pattern_matches(1);
    assert!(reactant_pattern_1.is_none());
}

#[test]
fn test_product_matches_get_pattern_matches() {
    // Test the get_pattern_matches helper method for products
    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#6]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    let product_pattern_0 = matches.product_matches.get_pattern_matches(0);
    assert!(product_pattern_0.is_some());
    assert_eq!(product_pattern_0.unwrap().len(), 1);
}

#[test]
fn test_pattern_matched_helper() {
    // Test the pattern_matched helper method
    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#6]>>[#7]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Reactant pattern should match
    assert!(matches.reactant_matches.pattern_matched(0));

    // Product pattern should not match (looking for nitrogen)
    assert!(!matches.product_matches.pattern_matched(0));
}

#[test]
fn test_multiple_patterns_or_logic() {
    // Test with multiple patterns (OR logic)
    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#6]|[#7]>>[#6]|[#8]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Should have two patterns each (carbon|nitrogen and carbon|oxygen)
    assert_eq!(matches.reactant_matches.pattern_matches.len(), 2);
    assert_eq!(matches.product_matches.pattern_matches.len(), 2);

    // Both reactant patterns should match (we have carbons)
    assert!(matches.reactant_matches.pattern_matched(0)); // [#6]
    assert!(!matches.reactant_matches.pattern_matched(1)); // [#7] - no nitrogen

    // Both product patterns should match
    assert!(matches.product_matches.pattern_matched(0)); // [#6]
    assert!(!matches.product_matches.pattern_matched(1)); // [#8] - no oxygen

    // Overall match is true because each side has one satisfied alternative.
    assert!(matches.is_complete_match);
}

#[test]
fn test_complex_reaction_matches() {
    // Test with more complex reaction
    let rxn = rxn("CC(C)>>CC=C");
    let query = parse_reaction_query("[#6]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    // Both should match
    assert!(matches.all_reactants_matched());
    assert!(matches.all_products_matched());
    assert!(matches.is_complete_match);

    // Reactant has one molecule matching the pattern
    assert_eq!(matches.reactant_matches.pattern_matches[0].len(), 1);
    assert!(
        !matches.reactant_matches.pattern_matches[0][0]
            .atom_indices
            .is_empty()
    );

    // Product has one molecule matching the pattern
    assert_eq!(matches.product_matches.pattern_matches[0].len(), 1);
    assert!(
        !matches.product_matches.pattern_matches[0][0]
            .atom_indices
            .is_empty()
    );
}

// ===== Phase 3: RDKit Compatibility Tooling =====

#[test]
fn test_rdkit_compat_check_legacy_format() {
    use crate::query::rdkit_compat::*;

    let result = check_rdkit_compatibility("[C:1][C:2]>>[C:2][C:1]").unwrap();
    assert!(result.is_compatible);
    assert_eq!(result.rdkit_format, RDKitFormat::Legacy);
    assert!(result.warnings.is_empty());
}

#[test]
fn test_rdkit_compat_check_agents_format() {
    use crate::query::rdkit_compat::*;

    let result = check_rdkit_compatibility("[C:1][C:2]>[Pd]>[C:2][C:1]").unwrap();
    assert!(result.is_compatible);
    assert_eq!(result.rdkit_format, RDKitFormat::WithAgents);
    assert!(result.warnings.is_empty());
}

#[test]
fn test_rdkit_compat_wildcard_warning() {
    use crate::query::rdkit_compat::*;

    let result = check_rdkit_compatibility("[*:1]>>[*:1]").unwrap();
    assert!(!result.is_compatible);
    assert!(result.warnings.iter().any(|w| w.contains("wildcard")));
}

#[test]
fn test_rdkit_compat_recursive_smarts_warning() {
    use crate::query::rdkit_compat::*;

    // Use simpler recursive SMARTS pattern
    let result = check_rdkit_compatibility("[#6;$([#6]~[#8])]>>[#6]").unwrap();
    assert!(!result.is_compatible);
    assert!(result.warnings.iter().any(|w| w.contains("recursive")));
}

#[test]
fn test_rdkit_compat_missing_delimiter() {
    use crate::query::rdkit_compat::*;

    let result = check_rdkit_compatibility("[C:1][C:2][C:1][C:2]");
    assert!(result.is_err());
}

#[test]
fn test_rdkit_format_detection_legacy() {
    use crate::query::rdkit_compat::*;

    let result = check_rdkit_compatibility("C>>C").unwrap();
    assert_eq!(result.rdkit_format, RDKitFormat::Legacy);
}

#[test]
fn test_rdkit_format_detection_agents() {
    use crate::query::rdkit_compat::*;

    let result = check_rdkit_compatibility("C>[Pd]>C").unwrap();
    assert_eq!(result.rdkit_format, RDKitFormat::WithAgents);
}

#[test]
fn test_rdkit_match_summary_complete() {
    use crate::query::rdkit_compat::*;

    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#6]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    let summary = match_to_rdkit_format(&matches);
    assert!(summary.matched);
    assert_eq!(summary.reactant_count, 1);
    assert_eq!(summary.product_count, 1);
    assert_eq!(summary.reactant_molecules_matched, 1);
    assert_eq!(summary.product_molecules_matched, 1);
}

#[test]
fn test_rdkit_match_summary_incomplete() {
    use crate::query::rdkit_compat::*;

    let rxn = rxn("CC>>CC");
    let query = parse_reaction_query("[#7]>>[#6]").unwrap();
    let matches = get_reaction_smarts_matches(&rxn, &query);

    let summary = match_to_rdkit_format(&matches);
    assert!(!summary.matched);
    assert_eq!(summary.reactant_molecules_matched, 0);
}

#[test]
fn test_rdkit_validate_conventions_valid() {
    use crate::query::rdkit_compat::*;

    let result = validate_rdkit_conventions("[C:1][C:2]>>[C:2][C:1]");
    assert!(result.is_ok());
}

#[test]
fn test_rdkit_validate_conventions_no_delimiter() {
    use crate::query::rdkit_compat::*;

    let result = validate_rdkit_conventions("[C:1][C:2][C:1][C:2]");
    assert!(result.is_err());
}

#[test]
fn test_rdkit_validate_conventions_triple_arrow() {
    use crate::query::rdkit_compat::*;

    let result = validate_rdkit_conventions("[C:1]>>>[C:1]");
    assert!(result.is_err());
}

#[test]
fn test_rdkit_validate_conventions_map_number_zero() {
    use crate::query::rdkit_compat::*;

    let result = validate_rdkit_conventions("[C:0]>>[C:0]");
    assert!(result.is_err());
    if let Err(msg) = result {
        assert!(msg.contains("outside RDKit range"));
    }
}

#[test]
fn test_rdkit_validate_conventions_map_number_out_of_range() {
    use crate::query::rdkit_compat::*;

    let result = validate_rdkit_conventions("[C:1000]>>[C:1000]");
    assert!(result.is_err());
    if let Err(msg) = result {
        assert!(msg.contains("outside RDKit range"));
    }
}

#[test]
fn test_rdkit_validate_conventions_map_number_valid() {
    use crate::query::rdkit_compat::*;

    let result = validate_rdkit_conventions("[C:1][N:999]>>[N:999][C:1]");
    assert!(result.is_ok());
}

#[test]
fn test_rdkit_compat_agents_pipes() {
    use crate::query::rdkit_compat::*;

    let result = check_rdkit_compatibility("[C:1]>[Pd]|[Ni]>[C:1]").unwrap();
    assert!(result.is_compatible);
    assert_eq!(result.rdkit_format, RDKitFormat::WithAgents);
}

#[test]
fn test_rdkit_compat_multiple_patterns() {
    use crate::query::rdkit_compat::*;

    let result = check_rdkit_compatibility("[#6]|[#7]>>[#8]|[#9]").unwrap();
    assert!(result.is_compatible);
}

// B7 Enhancement Tests: Batch querying and pattern libraries

#[test]
fn test_batch_query_reactions_empty() {
    let reactions: Vec<Reaction> = vec![];
    let smarts = "[C:1]>>[C:1]";

    let result = batch_query_reactions(&reactions, smarts);
    assert!(result.is_ok());

    let batch = result.unwrap();
    assert_eq!(batch.total_reactions, 0);
    assert_eq!(batch.matching_reactions, 0);
    assert_eq!(batch.match_percentage, 0.0);
}

#[test]
fn test_batch_query_limits_reject_oversized_batch() {
    let reactions = vec![parse_reaction("C>>C").unwrap()];
    let err = batch_query_reactions_with_limits(
        &reactions,
        "[C:1]>>[C:1]",
        &BatchQueryLimits {
            max_reactions: 0,
            ..Default::default()
        },
    );
    assert!(matches!(
        err,
        Err(ReactionQueryError::ResourceLimit {
            resource: "batch reactions",
            ..
        })
    ));
}

#[test]
fn test_batch_query_reactions_single_match() {
    let rxn = parse_reaction("C>>C").unwrap();
    let reactions = vec![rxn];
    let smarts = "[C:1]>>[C:1]";

    let result = batch_query_reactions(&reactions, smarts);
    assert!(result.is_ok());

    let batch = result.unwrap();
    assert_eq!(batch.total_reactions, 1);
    assert!(batch.matching_reactions > 0);
}

#[test]
fn test_batch_query_reactions_multiple() {
    let rxn1 = parse_reaction("C>>C").unwrap();
    let rxn2 = parse_reaction("CC>>C").unwrap();
    let rxn3 = parse_reaction("CCC>>CC").unwrap();

    let reactions = vec![rxn1, rxn2, rxn3];
    let smarts = "[C:1]>>[C:1]";

    let result = batch_query_reactions(&reactions, smarts);
    assert!(result.is_ok());

    let batch = result.unwrap();
    assert_eq!(batch.total_reactions, 3);
    assert_eq!(batch.matches.len(), 3);
}

#[test]
fn test_batch_query_results_indices() {
    let rxn1 = parse_reaction("C>>C").unwrap();
    let rxn2 = parse_reaction("CC>>C").unwrap();
    let rxn3 = parse_reaction("CCC>>CC").unwrap();

    let reactions = vec![rxn1, rxn2, rxn3];
    let smarts = "[C:1]>>[C:1]";

    let batch = batch_query_reactions(&reactions, smarts).unwrap();

    // Both matching and non-matching indices should return valid results
    let matching = batch.matching_indices();
    let non_matching = batch.non_matching_indices();

    assert_eq!(matching.len() + non_matching.len(), batch.total_reactions);
}

#[test]
fn test_reaction_pattern_library_new() {
    let library = ReactionPatternLibrary::new();
    assert!(library.is_empty());
    assert_eq!(library.len(), 0);
}

#[test]
fn test_reaction_pattern_library_add_pattern() {
    let mut library = ReactionPatternLibrary::new();

    let pattern = parse_reaction_smarts("[C:1]>>[C:1]").unwrap();
    library.add_pattern("simple_carbon".to_string(), pattern);

    assert_eq!(library.len(), 1);
    assert!(!library.is_empty());
    assert!(library.get("simple_carbon").is_some());
}

#[test]
fn test_reaction_pattern_library_add_from_smarts() {
    let mut library = ReactionPatternLibrary::new();

    let result = library.add_pattern_from_smarts(
        "acylation".to_string(),
        "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[N:3]",
    );
    assert!(result.is_ok());
    assert_eq!(library.len(), 1);
}

#[test]
fn test_reaction_pattern_library_multiple_patterns() {
    let mut library = ReactionPatternLibrary::new();

    let _ = library.add_pattern_from_smarts("reaction1".to_string(), "[C:1]>>[C:1]");
    let _ = library.add_pattern_from_smarts("reaction2".to_string(), "[N:1]>>[N:1]");
    let _ = library.add_pattern_from_smarts("reaction3".to_string(), "[O:1]>>[O:1]");

    assert_eq!(library.len(), 3);

    let names = library.pattern_names();
    assert_eq!(names.len(), 3);
    assert!(names.contains(&"reaction1".to_string()));
}

#[test]
fn test_batch_query_with_library() {
    let mut library = ReactionPatternLibrary::new();
    let _ = library.add_pattern_from_smarts("carbon_only".to_string(), "[C:1]>>[C:1]");

    let rxn1 = parse_reaction("C>>C").unwrap();
    let rxn2 = parse_reaction("CC>>C").unwrap();
    let reactions = vec![rxn1, rxn2];

    let results = batch_query_with_library(&reactions, &library);

    assert!(!results.is_empty());
    assert!(results.contains_key("carbon_only"));

    let batch = &results["carbon_only"];
    assert_eq!(batch.total_reactions, 2);
}

#[test]
fn test_query_reaction_simple() {
    let rxn = parse_reaction("C>>C").unwrap();
    let smarts = "[C:1]>>[C:1]";

    let result = query_reaction(&rxn, smarts);
    assert!(result.is_ok());

    let match_result = result.unwrap();
    // Check that it produced a valid result structure
    assert!(!match_result.reactant_matches.pattern_matches.is_empty());
}

#[test]
fn test_batch_query_match_percentage() {
    let rxn1 = parse_reaction("C>>C").unwrap();
    let rxn2 = parse_reaction("CC>>C").unwrap();
    let rxn3 = parse_reaction("CCC>>CC").unwrap();

    let reactions = vec![rxn1, rxn2, rxn3];
    let smarts = "[C:1]>>[C:1]";

    let batch = batch_query_reactions(&reactions, smarts).unwrap();

    // Match percentage should be between 0 and 100
    assert!(batch.match_percentage >= 0.0);
    assert!(batch.match_percentage <= 100.0);

    // Verify calculation
    let expected_percentage =
        (batch.matching_reactions as f64 / batch.total_reactions as f64) * 100.0;
    assert!((batch.match_percentage - expected_percentage).abs() < 0.01);
}

#[test]
fn test_reaction_pattern_library_default() {
    let library = ReactionPatternLibrary::default();
    assert!(library.is_empty());
}
