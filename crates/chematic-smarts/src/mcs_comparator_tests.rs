use super::*;
use chematic_smiles::parse;

#[test]
fn comparator_cross_element_chains_keep_the_full_common_topology() {
    for (a, b, size) in [("CCC", "NNN", 3), ("CCCO", "NCCN", 4), ("COC", "NNN", 3)] {
        let a = parse(a).unwrap();
        let b = parse(b).unwrap();
        for atom_compare in [AtomCompare::Any, AtomCompare::AnyHeavyAtom] {
            for refs in [[&a, &b], [&b, &a]] {
                let config = McsConfig {
                    atom_compare: atom_compare.clone(),
                    ..Default::default()
                };
                let outcome = find_mcs_with_config_checked(&refs, &config);
                assert!(matches!(outcome, McsOutcome::Exhaustive(_)));
                assert_eq!(outcome.into_query().atom_count(), size, "{atom_compare:?}");
            }
        }
    }
}

#[test]
fn comparator_cross_element_three_molecule_search_is_not_underbounded() {
    let mols = [
        parse("CCCC").unwrap(),
        parse("NNN").unwrap(),
        parse("OOO").unwrap(),
    ];
    for atom_compare in [AtomCompare::Any, AtomCompare::AnyHeavyAtom] {
        let config = McsConfig {
            atom_compare,
            ..Default::default()
        };
        assert_eq!(
            find_mcs_with_config(&mols.iter().collect::<Vec<_>>(), &config).atom_count(),
            3
        );
    }
}

#[test]
fn comparator_elements_remain_element_specific() {
    let a = parse("CCC").unwrap();
    let b = parse("NNN").unwrap();
    assert_eq!(find_mcs(&[&a, &b]).atom_count(), 0);
}

#[test]
fn comparator_heavy_mode_does_not_match_explicit_hydrogen() {
    let a = parse("[H][H]").unwrap();
    let b = parse("NN").unwrap();
    let config = McsConfig {
        atom_compare: AtomCompare::AnyHeavyAtom,
        ..Default::default()
    };
    assert_eq!(find_mcs_with_config(&[&a, &b], &config).atom_count(), 0);
    let config = McsConfig {
        atom_compare: AtomCompare::Any,
        ..Default::default()
    };
    assert_eq!(find_mcs_with_config(&[&a, &b], &config).atom_count(), 2);
}

#[test]
fn comparator_search_budget_still_reports_resource_limit() {
    let a = parse("CCCC").unwrap();
    let b = parse("NNNN").unwrap();
    let config = McsConfig {
        atom_compare: AtomCompare::Any,
        max_search_nodes: Some(1),
        ..Default::default()
    };
    assert!(matches!(
        find_mcs_with_config_checked(&[&a, &b], &config),
        McsOutcome::ResourceLimited(_)
    ));
}
