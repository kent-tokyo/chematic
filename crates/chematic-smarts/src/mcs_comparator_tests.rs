use super::*;
use chematic_smiles::parse;

#[test]
fn chiral_mcs_uses_mapped_neighbor_parity_not_smiles_tag_spelling() {
    let config = McsConfig {
        match_chiral_tag: true,
        ..Default::default()
    };
    let record: serde_json::Value = serde_json::from_str(include_str!(
        "../../../benchmarks/2026-10-22-mcs-stereo-neighbor-parity.json"
    ))
    .unwrap();
    let pairs = record["pairs"].as_array().unwrap();
    assert_eq!(pairs.len(), 7);
    for pair in pairs {
        let left = parse(pair["left"].as_str().unwrap()).unwrap();
        let right = parse(pair["right"].as_str().unwrap()).unwrap();
        for refs in [[&left, &right], [&right, &left]] {
            let outcome = find_mcs_with_config_checked(&refs, &config);
            assert!(matches!(outcome, McsOutcome::Exhaustive(_)));
            let query = outcome.into_query();
            assert_eq!(
                query.atom_count(),
                pair["expected"][0].as_u64().unwrap() as usize
            );
            assert_eq!(
                query.bonds.len(),
                pair["expected"][1].as_u64().unwrap() as usize
            );
        }
    }
}

#[test]
fn chiral_mcs_checks_each_molecule_in_a_multi_input_mapping() {
    let first = parse("F[C@H](Cl)Br").unwrap();
    let same_ordered_differently = parse("Cl[C@@H](F)Br").unwrap();
    let same_again = parse("Br[C@@H](Cl)F").unwrap();
    let config = McsConfig {
        match_chiral_tag: true,
        ..Default::default()
    };
    let outcome =
        find_mcs_with_config_checked(&[&first, &same_ordered_differently, &same_again], &config);
    assert!(matches!(outcome, McsOutcome::Exhaustive(_)));
    let query = outcome.into_query();
    assert_eq!([query.atom_count(), query.bonds.len()], [4, 3]);
}

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
    assert_eq!(find_mcs_with_config(&[&a, &a], &config).atom_count(), 0);
    let config = McsConfig {
        atom_compare: AtomCompare::Any,
        ..Default::default()
    };
    assert_eq!(find_mcs_with_config(&[&a, &b], &config).atom_count(), 2);
}

#[test]
fn rdkit_heavy_mode_matches_explicit_hydrogen_only_to_hydrogen() {
    let config = McsConfig {
        atom_compare: AtomCompare::AnyHeavyAtom,
        ..Default::default()
    };
    for (left_smiles, right_smiles, atoms, bonds) in [
        ("[2H]O[2H]", "[H]OC", 2, 1),
        ("[H][H]", "[H][H]", 2, 1),
        ("[H][H]", "NN", 0, 0),
        ("*C", "NC", 1, 0),
    ] {
        let left = parse(left_smiles).unwrap();
        let right = parse(right_smiles).unwrap();
        for refs in [[&left, &right], [&right, &left]] {
            let outcome = find_mcs_with_rdkit_any_heavy_checked(&refs, &config);
            assert!(matches!(outcome, McsOutcome::Exhaustive(_)));
            let query = outcome.into_query();
            assert_eq!(
                query.atom_count(),
                atoms,
                "{left_smiles} versus {right_smiles}"
            );
            assert_eq!(
                query.bonds.len(),
                bonds,
                "{left_smiles} versus {right_smiles}"
            );
        }
    }
}

#[test]
fn rdkit_heavy_mode_respects_explicit_isotope_matching() {
    let left = parse("[2H]O").unwrap();
    let right = parse("[H]N").unwrap();
    let config = McsConfig {
        atom_compare: AtomCompare::AnyHeavyAtom,
        match_isotope: true,
        ..Default::default()
    };
    let query = find_mcs_with_rdkit_any_heavy_checked(&[&left, &right], &config).into_query();
    assert_eq!(query.atom_count(), 1);
    assert_eq!(query.bonds.len(), 0);
}

#[test]
fn rdkit_heavy_mode_matches_pinned_mcs_count_oracle() {
    let rows =
        include_str!("../../../validation/results/mcs-comparator-bounds-2026-10-17/oracle.jsonl");
    let mut comparisons = 0;
    for line in rows.lines() {
        let row: serde_json::Value = serde_json::from_str(line).unwrap();
        if row["atom_mode"] != "any_heavy_atom" {
            continue;
        }
        let inputs: Vec<_> = row["inputs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| parse(value.as_str().unwrap()).unwrap())
            .collect();
        let refs: Vec<_> = inputs.iter().collect();
        let config = McsConfig {
            atom_compare: AtomCompare::AnyHeavyAtom,
            bond_compare: if row["bond_mode"] == "any" {
                BondCompare::Any
            } else {
                BondCompare::OrderOrAromatic
            },
            ..Default::default()
        };
        let outcome = find_mcs_with_rdkit_any_heavy_checked(&refs, &config);
        assert!(matches!(outcome, McsOutcome::Exhaustive(_)), "{line}");
        let query = outcome.into_query();
        assert_eq!(
            query.atom_count(),
            row["expected"][0].as_u64().unwrap() as usize,
            "{line}"
        );
        assert_eq!(
            query.bonds.len(),
            row["expected"][1].as_u64().unwrap() as usize,
            "{line}"
        );
        comparisons += 1;
    }
    assert_eq!(comparisons, 660);
}

#[test]
fn pinned_mcs_profile_counts_remain_exact() {
    let rows =
        include_str!("../../../validation/results/mcs-comparator-bounds-2026-10-17/oracle.jsonl");
    let mut comparisons = 0;
    for line in rows.lines() {
        let row: serde_json::Value = serde_json::from_str(line).unwrap();
        let inputs: Vec<_> = row["inputs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| parse(value.as_str().unwrap()).unwrap())
            .collect();
        let refs: Vec<_> = inputs.iter().collect();
        let config = McsConfig {
            atom_compare: match row["atom_mode"].as_str().unwrap() {
                "elements" => AtomCompare::Elements,
                "any" => AtomCompare::Any,
                "any_heavy_atom" => AtomCompare::AnyHeavyAtom,
                other => panic!("unexpected atom mode {other}"),
            },
            bond_compare: match row["bond_mode"].as_str().unwrap() {
                "order_or_aromatic" => BondCompare::OrderOrAromatic,
                "any" => BondCompare::Any,
                other => panic!("unexpected bond mode {other}"),
            },
            ..Default::default()
        };
        let outcome = find_mcs_with_config_checked(&refs, &config);
        assert!(matches!(outcome, McsOutcome::Exhaustive(_)), "{line}");
        let query = outcome.into_query();
        assert_eq!(
            [query.atom_count(), query.bonds.len()],
            [
                row["native_profile_expected"][0].as_u64().unwrap() as usize,
                row["native_profile_expected"][1].as_u64().unwrap() as usize,
            ],
            "{line}"
        );
        comparisons += 1;
    }
    assert_eq!(comparisons, 1980);
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

#[test]
fn complete_rings_only_keeps_smaller_valid_chain_after_larger_partial_ring() {
    let six_ring_and_chain = parse("C1CCCCC1.CCCC").unwrap();
    let five_ring_and_chain = parse("C1CCCC1.CCCC").unwrap();
    let config = McsConfig {
        complete_rings_only: true,
        ..Default::default()
    };
    for refs in [
        [&six_ring_and_chain, &five_ring_and_chain],
        [&five_ring_and_chain, &six_ring_and_chain],
    ] {
        let outcome = find_mcs_with_config_checked(&refs, &config);
        assert!(matches!(outcome, McsOutcome::Exhaustive(_)));
        let query = outcome.into_query();
        assert_eq!(query.atom_count(), 4);
        assert_eq!(query.bonds.len(), 3);
        assert!(!crate::find_matches(&query, refs[0]).is_empty());
        assert!(!crate::find_matches(&query, refs[1]).is_empty());
    }
}
