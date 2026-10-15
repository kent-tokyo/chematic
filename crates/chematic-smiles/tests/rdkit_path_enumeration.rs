//! Original synthetic MIT OR Apache-2.0 paths, pinned RDKit 2026.09.1.
use chematic_smiles::{
    RdkitPathLimits, RdkitPathParams, parse, rdkit_paths_of_length,
    rdkit_paths_of_length_with_limits,
};

#[test]
fn ordered_atom_expansion_and_bond_set_deduplication() {
    let mol = parse("C1CC1").unwrap();
    let p = RdkitPathParams::default();
    assert_eq!(
        rdkit_paths_of_length(&mol, 1, &p).unwrap(),
        [vec![0], vec![2], vec![1]]
    );
    assert_eq!(
        rdkit_paths_of_length(&mol, 2, &p).unwrap(),
        [vec![0, 1], vec![2, 1], vec![0, 2]]
    );
    let p = RdkitPathParams {
        use_bonds: false,
        ..p
    };
    assert_eq!(
        rdkit_paths_of_length(&mol, 3, &p).unwrap(),
        [vec![0, 1, 2], vec![0, 2, 1], vec![1, 0, 2]]
    );
}

#[test]
fn final_step_closes_ring_at_an_interior_atom() {
    let mol = parse("C1CCC1C").unwrap();
    let p = RdkitPathParams {
        rooted_at_atom: Some(4),
        ..Default::default()
    };
    assert_eq!(
        rdkit_paths_of_length(&mol, 5, &p).unwrap(),
        [vec![3, 4, 0, 1, 2]]
    );
    let p = RdkitPathParams {
        use_bonds: false,
        ..p
    };
    assert_eq!(
        rdkit_paths_of_length(&mol, 6, &p).unwrap(),
        [vec![4, 3, 0, 1, 2, 3]]
    );
    let mol = parse("CC(C)C").unwrap();
    assert!(
        rdkit_paths_of_length(&mol, 3, &Default::default())
            .unwrap()
            .is_empty()
    );
}

#[test]
fn shortest_paths_exclude_ring_detours() {
    let mol = parse("C1CCC1C").unwrap();
    let p = RdkitPathParams {
        only_shortest_paths: true,
        ..Default::default()
    };
    assert_eq!(
        rdkit_paths_of_length(&mol, 3, &p).unwrap(),
        [vec![0, 4, 3], vec![1, 2, 3]]
    );
    assert!(
        rdkit_paths_of_length(&parse("C1CC1").unwrap(), 3, &p)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn hydrogen_singletons_and_shortest_mode_follow_upstream() {
    let mol = parse("[2H]C").unwrap();
    let p = RdkitPathParams {
        use_bonds: false,
        ..Default::default()
    };
    assert_eq!(
        rdkit_paths_of_length(&mol, 1, &p).unwrap(),
        [vec![0], vec![1]]
    );
    assert!(rdkit_paths_of_length(&mol, 2, &p).unwrap().is_empty());
    for p in [
        RdkitPathParams { use_hs: true, ..p },
        RdkitPathParams {
            only_shortest_paths: true,
            ..p
        },
    ] {
        assert_eq!(rdkit_paths_of_length(&mol, 2, &p).unwrap(), [vec![0, 1]]);
    }
}

#[test]
fn empty_zero_impossible_and_outside_roots() {
    let mol = parse("C.C").unwrap();
    let zero = RdkitPathLimits {
        max_visits: 0,
        max_frontier: 0,
        max_results: 0,
        max_path_indices: 0,
    };
    for use_bonds in [false, true] {
        let p = RdkitPathParams {
            use_bonds,
            ..Default::default()
        };
        for length in [0, u32::MAX] {
            assert!(
                rdkit_paths_of_length_with_limits(&mol, length, &p, &zero)
                    .unwrap()
                    .is_empty()
            );
        }
        let p = RdkitPathParams {
            rooted_at_atom: Some(99),
            ..p
        };
        assert!(
            rdkit_paths_of_length_with_limits(&mol, 1, &p, &zero)
                .unwrap()
                .is_empty()
        );
    }
    assert!(
        rdkit_paths_of_length(
            &chematic_core::MoleculeBuilder::new().build(),
            1,
            &Default::default()
        )
        .unwrap()
        .is_empty()
    );
}

#[test]
fn budgets_are_exact_and_do_not_return_truncated_data() {
    let mol = parse("CCO").unwrap();
    for limits in [
        RdkitPathLimits {
            max_visits: 0,
            ..Default::default()
        },
        RdkitPathLimits {
            max_frontier: 2,
            ..Default::default()
        },
        RdkitPathLimits {
            max_results: 1,
            ..Default::default()
        },
        RdkitPathLimits {
            max_path_indices: 2,
            ..Default::default()
        },
    ] {
        assert!(
            rdkit_paths_of_length_with_limits(&mol, 1, &Default::default(), &limits)
                .unwrap_err()
                .to_string()
                .contains("resource limit")
        );
    }
    let p = RdkitPathParams {
        rooted_at_atom: Some(0),
        ..Default::default()
    };
    let limits = RdkitPathLimits {
        max_visits: 1,
        max_frontier: 1,
        max_results: 1,
        max_path_indices: 2,
    };
    assert_eq!(
        rdkit_paths_of_length_with_limits(&mol, 1, &p, &limits).unwrap(),
        [vec![0]]
    );
    let limits = RdkitPathLimits {
        max_path_indices: 1,
        ..limits
    };
    assert!(rdkit_paths_of_length_with_limits(&mol, 1, &p, &limits).is_err());
}

#[test]
fn limits_precede_hydrogen_suppression_and_source_is_immutable() {
    let text = std::iter::repeat_n("[H]C", 513)
        .collect::<Vec<_>>()
        .join(".");
    assert!(
        rdkit_paths_of_length(&parse(&text).unwrap(), 0, &Default::default())
            .unwrap_err()
            .to_string()
            .contains("1026 atoms")
    );
    let mol = parse("N->[Cu+2]").unwrap();
    let before = chematic_smiles::canonical_smiles(&mol);
    assert_eq!(
        rdkit_paths_of_length(&mol, 1, &Default::default()).unwrap(),
        [vec![0]]
    );
    assert_eq!(chematic_smiles::canonical_smiles(&mol), before);
}
