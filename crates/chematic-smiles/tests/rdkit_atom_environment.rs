//! Original MIT OR Apache-2.0 fixtures; pinned oracle manifest accompanies gate.
use chematic_core::{Atom, Element, MoleculeBuilder};
use chematic_smiles::{
    RdkitAtomEnvironmentLimits, parse, rdkit_atom_environment, rdkit_atom_environment_with_limits,
};
use std::collections::BTreeMap;

#[test]
fn radius_zero_and_linear_exhaustion_have_distinct_maps() {
    let m = parse("CCO").unwrap();
    let r = rdkit_atom_environment(&m, 0, 1, false, true).unwrap();
    assert!(r.bond_indices.is_empty());
    assert_eq!(r.atom_distances, BTreeMap::from([(1, 0)]));
    let r = rdkit_atom_environment(&m, 2, 0, false, true).unwrap();
    assert_eq!(r.bond_indices, [0, 1]);
    assert_eq!(r.atom_distances, BTreeMap::from([(0, 0), (1, 1), (2, 2)]));
    let r = rdkit_atom_environment(&m, 3, 0, false, true).unwrap();
    assert!(r.bond_indices.is_empty() && r.atom_distances.is_empty());
    let r = rdkit_atom_environment(&m, u32::MAX, 0, false, false).unwrap();
    assert_eq!(r.bond_indices, [0, 1]);
    assert_eq!(r.atom_distances.len(), 3);
}

#[test]
fn pending_ring_duplicates_preserve_upstream_enforce_size_behavior() {
    let m = parse("C1CCC1").unwrap();
    for radius in [2, 3] {
        let r = rdkit_atom_environment(&m, radius, 0, false, true).unwrap();
        assert_eq!(r.bond_indices, [0, 3, 1, 2]);
        assert_eq!(
            r.atom_distances,
            BTreeMap::from([(0, 0), (1, 1), (2, 2), (3, 1)])
        );
    }
    assert!(
        rdkit_atom_environment(&m, 4, 0, false, true)
            .unwrap()
            .atom_distances
            .is_empty()
    );
    assert_eq!(
        rdkit_atom_environment(&m, 2, 3, false, true)
            .unwrap()
            .bond_indices,
        [2, 3, 1, 0]
    );
}

#[test]
fn hydrogen_filter_checks_neighbor_and_not_root_element() {
    let m = parse("[2H]C").unwrap();
    assert_eq!(
        rdkit_atom_environment(&m, 1, 0, false, true)
            .unwrap()
            .bond_indices,
        [0]
    );
    assert!(
        rdkit_atom_environment(&m, 1, 1, false, true)
            .unwrap()
            .atom_distances
            .is_empty()
    );
    assert_eq!(
        rdkit_atom_environment(&m, 1, 1, true, true)
            .unwrap()
            .bond_indices,
        [0]
    );
    let m = parse("[H][H]").unwrap();
    assert!(
        rdkit_atom_environment(&m, 1, 0, false, true)
            .unwrap()
            .bond_indices
            .is_empty()
    );
}

#[test]
fn dative_and_disconnected_components_follow_graph_semantics() {
    let m = parse("N->[Cu].O").unwrap();
    assert_eq!(
        rdkit_atom_environment(&m, 1, 1, false, true)
            .unwrap()
            .bond_indices,
        [0]
    );
    assert!(
        rdkit_atom_environment(&m, 1, 2, false, true)
            .unwrap()
            .atom_distances
            .is_empty()
    );
    assert_eq!(
        rdkit_atom_environment(&m, 1, 2, false, false)
            .unwrap()
            .atom_distances,
        BTreeMap::from([(2, 0)])
    );
}

#[test]
fn frontier_and_visit_budgets_fail_without_partial_results() {
    let m = parse("CC(C)C").unwrap();
    for limits in [
        RdkitAtomEnvironmentLimits {
            max_visits: 2,
            max_frontier: 100,
        },
        RdkitAtomEnvironmentLimits {
            max_visits: 100,
            max_frontier: 2,
        },
    ] {
        assert!(
            rdkit_atom_environment_with_limits(&m, 1, 1, false, false, &limits)
                .unwrap_err()
                .to_string()
                .contains("resource limit")
        );
    }
    let limits = RdkitAtomEnvironmentLimits {
        max_visits: 0,
        max_frontier: 0,
    };
    assert_eq!(
        rdkit_atom_environment_with_limits(&m, 0, 1, false, true, &limits)
            .unwrap()
            .atom_distances,
        BTreeMap::from([(1, 0)])
    );
    // Radius 3 in this ring visits a duplicate frontier entry after four bonds.
    let m = parse("C1CCC1").unwrap();
    let limits = RdkitAtomEnvironmentLimits {
        max_visits: 4,
        max_frontier: 100,
    };
    assert!(rdkit_atom_environment_with_limits(&m, 3, 0, false, true, &limits).is_err());
}

#[test]
fn limits_precede_preparation_and_h_removal() {
    let m = parse(&vec!["[H]C"; 513].join(".")).unwrap();
    assert!(
        rdkit_atom_environment(&m, 0, 0, false, true)
            .unwrap_err()
            .to_string()
            .contains("1026 atoms")
    );
    let mut b = MoleculeBuilder::new();
    for _ in 0..65_536 {
        b.add_atom(Atom::new(Element::C));
    }
    assert!(
        rdkit_atom_environment(&b.build(), 0, 0, false, true)
            .unwrap_err()
            .to_string()
            .contains("resource limit")
    );
}

#[test]
fn exact_atom_cap_and_huge_radius_are_bounded_by_graph_exhaustion() {
    let m = parse(&"C".repeat(1024)).unwrap();
    let r = rdkit_atom_environment(&m, u32::MAX, 0, false, false).unwrap();
    assert_eq!(r.bond_indices, (0..1023).collect::<Vec<_>>());
    assert_eq!(r.atom_distances[&1023], 1023);
    assert!(
        rdkit_atom_environment(&m, u32::MAX, 0, false, true)
            .unwrap()
            .atom_distances
            .is_empty()
    );
}

#[test]
fn invalid_roots_query_bonds_and_empty_molecules_fail_explicitly() {
    let m = parse("[H]C([H])([H])[H]").unwrap();
    assert!(rdkit_atom_environment(&m, 0, 1, false, true).is_err());
    assert!(rdkit_atom_environment(&m, 0, usize::MAX, false, true).is_err());
    assert!(rdkit_atom_environment(&parse("C~N").unwrap(), 1, 0, true, true).is_err());
    assert!(rdkit_atom_environment(&MoleculeBuilder::new().build(), 0, 0, false, true).is_err());
}
