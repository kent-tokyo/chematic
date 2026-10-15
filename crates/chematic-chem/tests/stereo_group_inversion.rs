use chematic_chem::{StereoGroupInversionError as Error, invert_stereo_group};
use chematic_core::{
    AtomIdx, BondIdx, BondOrder, CipCode, DerivedSlot, Molecule, StereoGroup, StereoGroupKind,
};
use chematic_smiles::{canonical_smiles, parse};
fn grouped(smiles: &str, atoms: &[u32], bonds: &[u32]) -> Molecule {
    parse(smiles)
        .unwrap()
        .with_checked_stereo_groups(vec![(
            StereoGroup::new(
                StereoGroupKind::And(1),
                atoms.iter().map(|&a| AtomIdx(a)).collect(),
            ),
            bonds.iter().map(|&b| BondIdx(b)).collect(),
        )])
        .unwrap()
}
fn key(s: &str) -> String {
    canonical_smiles(&parse(s).unwrap())
}
#[test]
fn mixed_group_inverts_center_and_double_bond_together() {
    let m = grouped("C/C=C/[C@H](O)CC", &[3], &[1]);
    let original = canonical_smiles(&m);
    let out = diagnostic(invert_stereo_group(&m, 0)).unwrap();
    assert_eq!(canonical_smiles(&out), key("C/C=C\\[C@@H](O)CC"));
    assert_eq!(out.stereo_groups(), m.stereo_groups());
    assert_eq!(
        out.stereo_group_bond_indices(0),
        m.stereo_group_bond_indices(0)
    );
    assert_eq!(canonical_smiles(&m), original);
}
#[test]
fn shared_marker_solver_preserves_unselected_conjugated_bond_geometry() {
    for (bonds, expected) in [
        (vec![1], "C/C=C\\C=C\\C"),
        (vec![3], "C/C=C/C=C\\C"),
        (vec![1, 3], "C/C=C\\C=C/C"),
    ] {
        let m = grouped("C/C=C/C=C/C", &[], &bonds);
        assert_eq!(
            canonical_smiles(&diagnostic(invert_stereo_group(&m, 0)).unwrap()),
            key(expected)
        );
    }
}
#[test]
fn two_applications_restore_configuration_for_atom_bond_and_mixed_groups() {
    for (s, a, b) in [
        ("N[C@@H](C)C(=O)O", vec![1], vec![]),
        ("C/C=C/C=C/C", vec![], vec![1, 3]),
        ("C/C=C/[C@H](O)CC", vec![3], vec![1]),
    ] {
        let m = grouped(s, &a, &b);
        let twice = diagnostic(invert_stereo_group(
            &diagnostic(invert_stereo_group(&m, 0)).unwrap(),
            0,
        ))
        .unwrap();
        assert_eq!(canonical_smiles(&twice), canonical_smiles(&m));
        assert_eq!(twice.stereo_groups(), m.stereo_groups());
    }
}
#[test]
fn unselected_center_and_group_are_preserved() {
    let m = parse("C/C=C/[C@H](O)CC.F[C@](Cl)(Br)I").unwrap();
    let m = m
        .with_checked_stereo_groups(vec![
            (
                StereoGroup::new(StereoGroupKind::And(1), vec![AtomIdx(3)]),
                vec![BondIdx(1)],
            ),
            (
                StereoGroup::new(StereoGroupKind::Absolute, vec![AtomIdx(8)]),
                vec![],
            ),
        ])
        .unwrap();
    let out = diagnostic(invert_stereo_group(&m, 0)).unwrap();
    assert_eq!(out.atom(AtomIdx(8)).chirality, m.atom(AtomIdx(8)).chirality);
    assert_eq!(
        out.stereo_neighbor_order(AtomIdx(8)),
        m.stereo_neighbor_order(AtomIdx(8))
    );
    assert_eq!(out.stereo_groups(), m.stereo_groups());
}
#[test]
fn metadata_and_reference_frames_survive_while_cip_and_caches_are_cold() {
    let mut m = grouped("C/C=C/[C@H](O)CC", &[3], &[1]);
    m.set_tag(AtomIdx(0), Some(23));
    m.set_cip_code(AtomIdx(3), Some(CipCode::R));
    m.derived(DerivedSlot::Sssr, || 11u32);
    let out = diagnostic(invert_stereo_group(&m, 0)).unwrap();
    assert_eq!(out.atom_tag(AtomIdx(0)), m.atom_tag(AtomIdx(0)));
    for (idx, _) in m.atoms() {
        assert_eq!(out.stereo_neighbor_order(idx), m.stereo_neighbor_order(idx));
    }
    assert_eq!(out.atom(AtomIdx(3)).cip_code, None);
    assert_eq!(m.atom(AtomIdx(3)).cip_code, Some(CipCode::R));
    assert!(out.derived_if_computed::<u32>(DerivedSlot::Sssr).is_none());
    assert_eq!(
        *m.derived_if_computed::<u32>(DerivedSlot::Sssr).unwrap(),
        11
    );
}
#[test]
fn absolute_and_out_of_range_groups_have_no_relative_alternative() {
    let mut m = grouped("F[C@](Cl)(Br)I", &[1], &[]);
    assert!(matches!(
        diagnostic(invert_stereo_group(&m, 1)),
        Err(Error::GroupIndex { .. })
    ));
    m.set_stereo_groups(vec![StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomIdx(1)],
    )]);
    assert_eq!(
        diagnostic(invert_stereo_group(&m, 0)).err(),
        Some(Error::AbsoluteGroup)
    );
}
#[test]
fn unspecified_or_non_tetrahedral_atoms_are_refused() {
    for (s, idx) in [("FC(Cl)(Br)I", 1), ("Cl[Pt@SP1](Cl)(N)N", 1)] {
        let m = grouped(s, &[idx], &[]);
        assert_eq!(
            diagnostic(invert_stereo_group(&m, 0)).err(),
            Some(Error::UnsupportedAtom(AtomIdx(idx)))
        );
    }
}
#[test]
fn incomplete_or_mismatched_tetrahedral_frames_are_refused() {
    for frame in [vec![], vec![0, 2, 3, 99], vec![0, 2, 2, 4]] {
        let mut m = grouped("F[C@](Cl)(Br)I", &[1], &[]);
        m.set_stereo_neighbor_order(AtomIdx(1), frame);
        assert_eq!(
            diagnostic(invert_stereo_group(&m, 0)).err(),
            Some(Error::IncompleteFrame(AtomIdx(1)))
        );
    }
}
#[test]
fn unsupported_unspecified_and_small_ring_bonds_are_refused() {
    for (s, b, expected) in [
        ("CC", 0, Error::UnsupportedBond(BondIdx(0))),
        ("CC=CC", 1, Error::UnspecifiedBond(BondIdx(1))),
        ("C1/C=C/CCCC1", 1, Error::SmallRingBond(BondIdx(1))),
    ] {
        assert_eq!(
            diagnostic(invert_stereo_group(&grouped(s, &[], &[b]), 0)).err(),
            Some(expected)
        );
    }
}
#[test]
fn invalid_stored_group_references_fail_before_graph_access() {
    let mut m = grouped("F/C=C/Cl", &[], &[1]);
    m.set_stereo_groups_with_bonds(vec![(
        StereoGroup::new(StereoGroupKind::And(1), vec![]),
        vec![BondIdx(99)],
    )]);
    assert!(matches!(
        diagnostic(invert_stereo_group(&m, 0)),
        Err(Error::InvalidGroups(_))
    ));
}
#[test]
fn unsupported_stored_directions_fail_atomically() {
    let mut m = grouped("F/C=C/Cl", &[], &[1]);
    m.set_bond_direction(BondIdx(0), BondOrder::Up);
    let before = m.bonds().map(|(_, b)| b.clone()).collect::<Vec<_>>();
    assert_eq!(
        diagnostic(invert_stereo_group(&m, 0)).err(),
        Some(Error::UnrepresentableDirections)
    );
    assert_eq!(
        m.bonds().map(|(_, b)| b.clone()).collect::<Vec<_>>(),
        before
    );
}
#[test]
fn or_groups_use_the_same_full_alternative_and_explicit_h_frames_work() {
    let m = parse("[H][C@](F)(Cl)Br").unwrap();
    let m = m
        .with_checked_stereo_groups(vec![(
            StereoGroup::new(StereoGroupKind::Or(4), vec![AtomIdx(1)]),
            vec![],
        )])
        .unwrap();
    let out = diagnostic(invert_stereo_group(&m, 0)).unwrap();
    assert_ne!(out.atom(AtomIdx(1)).chirality, m.atom(AtomIdx(1)).chirality);
    assert_eq!(out.stereo_groups()[0].kind, StereoGroupKind::Or(4));
}

#[test]
fn eight_membered_ring_is_not_rejected_as_a_small_ring() {
    let m = grouped("C1/C=C/CCCCC1", &[], &[1]);
    let out = diagnostic(invert_stereo_group(&m, 0)).unwrap();
    assert_eq!(
        canonical_smiles(&diagnostic(invert_stereo_group(&out, 0)).unwrap()),
        canonical_smiles(&m)
    );
}
#[test]
fn conflicting_markers_are_rejected_without_repairing_input() {
    let mut m = grouped("F/C(Cl)=C(/Br)I", &[], &[2]);
    m.set_bond_order(BondIdx(1), BondOrder::Down);
    assert_eq!(
        diagnostic(invert_stereo_group(&m, 0)).err(),
        Some(Error::ConflictingDirections)
    );
    assert_eq!(m.bond(BondIdx(1)).order, BondOrder::Down);
}

#[test]
fn legacy_duplicate_atom_membership_is_preserved_and_flipped_once() {
    let mut m = parse("F[C@](Cl)(Br)I").unwrap();
    m.set_stereo_groups(vec![StereoGroup {
        kind: StereoGroupKind::And(1),
        atom_indices: vec![AtomIdx(1), AtomIdx(1)],
    }]);
    let out = diagnostic(invert_stereo_group(&m, 0)).unwrap();
    assert_eq!(out.stereo_groups(), m.stereo_groups());
    assert_ne!(out.atom(AtomIdx(1)).chirality, m.atom(AtomIdx(1)).chirality);
    assert_eq!(
        canonical_smiles(&diagnostic(invert_stereo_group(&out, 0)).unwrap()),
        canonical_smiles(&m)
    );
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
