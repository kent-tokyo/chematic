use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Chirality, DerivedSlot, Element, Molecule, MoleculeBuilder,
    RGroupLabel, STEREO_H_SENTINEL, StereoGroup, StereoGroupError, StereoGroupKind,
};

fn molecule() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    for _ in 0..5 {
        builder.add_atom(Atom::new(Element::C));
    }
    for i in 0..4 {
        builder
            .add_bond(AtomIdx(i), AtomIdx(i + 1), BondOrder::Single)
            .unwrap();
    }
    builder.set_stereo_neighbor_order(AtomIdx(2), vec![1, 3, STEREO_H_SENTINEL]);
    builder.set_bond_direction(BondIdx(1), BondOrder::Up);
    builder.set_bond_direction_anchor(BondIdx(1), AtomIdx(2));
    builder.set_r_group(AtomIdx(4), RGroupLabel::numbered(3).unwrap());
    builder.set_tag(AtomIdx(0), Some(37));
    let mut m = builder.build();
    m.set_chirality(AtomIdx(2), Chirality::Clockwise);
    m
}
fn group(kind: StereoGroupKind, atoms: &[u32], bonds: &[u32]) -> (StereoGroup, Vec<BondIdx>) {
    (
        StereoGroup::new(kind, atoms.iter().map(|&i| AtomIdx(i)).collect()),
        bonds.iter().map(|&i| BondIdx(i)).collect(),
    )
}
fn error(groups: Vec<(StereoGroup, Vec<BondIdx>)>) -> StereoGroupError {
    let error = molecule()
        .with_checked_stereo_groups(groups)
        .err()
        .expect("expected structural error");
    assert!(error.to_string().contains("group"));
    error
}
#[test]
fn mixed_bond_only_and_atom_only_groups_preserve_order_and_kind() {
    let m = molecule();
    let groups = vec![
        group(StereoGroupKind::And(7), &[2], &[1]),
        group(StereoGroupKind::Absolute, &[0], &[]),
        group(StereoGroupKind::Or(7), &[], &[3]),
    ];
    let result = m.with_checked_stereo_groups(groups.clone()).unwrap();
    for (i, (g, b)) in groups.iter().enumerate() {
        assert_eq!(&result.stereo_groups()[i], g);
        assert_eq!(result.stereo_group_bond_indices(i), Some(b.as_slice()));
    }
    assert!(m.stereo_groups().is_empty());
}
#[test]
fn checked_copy_retains_metadata_and_has_cold_caches() {
    let m = molecule();
    m.derived(DerivedSlot::Sssr, || 19u32);
    let result = m
        .with_checked_stereo_groups(vec![group(StereoGroupKind::And(1), &[2], &[1])])
        .unwrap();
    assert_eq!(
        result.atoms().collect::<Vec<_>>(),
        m.atoms().collect::<Vec<_>>()
    );
    assert_eq!(
        result.bonds().collect::<Vec<_>>(),
        m.bonds().collect::<Vec<_>>()
    );
    for (i, _) in m.atoms() {
        assert_eq!(
            result.neighbors(i).collect::<Vec<_>>(),
            m.neighbors(i).collect::<Vec<_>>()
        );
        assert_eq!(result.stereo_neighbor_order(i), m.stereo_neighbor_order(i));
        assert_eq!(result.atom_tag(i), m.atom_tag(i));
        assert_eq!(result.r_group_label(i), m.r_group_label(i));
    }
    assert_eq!(
        result.bond_direction(BondIdx(1)),
        m.bond_direction(BondIdx(1))
    );
    assert_eq!(
        result.bond_direction_anchor(BondIdx(1)),
        m.bond_direction_anchor(BondIdx(1))
    );
    assert!(
        result
            .derived_if_computed::<u32>(DerivedSlot::Sssr)
            .is_none()
    );
    assert_eq!(
        *m.derived_if_computed::<u32>(DerivedSlot::Sssr).unwrap(),
        19
    );
}
#[test]
fn duplicates_within_one_group_deduplicate_in_first_seen_order() {
    let raw = StereoGroup {
        kind: StereoGroupKind::And(2),
        atom_indices: vec![AtomIdx(3), AtomIdx(1), AtomIdx(3)],
    };
    let result = molecule()
        .with_checked_stereo_groups(vec![(raw, vec![BondIdx(3), BondIdx(1), BondIdx(3)])])
        .unwrap();
    assert_eq!(
        result.stereo_groups()[0].atom_indices,
        vec![AtomIdx(3), AtomIdx(1)]
    );
    assert_eq!(
        result.stereo_group_bond_indices(0),
        Some([BondIdx(3), BondIdx(1)].as_slice())
    );
}
#[test]
fn rejects_empty_group_but_empty_list_clears_metadata() {
    assert_eq!(
        error(vec![group(StereoGroupKind::Absolute, &[], &[])]),
        StereoGroupError::EmptyGroup { group: 0 }
    );
    let m = molecule()
        .with_checked_stereo_groups(vec![group(StereoGroupKind::And(1), &[2], &[1])])
        .unwrap();
    let cleared = m.with_checked_stereo_groups(vec![]).unwrap();
    assert!(cleared.stereo_groups().is_empty());
    assert_eq!(m.stereo_groups().len(), 1);
}
#[test]
fn rejects_zero_relative_group_numbers() {
    for kind in [StereoGroupKind::And(0), StereoGroupKind::Or(0)] {
        assert_eq!(
            error(vec![group(kind, &[0], &[])]),
            StereoGroupError::InvalidGroupNumber { group: 0 }
        );
    }
}
#[test]
fn rejects_duplicate_keys_even_when_members_are_disjoint() {
    for kind in [
        StereoGroupKind::Absolute,
        StereoGroupKind::And(2),
        StereoGroupKind::Or(2),
    ] {
        assert_eq!(
            error(vec![group(kind.clone(), &[0], &[]), group(kind, &[1], &[])]),
            StereoGroupError::DuplicateGroup {
                first: 0,
                second: 1
            }
        );
    }
}
#[test]
fn invalid_indices_fail_without_accessing_graph_memory() {
    for index in [5, u32::MAX] {
        assert_eq!(
            error(vec![group(StereoGroupKind::Absolute, &[index], &[])]),
            StereoGroupError::InvalidAtom {
                group: 0,
                atom: AtomIdx(index)
            }
        );
    }
    for index in [4, u32::MAX] {
        assert_eq!(
            error(vec![group(StereoGroupKind::Absolute, &[], &[index])]),
            StereoGroupError::InvalidBond {
                group: 0,
                bond: BondIdx(index)
            }
        );
    }
}
#[test]
fn shared_atom_or_bond_cannot_have_conflicting_groups() {
    assert_eq!(
        error(vec![
            group(StereoGroupKind::And(1), &[1], &[]),
            group(StereoGroupKind::Or(1), &[1], &[])
        ]),
        StereoGroupError::AtomOverlap {
            first: 0,
            second: 1,
            atom: AtomIdx(1)
        }
    );
    assert_eq!(
        error(vec![
            group(StereoGroupKind::And(1), &[], &[1]),
            group(StereoGroupKind::Absolute, &[], &[1])
        ]),
        StereoGroupError::BondOverlap {
            first: 0,
            second: 1,
            bond: BondIdx(1)
        }
    );
}
#[test]
fn failed_replacement_leaves_original_groups_untouched() {
    let m = molecule()
        .with_checked_stereo_groups(vec![group(StereoGroupKind::And(1), &[1], &[1])])
        .unwrap();
    assert!(
        m.with_checked_stereo_groups(vec![group(StereoGroupKind::Or(1), &[99], &[])])
            .is_err()
    );
    assert_eq!(m.stereo_groups()[0].kind, StereoGroupKind::And(1));
    assert_eq!(
        m.stereo_group_bond_indices(0),
        Some([BondIdx(1)].as_slice())
    );
}
