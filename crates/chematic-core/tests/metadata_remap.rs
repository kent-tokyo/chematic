use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Element, Molecule, MoleculeBuilder, RGroupLabel,
    STEREO_H_SENTINEL, STEREO_VACANT_SENTINEL, StereoGroup, StereoGroupKind,
};

fn four_atom_chain() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    for element in [Element::C, Element::N, Element::O, Element::F] {
        builder.add_atom(Atom::new(element));
    }
    builder
        .add_bond(AtomIdx(0), AtomIdx(1), BondOrder::Single)
        .unwrap();
    builder
        .add_bond(AtomIdx(1), AtomIdx(2), BondOrder::Double)
        .unwrap();
    builder
        .add_bond(AtomIdx(2), AtomIdx(3), BondOrder::Single)
        .unwrap();
    builder.build()
}

#[test]
fn removing_an_atom_remaps_every_indexed_metadata_sidecar() {
    let mut mol = four_atom_chain();
    mol.set_bond_direction(BondIdx(0), BondOrder::Up);
    mol.set_bond_direction(BondIdx(2), BondOrder::Down);
    mol.set_bond_direction_anchor(BondIdx(0), AtomIdx(0));
    mol.set_bond_direction_anchor(BondIdx(2), AtomIdx(3));
    mol.set_bond_stereo_unspecified(BondIdx(0));
    mol.set_bond_stereo_unspecified(BondIdx(2));
    mol.set_stereo_neighbor_order(
        AtomIdx(3),
        vec![AtomIdx(2).0, STEREO_H_SENTINEL, STEREO_VACANT_SENTINEL],
    );
    mol.set_r_group(AtomIdx(3), RGroupLabel::numbered(7).unwrap());
    mol.add_stereo_group_with_bonds(
        StereoGroup::new(StereoGroupKind::And(2), vec![AtomIdx(0), AtomIdx(3)]),
        vec![BondIdx(0), BondIdx(2)],
    );

    let remap = mol.remove_atom(AtomIdx(1));
    assert_eq!(
        remap,
        vec![Some(AtomIdx(0)), None, Some(AtomIdx(1)), Some(AtomIdx(2))]
    );
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(mol.bond_direction(BondIdx(0)), Some(BondOrder::Down));
    assert_eq!(mol.bond_direction_anchor(BondIdx(0)), Some(AtomIdx(2)));
    assert!(mol.bond_stereo_is_unspecified(BondIdx(0)));
    assert_eq!(
        mol.stereo_neighbor_order(AtomIdx(2)),
        Some(&[AtomIdx(1).0, STEREO_H_SENTINEL, STEREO_VACANT_SENTINEL][..])
    );
    assert_eq!(
        mol.r_group_label(AtomIdx(2)).and_then(RGroupLabel::number),
        Some(7)
    );
    assert_eq!(
        mol.stereo_groups()[0].atom_indices,
        vec![AtomIdx(0), AtomIdx(2)]
    );
    assert_eq!(mol.stereo_group_bond_indices(0), Some(&[BondIdx(0)][..]));
}

#[test]
fn removing_a_bond_drops_its_metadata_and_shifts_later_entries() {
    let mut mol = four_atom_chain();
    for (bond, direction, anchor) in [
        (BondIdx(0), BondOrder::Up, AtomIdx(0)),
        (BondIdx(1), BondOrder::Down, AtomIdx(1)),
        (BondIdx(2), BondOrder::Up, AtomIdx(2)),
    ] {
        mol.set_bond_direction(bond, direction);
        mol.set_bond_direction_anchor(bond, anchor);
        mol.set_bond_stereo_unspecified(bond);
    }
    mol.set_stereo_groups_with_bonds(vec![
        (
            StereoGroup::new(StereoGroupKind::Absolute, vec![AtomIdx(0)]),
            vec![BondIdx(0), BondIdx(1), BondIdx(2), BondIdx(2)],
        ),
        (
            StereoGroup::new(StereoGroupKind::Or(4), Vec::new()),
            vec![BondIdx(1)],
        ),
    ]);

    mol.remove_bond(BondIdx(1));
    assert_eq!(mol.bond_count(), 2);
    assert_eq!(mol.bond_direction(BondIdx(0)), Some(BondOrder::Up));
    assert_eq!(mol.bond_direction(BondIdx(1)), Some(BondOrder::Up));
    assert_eq!(mol.bond_direction_anchor(BondIdx(1)), Some(AtomIdx(2)));
    assert!(mol.bond_stereo_is_unspecified(BondIdx(0)));
    assert!(mol.bond_stereo_is_unspecified(BondIdx(1)));
    assert_eq!(
        mol.stereo_group_bond_indices(0),
        Some(&[BondIdx(0), BondIdx(1)][..])
    );
    assert_eq!(mol.stereo_groups().len(), 1);

    mol.clear_bond_direction(BondIdx(0));
    mol.clear_bond_stereo_unspecified(BondIdx(0));
    assert_eq!(mol.bond_direction(BondIdx(0)), None);
    assert!(!mol.bond_stereo_is_unspecified(BondIdx(0)));
}
