use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Element, MoleculeBuilder, StereoGroup, StereoGroupKind,
};

#[test]
fn immutable_removals_remap_unspecified_bonds_without_mutating_source() {
    let mut builder = MoleculeBuilder::new();
    for element in [Element::C, Element::N, Element::O, Element::F] {
        builder.add_atom(Atom::new(element));
    }
    for (a, b, order) in [
        (0, 1, BondOrder::Single),
        (1, 2, BondOrder::Double),
        (2, 3, BondOrder::Single),
    ] {
        builder.add_bond(AtomIdx(a), AtomIdx(b), order).unwrap();
    }
    let mut source = builder.build();
    for i in 0..3 {
        source.set_bond_stereo_unspecified(BondIdx(i));
    }
    source.add_stereo_group_with_bonds(
        StereoGroup::new(StereoGroupKind::Or(1), vec![]),
        vec![BondIdx(0), BondIdx(2)],
    );
    let (without_atom, remap) = source.with_atom_removed(AtomIdx(0));
    assert_eq!(remap[0], None);
    assert_eq!(without_atom.bond_count(), 2);
    assert!(without_atom.bond_stereo_is_unspecified(BondIdx(0)));
    assert!(without_atom.bond_stereo_is_unspecified(BondIdx(1)));
    assert_eq!(
        without_atom.stereo_group_bond_indices(0),
        Some(&[BondIdx(1)][..])
    );
    let without_bond = source.with_bond_removed(BondIdx(0));
    assert_eq!(without_bond.bond_count(), 2);
    assert!(without_bond.bond_stereo_is_unspecified(BondIdx(0)));
    assert!(without_bond.bond_stereo_is_unspecified(BondIdx(1)));
    assert_eq!(
        without_bond.stereo_group_bond_indices(0),
        Some(&[BondIdx(1)][..])
    );
    assert_eq!(source.bond_count(), 3);
    assert_eq!(
        source.stereo_group_bond_indices(0),
        Some(&[BondIdx(0), BondIdx(2)][..])
    );
}
