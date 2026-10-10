use crate::{
    Atom, AtomIdx, BondIdx, BondOrder, Chirality, Element, Molecule, MoleculeBuilder, RGroupLabel,
    STEREO_H_SENTINEL, StereoGroup, StereoGroupKind,
};

fn annotated_chain() -> Molecule {
    let mut b = MoleculeBuilder::new();
    for _ in 0..6 {
        b.add_atom(Atom::new(Element::C));
    }
    for i in 0..5 {
        b.add_bond(AtomIdx(i), AtomIdx(i + 1), BondOrder::Single)
            .unwrap();
    }
    b.set_stereo_neighbor_order(AtomIdx(2), vec![1, 3, STEREO_H_SENTINEL]);
    b.set_bond_direction(BondIdx(3), BondOrder::Up);
    b.set_bond_direction_anchor(BondIdx(3), AtomIdx(4));
    b.set_r_group(AtomIdx(5), RGroupLabel::numbered(17).unwrap());
    b.set_tag(AtomIdx(4), Some(29));
    b.build()
}

#[test]
fn atom_removal_remaps_stereo_directions_tags_and_rgroups_consistently() {
    let original = annotated_chain();
    for removed in 0..6 {
        let (copy, map) = original.with_atom_removed(AtomIdx(removed));
        let mut changed = original.clone();
        assert_eq!(changed.remove_atom(AtomIdx(removed)), map);
        assert_eq!(copy.atom_count(), 5);
        assert_eq!(
            copy.bond_count(),
            if removed == 0 || removed == 5 { 4 } else { 3 }
        );
        for (old, &mapped) in map.iter().enumerate() {
            if let Some(new) = mapped {
                for result in [&copy, &changed] {
                    assert_eq!(result.atom(new), original.atom(AtomIdx(old as u32)));
                    assert_eq!(result.atom_tag(new), original.atom_tag(AtomIdx(old as u32)));
                    assert_eq!(
                        result.r_group_label(new),
                        original.r_group_label(AtomIdx(old as u32))
                    );
                }
            } else {
                assert_eq!(old, removed as usize);
            }
        }
        if let Some(center) = map[2] {
            let expected: Vec<u32> = [1usize, 3]
                .iter()
                .filter_map(|i| map[*i].map(|a| a.0))
                .chain([STEREO_H_SENTINEL])
                .collect();
            assert_eq!(
                copy.stereo_neighbor_order(center),
                Some(expected.as_slice())
            );
            assert_eq!(
                changed.stereo_neighbor_order(center),
                Some(expected.as_slice())
            );
        }
        if let (Some(a), Some(c)) = (map[3], map[4]) {
            for result in [&copy, &changed] {
                let (bond, _) = result.bond_between(a, c).unwrap();
                assert_eq!(result.bond_direction(bond), Some(BondOrder::Up));
                assert_eq!(result.bond_direction_anchor(bond), Some(c));
            }
        }
    }
    assert_eq!(original.atom_count(), 6);
}

#[test]
fn bond_removal_preserves_surviving_direction_anchors_and_atom_metadata() {
    let original = annotated_chain();
    let copy = original.with_bond_removed(BondIdx(0));
    assert_eq!(copy.bond_count(), 4);
    assert_eq!(copy.atom_tag(AtomIdx(4)).unwrap().get(), 29);
    assert_eq!(copy.r_group_label(AtomIdx(5)), RGroupLabel::numbered(17));
    assert_eq!(
        copy.stereo_neighbor_order(AtomIdx(2)),
        original.stereo_neighbor_order(AtomIdx(2))
    );
    assert_eq!(copy.bond_direction(BondIdx(2)), Some(BondOrder::Up));
    assert_eq!(copy.bond_direction_anchor(BondIdx(2)), Some(AtomIdx(4)));
    let fragments = copy.fragments_with_source_atoms();
    assert_eq!(fragments.len(), 2);
    let (chain, source) = fragments.iter().find(|(_, s)| s.len() == 5).unwrap();
    assert_eq!(
        source,
        &vec![AtomIdx(1), AtomIdx(2), AtomIdx(3), AtomIdx(4), AtomIdx(5)]
    );
    assert_eq!(chain.r_group_label(AtomIdx(4)), RGroupLabel::numbered(17));
}

#[test]
fn wildcard_and_dative_edits_keep_atom_identity_and_validate_duplicate_edges() {
    let mut mol = annotated_chain();
    mol.set_chirality(AtomIdx(5), Chirality::Clockwise);
    mol.set_atom_aromatic(AtomIdx(5), true);
    mol.set_wildcard(AtomIdx(5));
    assert!(mol.atom(AtomIdx(5)).wildcard);
    assert!(!mol.atom(AtomIdx(5)).aromatic);
    assert_eq!(mol.atom(AtomIdx(5)).chirality, Chirality::None);
    assert_eq!(mol.r_group_label(AtomIdx(5)), None);
    for donor in [AtomIdx(1), AtomIdx(0)] {
        mol.set_dative_bond(BondIdx(0), donor);
        assert_eq!(mol.bond(BondIdx(0)).atom1, donor);
        assert_eq!(mol.bond(BondIdx(0)).order, BondOrder::Dative);
        assert!(
            mol.neighbors(AtomIdx(0))
                .any(|(a, b)| a == AtomIdx(1) && b == BondIdx(0))
        );
    }
    assert!(
        mol.add_bond(AtomIdx(99), AtomIdx(0), BondOrder::Single)
            .is_err()
    );
    assert!(
        mol.add_bond(AtomIdx(0), AtomIdx(1), BondOrder::Single)
            .is_err()
    );
    let mut builder = MoleculeBuilder::from_molecule(&mol);
    builder.add_stereo_group(StereoGroup::new(
        StereoGroupKind::Or(3),
        vec![AtomIdx(2), AtomIdx(2)],
    ));
    let rebuilt = builder.build();
    assert_eq!(rebuilt.stereo_groups()[0].atom_indices, vec![AtomIdx(2)]);
}
