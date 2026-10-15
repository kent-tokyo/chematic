use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Element, Molecule, MoleculeBuilder, RGroupLabel,
    RenumberError, STEREO_H_SENTINEL, StereoGroup, StereoGroupKind,
};
use chematic_smiles::{canonical_smiles, parse};

fn identity(m: &Molecule) -> Vec<AtomIdx> {
    (0..m.atom_count()).map(|i| AtomIdx(i as u32)).collect()
}
fn reverse(m: &Molecule) -> Vec<AtomIdx> {
    identity(m).into_iter().rev().collect()
}
fn assert_indexed_graph(source: &Molecule, output: &Molecule, order: &[AtomIdx]) {
    let mut inverse = vec![AtomIdx(0); order.len()];
    for (new, &old) in order.iter().enumerate() {
        inverse[old.0 as usize] = AtomIdx(new as u32);
        assert_eq!(output.atom(AtomIdx(new as u32)), source.atom(old));
    }
    assert_eq!(output.bond_count(), source.bond_count());
    for (idx, bond) in source.bonds() {
        let b = output.bond(idx);
        assert_eq!(
            (b.atom1, b.atom2, b.order),
            (
                inverse[bond.atom1.0 as usize],
                inverse[bond.atom2.0 as usize],
                bond.order
            )
        );
    }
}

#[test]
fn issue9629_identity_and_reverse_keep_one_stereocenter() {
    let m = parse("C1CCN2CCCC2C1").unwrap();
    assert_eq!(chematic_chem::num_stereocenters(&m), 1);
    let before = canonical_smiles(&m);
    for order in [identity(&m), reverse(&m)] {
        let out = diagnostic(m.renumber_atoms(&order)).unwrap();
        assert_indexed_graph(&m, &out, &order);
        assert_eq!(chematic_chem::num_stereocenters(&out), 1);
        assert_eq!(chematic_chem::num_unspecified_stereocenters(&out), 1);
        assert_eq!(canonical_smiles(&out), before);
    }
}

fn permutations(order: &mut [AtomIdx], start: usize, visit: &mut impl FnMut(&[AtomIdx])) {
    if start == order.len() {
        visit(order);
        return;
    }
    for i in start..order.len() {
        order.swap(start, i);
        permutations(order, start + 1, visit);
        order.swap(start, i);
    }
}

#[test]
fn all_120_tetrahedral_permutations_preserve_configuration() {
    let m = parse("F[C@](Cl)(Br)I").unwrap();
    let before = canonical_smiles(&m);
    let mut count = 0;
    permutations(&mut identity(&m), 0, &mut |order| {
        let out = diagnostic(m.renumber_atoms(order)).unwrap();
        assert_indexed_graph(&m, &out, order);
        assert_eq!(canonical_smiles(&out), before);
        assert_eq!(chematic_chem::num_stereocenters(&out), 1);
        count += 1;
    });
    assert_eq!(count, 120);
}

#[test]
fn explicit_implicit_isotopic_h_and_atom_maps_survive() {
    for smiles in [
        "N[C@@H:7](C)C(=O)O",
        "[H][C@:9](F)(Cl)Br",
        "[2H][C@](F)(Cl)Br",
    ] {
        let m = parse(smiles).unwrap();
        let order = reverse(&m);
        let out = diagnostic(m.renumber_atoms(&order)).unwrap();
        assert_indexed_graph(&m, &out, &order);
        assert_eq!(canonical_smiles(&out), canonical_smiles(&m));
        let inverse: Vec<_> = (0..m.atom_count())
            .map(|i| (m.atom_count() - 1 - i) as u32)
            .collect();
        for (a, _) in m.atoms() {
            if let Some(frame) = m.stereo_neighbor_order(a) {
                let expected: Vec<_> = frame
                    .iter()
                    .map(|&n| {
                        if n == STEREO_H_SENTINEL {
                            n
                        } else {
                            inverse[n as usize]
                        }
                    })
                    .collect();
                assert_eq!(
                    out.stereo_neighbor_order(AtomIdx(inverse[a.0 as usize])),
                    Some(expected.as_slice())
                );
            }
        }
    }
}

#[test]
fn axial_ez_and_coordination_stereo_survive() {
    for s in [
        "F/C=C/Cl",
        "F/C=C\\Cl",
        "CC=[C@AL1]=C(C)F",
        "Cl[Pt@SP1](Cl)(N)N",
        "F[P@TB1](Cl)(Br)(I)N",
        "F[Co@OH1](Cl)(Br)(I)(N)O",
    ] {
        let m = parse(s).unwrap();
        let order = reverse(&m);
        let out = diagnostic(m.renumber_atoms(&order)).unwrap();
        assert_indexed_graph(&m, &out, &order);
        assert_eq!(canonical_smiles(&out), canonical_smiles(&m), "{s}");
    }
}

#[test]
fn cyclic_numberings_keep_bicyclic_descriptor_counts() {
    for s in [
        "C1CCN2CCCC2C1",
        "C1CCN2CCC2C1",
        "C1CN2CCC1CC2",
        "C1CC2CCC1N2",
    ] {
        let m = parse(s).unwrap();
        let n = m.atom_count();
        let expected = chematic_chem::num_stereocenters(&m);
        for start in 0..n {
            let order: Vec<_> = (0..n).map(|i| AtomIdx(((start + i) % n) as u32)).collect();
            let out = diagnostic(m.renumber_atoms(&order)).unwrap();
            assert_eq!(
                chematic_chem::num_stereocenters(&out),
                expected,
                "{s} {start}"
            );
        }
    }
}

#[test]
fn remaps_tags_r_labels_and_enhanced_groups_but_keeps_bond_indices() {
    let mut m = parse("F[C@](Cl)(Br)I.*").unwrap();
    m.set_tag(AtomIdx(1), Some(431));
    m.set_r_group(AtomIdx(5), RGroupLabel::numbered(7).unwrap());
    m.set_stereo_groups_with_bonds(vec![(
        StereoGroup::new(StereoGroupKind::Or(3), vec![AtomIdx(1)]),
        vec![BondIdx(2)],
    )]);
    let out = diagnostic(m.renumber_atoms(&reverse(&m))).unwrap();
    assert_eq!(out.atom_tag(AtomIdx(4)).unwrap().get(), 431);
    assert_eq!(out.r_group_label(AtomIdx(0)).unwrap().number(), Some(7));
    assert_eq!(
        out.stereo_groups()[0],
        StereoGroup::new(StereoGroupKind::Or(3), vec![AtomIdx(4)])
    );
    assert_eq!(out.stereo_group_bond_indices(0), Some(&[BondIdx(2)][..]));
}

#[test]
fn preserves_direction_anchors_ring_closures_and_dative_orientation() {
    let mut m = parse("N->[Cu].C1CC1").unwrap();
    m.set_bond_direction(BondIdx(0), BondOrder::Up);
    m.set_bond_direction_anchor(BondIdx(0), AtomIdx(0));
    let n = m.atom_count();
    let out = diagnostic(m.renumber_atoms(&reverse(&m))).unwrap();
    assert_eq!(out.bond(BondIdx(0)).atom1, AtomIdx((n - 1) as u32));
    assert_eq!(out.bond(BondIdx(0)).atom2, AtomIdx((n - 2) as u32));
    assert_eq!(out.bond(BondIdx(0)).order, BondOrder::Dative);
    assert_eq!(out.bond_direction(BondIdx(0)), Some(BondOrder::Up));
    assert_eq!(
        out.bond_direction_anchor(BondIdx(0)),
        Some(AtomIdx((n - 1) as u32))
    );
    for (b, _) in m.bonds() {
        assert_eq!(
            out.smiles_ring_closure_metadata(b),
            m.smiles_ring_closure_metadata(b)
        );
    }
}

#[test]
fn derived_caches_start_empty_even_after_identity_permutation() {
    use chematic_core::derived_cache::DerivedSlot;
    let m = parse("C1CC1").unwrap();
    m.seed_derived(
        DerivedSlot::CanonicalRanks,
        std::sync::Arc::new(vec![42u32]),
    );
    assert_eq!(m.derived_rdkit_num_rings(|| 7), 7);
    let out = diagnostic(m.renumber_atoms(&identity(&m))).unwrap();
    assert!(
        out.derived_if_computed::<Vec<u32>>(DerivedSlot::CanonicalRanks)
            .is_none()
    );
    assert_eq!(out.derived_rdkit_num_rings(|| 1), 1);
    assert_eq!(m.derived_rdkit_num_rings(|| 99), 7);
    assert_eq!(
        *m.derived_if_computed::<Vec<u32>>(DerivedSlot::CanonicalRanks)
            .unwrap(),
        vec![42]
    );
}

#[test]
fn inverse_permutation_restores_indexed_graph_and_frames() {
    let m = parse("N[C@@H](C)C(=O)O.F/C=C/Cl").unwrap();
    let n = m.atom_count();
    let order: Vec<_> = (0..n).map(|i| AtomIdx(((i + 3) % n) as u32)).collect();
    let mut inv = vec![AtomIdx(0); n];
    for (new, old) in order.iter().enumerate() {
        inv[old.0 as usize] = AtomIdx(new as u32);
    }
    let out = m
        .renumber_atoms(&order)
        .unwrap()
        .renumber_atoms(&inv)
        .unwrap();
    assert_indexed_graph(&m, &out, &identity(&m));
    for (a, _) in m.atoms() {
        assert_eq!(out.stereo_neighbor_order(a), m.stereo_neighbor_order(a));
    }
}

#[test]
fn invalid_orders_are_typed_errors_and_empty_graph_is_valid() {
    let m = parse("CCO").unwrap();
    assert!(matches!(
        diagnostic(m.renumber_atoms(&[AtomIdx(0)])),
        Err(RenumberError::Length {
            expected: 3,
            actual: 1
        })
    ));
    assert_eq!(
        diagnostic(m.renumber_atoms(&[AtomIdx(0), AtomIdx(1), AtomIdx(3)])).err(),
        Some(RenumberError::InvalidAtom(AtomIdx(3)))
    );
    assert_eq!(
        diagnostic(m.renumber_atoms(&[AtomIdx(0), AtomIdx(1), AtomIdx(0)])).err(),
        Some(RenumberError::DuplicateAtom(AtomIdx(0)))
    );
    assert!(MoleculeBuilder::new().build().renumber_atoms(&[]).is_ok());
}

#[test]
fn invalid_caller_metadata_fails_without_panicking() {
    let m = parse("CC").unwrap();
    let mut variants = Vec::new();
    let mut v = m.clone();
    v.set_stereo_neighbor_order(AtomIdx(0), vec![99]);
    variants.push(v);
    let mut v = m.clone();
    v.set_stereo_neighbor_order(AtomIdx(99), vec![0]);
    variants.push(v);
    let mut v = m.clone();
    v.add_stereo_group(StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomIdx(99)],
    ));
    variants.push(v);
    let mut v = m.clone();
    v.set_bond_direction(BondIdx(99), BondOrder::Up);
    variants.push(v);
    let mut v = m.clone();
    v.set_bond_direction_anchor(BondIdx(0), AtomIdx(99));
    variants.push(v);
    let mut v = m.clone();
    v.add_stereo_group_with_bonds(
        StereoGroup::new(StereoGroupKind::Absolute, vec![]),
        vec![BondIdx(99)],
    );
    variants.push(v);
    for source in variants {
        assert!(matches!(
            diagnostic(source.renumber_atoms(&identity(&source))),
            Err(RenumberError::InvalidMetadata { .. })
        ));
    }
}

#[test]
fn builder_stereo_without_recorded_frame_keeps_ordered_adjacency() {
    let mut b = MoleculeBuilder::new();
    let mut atom = Atom::new(Element::C);
    atom.chirality = chematic_core::Chirality::Clockwise;
    b.add_atom(atom);
    for e in [Element::F, Element::CL, Element::BR, Element::I] {
        let n = b.add_atom(Atom::new(e));
        b.add_bond(AtomIdx(0), n, BondOrder::Single).unwrap();
    }
    let m = b.build();
    let out = diagnostic(m.renumber_atoms(&reverse(&m))).unwrap();
    assert_eq!(
        out.neighbors(AtomIdx(4))
            .map(|(a, _)| a.0)
            .collect::<Vec<_>>(),
        vec![3, 2, 1, 0]
    );
    assert_eq!(canonical_smiles(&out), canonical_smiles(&m));
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
