//! Original MIT OR Apache-2.0 meso/full-alternative graph proof contracts.
use chematic_chem::{MesoStereoGroupError as Error, reduce_meso_stereo_group as reduce};
use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, CipCode, DerivedSlot, Element, Molecule, MoleculeBuilder,
    StereoGroup, StereoGroupKind as Kind,
};
use chematic_smiles::parse;
fn group(kind: Kind, atoms: &[u32], bonds: &[u32]) -> (StereoGroup, Vec<BondIdx>) {
    (
        StereoGroup::new(kind, atoms.iter().copied().map(AtomIdx).collect()),
        bonds.iter().copied().map(BondIdx).collect(),
    )
}
fn grouped(s: &str, a: &[u32], b: &[u32]) -> Molecule {
    parse(s)
        .unwrap()
        .with_checked_stereo_groups(vec![group(Kind::And(1), a, b)])
        .unwrap()
}
fn meso() -> Molecule {
    grouped("C[C@H](O)[C@H](O)C", &[1, 3], &[])
}
#[test]
fn atom_only_meso_group_becomes_absolute_with_a_witness() {
    let m = meso();
    let r = diagnostic(reduce(&m, 0, 1000000)).unwrap();
    assert!(r.reduced);
    assert_eq!(r.checked_pairs, 1);
    assert_eq!(r.witnesses.len(), 1);
    assert_eq!(r.molecule.stereo_groups()[0].kind, Kind::Absolute);
    assert_eq!(
        r.witnesses[0].atom_permutation,
        vec![
            AtomIdx(5),
            AtomIdx(3),
            AtomIdx(4),
            AtomIdx(1),
            AtomIdx(2),
            AtomIdx(0)
        ]
    );
}
#[test]
fn mixed_atom_and_bond_group_preserves_its_whole_correlation() {
    let m = grouped("C/C=C/[C@H](C)/C=C\\C", &[3], &[1, 5]);
    let r = diagnostic(reduce(&m, 0, 1000000)).unwrap();
    assert!(r.reduced);
    assert_eq!(
        r.molecule.stereo_group_bond_indices(0).unwrap(),
        [BondIdx(1), BondIdx(5)]
    );
    assert_eq!(r.molecule.stereo_groups()[0].atom_indices, [AtomIdx(3)]);
    assert_eq!(r.molecule.stereo_groups()[0].kind, Kind::Absolute);
}
#[test]
fn true_racemate_and_center_only_alternative_are_retained() {
    for m in [
        grouped("C[C@H](O)[C@@H](O)C", &[1, 3], &[]),
        grouped("C/C=C/[C@H](C)/C=C\\C", &[3], &[]),
    ] {
        let r = diagnostic(reduce(&m, 0, 1000000)).unwrap();
        assert!(!r.reduced);
        assert_eq!(r.molecule.stereo_groups(), m.stereo_groups());
        assert!(r.witnesses.is_empty());
    }
}
#[test]
fn bond_only_group_reduces_only_when_geometry_maps() {
    for (s, reduced) in [("C/C=C/C=C\\C", true), ("C/C=C/C=C/C", false)] {
        let m = grouped(s, &[], &[1, 3]);
        assert_eq!(diagnostic(reduce(&m, 0, 1000000)).unwrap().reduced, reduced);
    }
}
#[test]
fn or_groups_use_the_same_full_alternative_proof() {
    let m = parse("C[C@H](O)[C@H](O)C")
        .unwrap()
        .with_checked_stereo_groups(vec![group(Kind::Or(7), &[1, 3], &[])])
        .unwrap();
    assert!(diagnostic(reduce(&m, 0, 1000000)).unwrap().reduced);
}
#[test]
fn existing_absolute_members_are_merged_and_other_groups_remain_relative() {
    let m = parse("C[C@H](O)[C@H](O)C.N[C@H](F)Cl.[C@H](F)(Cl)Br")
        .unwrap()
        .with_checked_stereo_groups(vec![
            group(Kind::Absolute, &[7], &[]),
            group(Kind::And(8), &[1, 3], &[]),
            group(Kind::Or(9), &[10], &[]),
        ])
        .unwrap();
    let r = diagnostic(reduce(&m, 1, 1000000)).unwrap();
    assert!(r.reduced);
    assert_eq!(r.checked_pairs, 2);
    assert_eq!(
        r.witnesses.iter().map(|w| w.mask).collect::<Vec<_>>(),
        [0, 2]
    );
    assert_eq!(r.molecule.stereo_groups().len(), 2);
    assert_eq!(
        r.molecule.stereo_groups()[0].atom_indices,
        [AtomIdx(7), AtomIdx(1), AtomIdx(3)]
    );
    assert_eq!(r.molecule.stereo_groups()[1].kind, Kind::Or(9));
}
#[test]
fn every_other_group_context_is_proved_with_one_shared_budget() {
    let s = format!("C[C@H](O)[C@H](O)C.{}", ["[C@H](F)(Cl)Br"; 5].join("."));
    let mut groups = vec![group(Kind::And(1), &[1, 3], &[])];
    groups.extend((0..5).map(|i| group(Kind::Or(i + 2), &[6 + i * 4], &[])));
    let m = parse(&s)
        .unwrap()
        .with_checked_stereo_groups(groups)
        .unwrap();
    let r = diagnostic(reduce(&m, 0, 1000000)).unwrap();
    assert!(r.reduced);
    assert_eq!(r.checked_pairs, 32);
    assert_eq!(r.witnesses.len(), 32);
    assert!(matches!(
        diagnostic(reduce(&m, 0, r.search_visits - 1)),
        Err(Error::BudgetExceeded { .. })
    ));
    assert_eq!(m.stereo_groups().len(), 6);
}
#[test]
fn isotope_map_and_tag_distinctions_are_conservatively_preserved() {
    for s in [
        "[13CH3][C@H](O)[C@H](O)C",
        "[CH3:11][C@H](O)[C@H](O)[CH3:12]",
    ] {
        assert!(
            !diagnostic(reduce(&grouped(s, &[1, 3], &[]), 0, 1000000))
                .unwrap()
                .reduced
        );
    }
    let mut m = meso();
    m.set_tag(AtomIdx(0), Some(31));
    let r = diagnostic(reduce(&m, 0, 1000000)).unwrap();
    assert!(!r.reduced);
    assert_eq!(r.molecule.atom_tag(AtomIdx(0)), m.atom_tag(AtomIdx(0)));
}
#[test]
fn drawing_graph_frames_and_source_metadata_are_unchanged() {
    let mut m = meso();
    m.set_cip_code(AtomIdx(1), Some(CipCode::S));
    m.derived(DerivedSlot::Sssr, || 77u32);
    let r = diagnostic(reduce(&m, 0, 1000000)).unwrap();
    assert!(r.reduced);
    for (i, a) in m.atoms() {
        assert_eq!(a.chirality, r.molecule.atom(i).chirality);
        assert_eq!(
            m.stereo_neighbor_order(i),
            r.molecule.stereo_neighbor_order(i)
        );
    }
    for (i, b) in m.bonds() {
        assert_eq!(b, r.molecule.bond(i));
    }
    assert_eq!(r.molecule.atom(AtomIdx(1)).cip_code, None);
    assert_eq!(m.atom(AtomIdx(1)).cip_code, Some(CipCode::S));
    assert!(
        r.molecule
            .derived_if_computed::<u32>(DerivedSlot::Sssr)
            .is_none()
    );
    assert_eq!(
        *m.derived_if_computed::<u32>(DerivedSlot::Sssr).unwrap(),
        77
    );
}
#[test]
fn tiny_budget_is_unknown_and_not_a_negative_or_partial_reduction() {
    assert_eq!(
        diagnostic(reduce(&meso(), 0, 1)).err(),
        Some(Error::BudgetExceeded {
            visits: 1,
            limit: 1
        })
    );
    for limit in [0, 10000001, u64::MAX] {
        assert_eq!(
            diagnostic(reduce(&meso(), 0, limit)).err(),
            Some(Error::InvalidVisitLimit)
        );
    }
}
#[test]
fn absolute_and_invalid_indices_are_refused() {
    assert!(matches!(
        diagnostic(reduce(&meso(), 1, 1000000)),
        Err(Error::GroupIndex { .. })
    ));
    let m = parse("C[C@H](O)[C@H](O)C")
        .unwrap()
        .with_checked_stereo_groups(vec![group(Kind::Absolute, &[1, 3], &[])])
        .unwrap();
    assert_eq!(
        diagnostic(reduce(&m, 0, 1000000)).err(),
        Some(Error::AbsoluteGroup)
    );
}
#[test]
fn malformed_stored_references_fail_before_graph_indexing() {
    let mut m = meso();
    m.set_stereo_groups(vec![StereoGroup::new(Kind::And(1), vec![AtomIdx(999)])]);
    assert!(matches!(
        diagnostic(reduce(&m, 0, 1000000)),
        Err(Error::Alternatives(_))
    ));
}
#[test]
fn unassigned_tetrahedral_sites_and_unsupported_global_stereo_fail_closed() {
    let m = grouped("C[C@H](O)[C@H](O)C.FC(Cl)(Br)I", &[1, 3], &[]);
    assert_eq!(
        diagnostic(reduce(&m, 0, 1000000)).err(),
        Some(Error::UnspecifiedTetrahedralAtom(AtomIdx(7)))
    );
    let m = grouped("C[C@H](O)[C@H](O)C.Cl[Pt@SP1](Cl)(N)N", &[1, 3], &[]);
    assert_eq!(
        diagnostic(reduce(&m, 0, 1000000)).err(),
        Some(Error::UnsupportedAtom(AtomIdx(7)))
    );
}
#[test]
fn unsupported_bond_semantics_and_direction_stashes_are_refused() {
    let m = grouped("C[C@H](O)[C@H](O)C.C~C", &[1, 3], &[]);
    assert!(matches!(
        diagnostic(reduce(&m, 0, 1000000)),
        Err(Error::UnsupportedBond(_))
    ));
    let mut m = meso();
    m.set_bond_direction(BondIdx(0), BondOrder::Up);
    assert_eq!(
        diagnostic(reduce(&m, 0, 1000000)).err(),
        Some(Error::UnsupportedDirections)
    );
}
#[test]
fn unsupported_graph_sizes_fail_before_enumeration() {
    let mut b = MoleculeBuilder::new();
    for _ in 0..65 {
        b.add_atom(Atom::new(Element::C));
    }
    assert_eq!(
        diagnostic(reduce(&b.build(), 0, 1000000)).err(),
        Some(Error::GraphLimit)
    );
    let mut b = MoleculeBuilder::new();
    let atoms: Vec<_> = (0..17).map(|_| b.add_atom(Atom::new(Element::C))).collect();
    for i in 0..17 {
        for j in i + 1..17 {
            b.add_bond(atoms[i], atoms[j], BondOrder::Single).unwrap();
        }
    }
    assert_eq!(
        diagnostic(reduce(&b.build(), 0, 1000000)).err(),
        Some(Error::GraphLimit)
    );
}
#[test]
fn incomplete_unselected_absolute_frame_is_not_ignored() {
    let mut m = parse("C[C@H](O)[C@H](O)C.N[C@H](F)Cl")
        .unwrap()
        .with_checked_stereo_groups(vec![
            group(Kind::And(1), &[1, 3], &[]),
            group(Kind::Absolute, &[7], &[]),
        ])
        .unwrap();
    m.set_stereo_neighbor_order(AtomIdx(7), vec![]);
    assert_eq!(
        diagnostic(reduce(&m, 0, 1000000)).err(),
        Some(Error::IncompleteFrame(AtomIdx(7)))
    );
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
