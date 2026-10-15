//! Original MIT OR Apache-2.0 examples; see validation/rdkit_issues/stereo/group-alternatives.json.
use chematic_chem::{
    StereoGroupAlternativesError as Error, enumerate_stereo_group_alternatives as enumerate,
};
use chematic_core::ez_markers::read_ez_facts;
use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, CipCode, DerivedSlot, Element, Molecule, MoleculeBuilder,
    StereoGroup, StereoGroupKind as Kind,
};
use chematic_smiles::parse;
fn group(kind: Kind, atoms: &[u32], bonds: &[u32]) -> (StereoGroup, Vec<BondIdx>) {
    (
        StereoGroup::new(kind, atoms.iter().map(|&i| AtomIdx(i)).collect()),
        bonds.iter().map(|&i| BondIdx(i)).collect(),
    )
}
fn mixed() -> Molecule {
    parse("C/C=C/[C@H](F)C[C@@H](Cl)/C=C/Br")
        .unwrap()
        .with_checked_stereo_groups(vec![
            group(Kind::And(1), &[3], &[1]),
            group(Kind::Or(2), &[6], &[8]),
        ])
        .unwrap()
}
fn parity(m: &Molecule, bond: u32) -> bool {
    let (facts, conflict) = read_ez_facts(m).unwrap();
    assert!(!conflict);
    let f = facts.iter().find(|f| f.double == BondIdx(bond)).unwrap();
    f.ends[0].up ^ f.ends[1].up
}
#[test]
fn two_mixed_groups_expand_by_full_correlated_mask() {
    let m = mixed();
    let result = diagnostic(enumerate(&m, 64)).unwrap();
    assert_eq!(result.relative_group_indices, [0, 1]);
    assert_eq!(result.molecules.len(), 4);
    for (mask, out) in result.molecules.iter().enumerate() {
        for (bit, atom, bond) in [(0, 3, 1), (1, 6, 8)] {
            let flip = mask & (1 << bit) != 0;
            assert_eq!(
                out.atom(AtomIdx(atom)).chirality != m.atom(AtomIdx(atom)).chirality,
                flip
            );
            assert_eq!(parity(out, bond) != parity(&m, bond), flip);
        }
        assert_eq!(out.stereo_groups(), m.stereo_groups());
    }
}
#[test]
fn conjugated_groups_keep_independent_geometry_with_shared_markers() {
    let m = parse("C/C=C/C=C/C")
        .unwrap()
        .with_checked_stereo_groups(vec![
            group(Kind::And(4), &[], &[1]),
            group(Kind::And(7), &[], &[3]),
        ])
        .unwrap();
    let r = diagnostic(enumerate(&m, 4)).unwrap();
    assert_eq!(r.molecules.len(), 4);
    for (mask, out) in r.molecules.iter().enumerate() {
        for (bit, bond) in [(0, 1), (1, 3)] {
            assert_eq!(parity(out, bond) ^ parity(&m, bond), mask & (1 << bit) != 0);
        }
    }
}
#[test]
fn abs_groups_are_skipped_in_bit_positions_and_remain_fixed() {
    let m = parse("F[C@](Cl)(Br)I.N[C@H](F)[C@@H](Cl)Br")
        .unwrap()
        .with_checked_stereo_groups(vec![
            group(Kind::Absolute, &[1], &[]),
            group(Kind::Or(1), &[6], &[]),
            group(Kind::And(1), &[8], &[]),
        ])
        .unwrap();
    let r = diagnostic(enumerate(&m, 4)).unwrap();
    assert_eq!(r.relative_group_indices, [1, 2]);
    for out in &r.molecules {
        assert_eq!(out.atom(AtomIdx(1)).chirality, m.atom(AtomIdx(1)).chirality);
    }
}
#[test]
fn metadata_graph_frames_and_input_cache_survive() {
    let mut m = mixed();
    m.set_tag(AtomIdx(0), Some(99));
    m.set_cip_code(AtomIdx(3), Some(CipCode::R));
    m.derived(DerivedSlot::Sssr, || 42u32);
    let r = diagnostic(enumerate(&m, 64)).unwrap();
    for out in &r.molecules {
        assert_eq!(out.atom_tag(AtomIdx(0)), m.atom_tag(AtomIdx(0)));
        assert_eq!(out.atom(AtomIdx(3)).cip_code, None);
        assert!(out.derived_if_computed::<u32>(DerivedSlot::Sssr).is_none());
        for (idx, _) in m.atoms() {
            assert_eq!(out.stereo_neighbor_order(idx), m.stereo_neighbor_order(idx));
        }
        for (idx, bond) in m.bonds() {
            assert_eq!(
                (out.bond(idx).atom1, out.bond(idx).atom2),
                (bond.atom1, bond.atom2)
            );
        }
        assert_eq!(out.stereo_groups(), m.stereo_groups());
        for i in 0..m.stereo_groups().len() {
            assert_eq!(
                out.stereo_group_bond_indices(i),
                m.stereo_group_bond_indices(i)
            );
        }
    }
    assert_eq!(m.atom(AtomIdx(3)).cip_code, Some(CipCode::R));
    assert_eq!(
        *m.derived_if_computed::<u32>(DerivedSlot::Sssr).unwrap(),
        42
    );
}
#[test]
fn no_groups_and_abs_only_return_one_independent_representative() {
    for m in [
        parse("CC(F)Cl").unwrap(),
        parse("F[C@](Cl)(Br)I")
            .unwrap()
            .with_checked_stereo_groups(vec![group(Kind::Absolute, &[1], &[])])
            .unwrap(),
    ] {
        let r = diagnostic(enumerate(&m, 1)).unwrap();
        assert!(r.relative_group_indices.is_empty());
        assert_eq!(r.molecules.len(), 1);
        assert_eq!(r.molecules[0].atom_count(), m.atom_count());
        for (i, a) in m.atoms() {
            assert_eq!(r.molecules[0].atom(i).chirality, a.chirality);
        }
    }
}
#[test]
fn all_masks_remain_even_if_symmetry_makes_duplicates() {
    let m = parse("C[C@H](O)[C@H](O)C")
        .unwrap()
        .with_checked_stereo_groups(vec![group(Kind::And(1), &[1, 3], &[])])
        .unwrap();
    assert_eq!(diagnostic(enumerate(&m, 2)).unwrap().molecules.len(), 2);
}
#[test]
fn budget_is_checked_before_generating_and_never_truncates() {
    let m = mixed();
    assert_eq!(
        diagnostic(enumerate(&m, 3)).err(),
        Some(Error::AlternativeLimit {
            required: 4,
            limit: 3
        })
    );
    for limit in [0, 65, usize::MAX] {
        assert_eq!(
            diagnostic(enumerate(&m, limit)).err(),
            Some(Error::InvalidLimit)
        );
    }
    assert_eq!(diagnostic(enumerate(&m, 4)).unwrap().molecules.len(), 4);
}
fn independent_groups(n: usize) -> Molecule {
    let s = std::iter::repeat_n("[C@H](F)(Cl)Br", n)
        .collect::<Vec<_>>()
        .join(".");
    parse(&s)
        .unwrap()
        .with_checked_stereo_groups(
            (0..n)
                .map(|i| group(Kind::And(i as u32 + 1), &[(i * 4) as u32], &[]))
                .collect(),
        )
        .unwrap()
}
#[test]
fn six_groups_produce_64_and_seven_fail() {
    assert_eq!(
        diagnostic(enumerate(&independent_groups(6), 64))
            .unwrap()
            .molecules
            .len(),
        64
    );
    assert_eq!(
        diagnostic(enumerate(&independent_groups(7), 64)).err(),
        Some(Error::TooManyRelativeGroups {
            actual: 7,
            limit: 6
        })
    );
}
#[test]
fn unsupported_later_group_returns_original_index_and_leaves_source_intact() {
    let m = parse("F[C@](Cl)(Br)I.CC=CC")
        .unwrap()
        .with_checked_stereo_groups(vec![
            group(Kind::Absolute, &[1], &[]),
            group(Kind::And(1), &[], &[5]),
        ])
        .unwrap();
    assert!(matches!(
        diagnostic(enumerate(&m, 2)),
        Err(Error::Inversion { group: 1, .. })
    ));
    assert_eq!(m.stereo_group_bond_indices(1).unwrap(), [BondIdx(5)]);
}
#[test]
fn invalid_stored_groups_are_not_indexed_unsafely() {
    let mut m = mixed();
    m.set_stereo_groups(vec![StereoGroup::new(Kind::And(1), vec![AtomIdx(999)])]);
    assert!(matches!(
        diagnostic(enumerate(&m, 64)),
        Err(Error::InvalidGroups(_))
    ));
}
#[test]
fn empty_graph_is_one_empty_assignment() {
    let r = diagnostic(enumerate(&MoleculeBuilder::new().build(), 1)).unwrap();
    assert_eq!(r.molecules.len(), 1);
    assert_eq!(r.molecules[0].atom_count(), 0);
}
#[test]
fn oversized_graphs_fail_before_copying_or_search() {
    let mut builder = MoleculeBuilder::new();
    for _ in 0..1025 {
        builder.add_atom(Atom::new(Element::C));
    }
    assert_eq!(
        diagnostic(enumerate(&builder.build(), 64)).err(),
        Some(Error::AtomLimit {
            actual: 1025,
            limit: 1024
        })
    );
    let mut builder = MoleculeBuilder::new();
    let atoms: Vec<_> = (0..92)
        .map(|_| builder.add_atom(Atom::new(Element::C)))
        .collect();
    for i in 0..92 {
        for j in i + 1..92 {
            builder
                .add_bond(atoms[i], atoms[j], BondOrder::Single)
                .unwrap();
        }
    }
    assert_eq!(
        diagnostic(enumerate(&builder.build(), 64)).err(),
        Some(Error::BondLimit {
            actual: 4186,
            limit: 4096
        })
    );
}
#[test]
fn legacy_duplicate_members_are_retained_but_flipped_once() {
    let mut m = mixed();
    m.set_stereo_groups_with_bonds(vec![group(Kind::And(1), &[3, 3], &[1])]);
    let r = diagnostic(enumerate(&m, 2)).unwrap();
    assert_eq!(r.molecules[1].stereo_groups(), m.stereo_groups());
    assert_ne!(
        r.molecules[1].atom(AtomIdx(3)).chirality,
        m.atom(AtomIdx(3)).chirality
    );
}
#[test]
fn malformed_relative_frames_fail_even_if_other_centers_are_valid() {
    let mut m = mixed();
    m.set_stereo_neighbor_order(AtomIdx(6), vec![]);
    assert!(matches!(
        diagnostic(enumerate(&m, 4)),
        Err(Error::Inversion { group: 1, .. })
    ));
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
