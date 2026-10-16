use chematic_chem::{
    TetrahedralCleanupError as Error, remove_redundant_tetrahedral_stereo as clean,
};
use chematic_core::{
    Atom, AtomIdx, BondOrder, Chirality, Element, MoleculeBuilder, StereoGroup, StereoGroupKind,
};
use chematic_smiles::{canonical_smiles, parse};

fn cases() -> serde_json::Value {
    serde_json::from_str(include_str!(
        "../../../validation/rdkit_issues/tetrahedral-cleanup/cases.json"
    ))
    .unwrap()
}
#[test]
fn original_controls_and_reported_bridge_have_expected_decisions() {
    for c in cases()["cases"].as_array().unwrap() {
        let m = parse(c["smiles"].as_str().unwrap()).unwrap();
        let a = AtomIdx(c["atom_index"].as_u64().unwrap() as u32);
        let r = clean(&m, a, 1_000_000).unwrap();
        assert_eq!(r.removed, c["removed"].as_bool().unwrap(), "{}", c["id"]);
        assert!(r.search_visits > 0);
        assert_eq!(r.witness.is_some(), r.removed);
        if let Some(w) = r.witness {
            let mut v = w.clone();
            v.sort_unstable();
            assert_eq!(
                v,
                (0..m.atom_count())
                    .map(|i| AtomIdx(i as u32))
                    .collect::<Vec<_>>()
            );
            assert_eq!(w[a.0 as usize], a);
        }
        assert_ne!(m.atom(a).chirality, Chirality::None);
    }
}
#[test]
fn source_graph_and_unselected_frames_are_unchanged_and_cip_is_cleared() {
    let mut m = parse("C=C1C[C@H]2CC(=C)C[C@@H](C1)[C@H]2C").unwrap();
    for (idx, _) in m.clone().atoms() {
        m.set_tag(idx, Some(7));
        m.set_cip_code(idx, Some(chematic_core::CipCode::R));
    }
    let before = canonical_smiles(&m);
    let r = clean(&m, AtomIdx(10), 1_000_000).unwrap();
    assert!(r.removed);
    assert_eq!(canonical_smiles(&m), before);
    assert_eq!(m.atom_count(), r.molecule.atom_count());
    assert_eq!(m.bond_count(), r.molecule.bond_count());
    for (i, a) in m.atoms() {
        let b = r.molecule.atom(i);
        assert_eq!(
            (
                a.element,
                a.charge,
                a.isotope,
                a.hydrogen_count,
                a.aromatic,
                a.atom_map
            ),
            (
                b.element,
                b.charge,
                b.isotope,
                b.hydrogen_count,
                b.aromatic,
                b.atom_map
            )
        );
        assert_eq!(m.atom_tag(i), r.molecule.atom_tag(i));
        assert_eq!(b.cip_code, None);
        if i == AtomIdx(10) {
            assert_eq!(b.chirality, Chirality::None);
            assert!(r.molecule.stereo_neighbor_order(i).is_none());
        } else {
            assert_eq!(a.chirality, b.chirality);
            assert_eq!(
                m.stereo_neighbor_order(i),
                r.molecule.stereo_neighbor_order(i)
            );
        }
        assert!(a.cip_code.is_some());
    }
    for (i, b) in m.bonds() {
        let c = r.molecule.bond(i);
        assert_eq!((b.atom1, b.atom2, b.order), (c.atom1, c.atom2, c.order));
    }
}
#[test]
fn source_labels_can_block_an_odd_ligand_exchange() {
    let mut m = parse("F[C@H](F)Cl").unwrap();
    assert!(clean(&m, AtomIdx(1), 1_000_000).unwrap().removed);
    m.set_tag(AtomIdx(0), Some(3));
    let r = clean(&m, AtomIdx(1), 1_000_000).unwrap();
    assert!(!r.removed);
    assert_eq!(canonical_smiles(&m), canonical_smiles(&r.molecule));
    assert_eq!(m.atom_tag(AtomIdx(0)), r.molecule.atom_tag(AtomIdx(0)));
}
#[test]
fn frame_removal_is_local_and_does_not_erase_other_configuration() {
    let mut m = parse("C[C@@H]([C@@H](C)O)[C@@H](C)O").unwrap();
    let other = m.stereo_neighbor_order(AtomIdx(2)).unwrap().to_vec();
    m.set_chirality(AtomIdx(1), Chirality::None);
    m.clear_stereo_neighbor_order(AtomIdx(1));
    assert!(m.stereo_neighbor_order(AtomIdx(1)).is_none());
    assert_eq!(m.stereo_neighbor_order(AtomIdx(2)).unwrap(), other);
    m.clear_stereo_neighbor_order(AtomIdx(1));
    assert!(m.stereo_neighbor_order(AtomIdx(1)).is_none());
}
#[test]
fn budget_exhaustion_never_returns_retention_or_a_partial_cleanup() {
    let m = parse("F[C@H](F)Cl").unwrap();
    let before = canonical_smiles(&m);
    let e = clean(&m, AtomIdx(1), 1).err().unwrap();
    assert!(matches!(
        e,
        Error::Proof(chematic_chem::MesoStereoGroupError::BudgetExceeded { .. })
    ));
    assert_eq!(canonical_smiles(&m), before);
    for limit in [0, 10_000_001, u64::MAX] {
        assert!(matches!(
            clean(&m, AtomIdx(1), limit),
            Err(Error::InvalidVisitLimit)
        ));
    }
}
#[test]
fn invalid_targets_and_enhanced_groups_are_refused_before_proof() {
    let m = parse("F[C@H](F)Cl").unwrap();
    assert!(matches!(
        clean(&m, AtomIdx(99), 1000),
        Err(Error::InvalidAtom { .. })
    ));
    assert!(matches!(
        clean(&m, AtomIdx(0), 1000),
        Err(Error::UnspecifiedAtom(_))
    ));
    let g = m
        .with_checked_stereo_groups(vec![(
            StereoGroup::new(StereoGroupKind::And(1), vec![AtomIdx(1)]),
            vec![],
        )])
        .unwrap();
    assert!(matches!(
        clean(&g, AtomIdx(1), 1000),
        Err(Error::EnhancedGroups)
    ));
}
#[test]
fn incomplete_unspecified_axial_and_unsupported_graph_semantics_are_refused() {
    let mut incomplete = parse("F[C@H](F)Cl").unwrap();
    incomplete.clear_stereo_neighbor_order(AtomIdx(1));
    assert!(matches!(
        clean(&incomplete, AtomIdx(1), 10000),
        Err(Error::Proof(
            chematic_chem::MesoStereoGroupError::IncompleteFrame(_)
        ))
    ));
    for (s, a) in [
        ("F[C@H](F)C(Cl)(Br)I", 1),
        ("F[C@H](F)*", 1),
        ("FC=[C@AL1]=C(Cl)Br", 2),
        ("F[C@H](F)N->[Cu]", 1),
        ("C/[C@H](C)F", 1),
    ] {
        let m = parse(s).unwrap();
        assert!(clean(&m, AtomIdx(a), 10000).is_err(), "accepted {s}");
    }
}
#[test]
fn graph_caps_precede_preparation_and_the_atom_edge_is_accepted() {
    let mut b = MoleculeBuilder::new();
    for _ in 0..65 {
        b.add_atom(Atom::new(Element::C));
    }
    assert!(matches!(
        clean(&b.build(), AtomIdx(0), 1000),
        Err(Error::GraphLimit)
    ));
    let mut b = MoleculeBuilder::new();
    for _ in 0..17 {
        b.add_atom(Atom::new(Element::C));
    }
    for i in 0..17 {
        for j in i + 1..17 {
            b.add_bond(AtomIdx(i), AtomIdx(j), BondOrder::Single)
                .unwrap();
        }
    }
    assert!(matches!(
        clean(&b.build(), AtomIdx(0), 1000),
        Err(Error::GraphLimit)
    ));
    let tail = (1..=60).map(|i| format!(".[CH4:{i}]")).collect::<String>();
    let m = parse(&format!("F[C@H](F)Cl{tail}")).unwrap();
    assert_eq!(m.atom_count(), 64);
    assert!(clean(&m, AtomIdx(1), 1_000_000).unwrap().removed);
}
#[test]
fn public_error_messages_are_actionable() {
    for e in [
        Error::InvalidAtom {
            atom: AtomIdx(4),
            atoms: 3,
        },
        Error::UnspecifiedAtom(AtomIdx(0)),
        Error::EnhancedGroups,
        Error::GraphLimit,
        Error::InvalidVisitLimit,
        Error::Proof(chematic_chem::MesoStereoGroupError::UnsupportedDirections),
    ] {
        assert!(!e.to_string().is_empty());
    }
}
