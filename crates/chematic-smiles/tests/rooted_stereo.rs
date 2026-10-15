use chematic_core::{
    AtomIdx, BondOrder, Chirality, Molecule, MoleculeBuilder, STEREO_H_SENTINEL, StereoGroup,
    StereoGroupKind, remap_tetrahedral_parity,
};
use chematic_smiles::{ROOTED_SMILES_MAX_ATOMS, RootedSmilesError, parse, write_rooted_checked};
#[path = "support/rooted_bridgehead.rs"]
mod fixture;

fn check_frames(before: &Molecule, text: &str, order: &[AtomIdx]) {
    let after = parse(text).unwrap();
    assert_eq!(after.atom_count(), before.atom_count());
    assert_eq!(after.bond_count(), before.bond_count());
    let mut sorted = order.to_vec();
    sorted.sort_unstable_by_key(|i| i.0);
    assert_eq!(
        sorted,
        (0..before.atom_count())
            .map(|i| AtomIdx(i as u32))
            .collect::<Vec<_>>()
    );
    for (new, &old) in order.iter().enumerate() {
        let original = before.atom(old);
        let written = after.atom(AtomIdx(new as u32));
        assert_eq!(
            (written.element, written.atom_map, written.charge),
            (original.element, original.atom_map, original.charge)
        );
        if !original.chirality.is_tetrahedral() {
            continue;
        }
        let old_frame: [u32; 4] = before
            .stereo_neighbor_order(old)
            .unwrap()
            .try_into()
            .unwrap();
        let mapped: Vec<_> = after
            .stereo_neighbor_order(AtomIdx(new as u32))
            .unwrap()
            .iter()
            .map(|&i| {
                if i == STEREO_H_SENTINEL {
                    i
                } else {
                    order[i as usize].0
                }
            })
            .collect();
        let odd = remap_tetrahedral_parity(old_frame, mapped.try_into().unwrap()).unwrap();
        let expected = match (original.chirality, odd) {
            (Chirality::Clockwise, true) => Chirality::CounterClockwise,
            (Chirality::CounterClockwise, true) => Chirality::Clockwise,
            (tag, _) => tag,
        };
        assert_eq!(
            written.chirality, expected,
            "rooted text {text}, old atom {}",
            old.0
        );
    }
}

#[test]
fn constructed_implicit_h_bridgeheads_preserve_frames_at_every_root() {
    for explicit in [false, true] {
        for flips in 0..4 {
            for rotation in 0..8 {
                for reverse in [false, true] {
                    let mol = fixture::bridgehead(explicit, flips, rotation, reverse);
                    for root in 0..7 {
                        let (text, order) = write_rooted_checked(&mol, root).unwrap();
                        assert_eq!(order[0], AtomIdx(root as u32));
                        check_frames(&mol, &text, &order);
                    }
                }
            }
        }
    }
}

#[test]
fn parsed_centers_isotopes_disconnected_components_and_ez_keep_their_atoms() {
    for text in [
        "[13CH3:1][C@H:2]([OH:3])[F:4]",
        "[C@H:1]1([F:2])[CH2:3][CH2:4][CH2:5]1.[Cl-:6]",
        "[F:1]/[CH:2]=[CH:3]/[Cl:4]",
        "[N:1]->[Cu:2]",
    ] {
        let mol = parse(text).unwrap();
        for root in 0..mol.atom_count() {
            let (s, order) = write_rooted_checked(&mol, root).unwrap();
            check_frames(&mol, &s, &order);
            let reparsed = parse(&s).unwrap();
            for (i, &old) in order.iter().enumerate() {
                assert_eq!(
                    mol.atom(old).isotope,
                    reparsed.atom(AtomIdx(i as u32)).isotope
                );
            }
        }
    }
}

#[test]
fn missing_duplicate_foreign_and_missing_h_frames_are_refused() {
    let mol = parse("C[C@H](O)F").unwrap();
    let mut b = MoleculeBuilder::from_molecule(&mol);
    b.clear_stereo_neighbor_order(AtomIdx(1));
    assert!(matches!(
        write_rooted_checked(&b.build(), 0),
        Err(RootedSmilesError::MissingStereoFrame { atom: 1 })
    ));
    for order in [
        vec![0, 2, 3],
        vec![0, 2, 2, STEREO_H_SENTINEL],
        vec![0, 2, 99, STEREO_H_SENTINEL],
        vec![0, 2, 3, 1],
    ] {
        let mut m = mol.clone();
        m.set_stereo_neighbor_order(AtomIdx(1), order);
        assert!(matches!(
            write_rooted_checked(&m, 0),
            Err(RootedSmilesError::InvalidStereoFrame { atom: 1 })
        ));
    }
}

#[test]
fn other_stereo_and_query_chemistry_are_refused_without_guessing() {
    for text in ["F[Pt@SP1](Cl)(Br)I", "FC=[C@AL1]=C(Cl)Br", "C[S@](=O)CC"] {
        let m = parse(text).unwrap();
        assert!(write_rooted_checked(&m, 0).is_err(), "{text}");
    }
    let mut m = parse("C[C@H](O)F").unwrap();
    m.add_stereo_group(StereoGroup::new(StereoGroupKind::And(1), vec![AtomIdx(1)]));
    assert!(matches!(
        write_rooted_checked(&m, 0),
        Err(RootedSmilesError::UnsupportedStereo { atom: None })
    ));
    for order in [
        BondOrder::Zero,
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ] {
        let mut m = parse("CC").unwrap();
        m.set_bond_order(chematic_core::BondIdx(0), order);
        assert!(matches!(
            write_rooted_checked(&m, 0),
            Err(RootedSmilesError::UnsupportedBond { .. })
        ));
    }
}

#[test]
fn invalid_roots_and_atom_limits_refuse_before_traversal() {
    for m in [MoleculeBuilder::new().build(), parse("CC").unwrap()] {
        for root in [m.atom_count(), usize::MAX] {
            assert!(matches!(
                write_rooted_checked(&m, root),
                Err(RootedSmilesError::InvalidRoot { .. })
            ));
        }
    }
    let at_cap = parse(&"C".repeat(ROOTED_SMILES_MAX_ATOMS)).unwrap();
    assert!(write_rooted_checked(&at_cap, 0).is_ok());
    let above = parse(&"C".repeat(ROOTED_SMILES_MAX_ATOMS + 1)).unwrap();
    assert!(matches!(
        write_rooted_checked(&above, 0),
        Err(RootedSmilesError::ResourceLimit { .. })
    ));
}
