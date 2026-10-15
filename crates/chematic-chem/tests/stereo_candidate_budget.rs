use chematic_chem::{add_hydrogens, enumerate_stereoisomers};
use chematic_core::{AtomIdx, BondIdx, RGroupLabel, StereoGroup, StereoGroupKind};
use chematic_smiles::{canonical_smiles, parse};
use std::collections::HashSet;

fn keys(s: &str) -> HashSet<String> {
    enumerate_stereoisomers(&parse(s).unwrap())
        .iter()
        .map(canonical_smiles)
        .collect()
}

#[test]
fn aminoadamantane_has_four_unique_stereochemical_assignments() {
    let v = keys("NC12CC3CC(CC(C3)C1)C2");
    assert_eq!(v.len(), 4);
    for s in v {
        assert_eq!(canonical_smiles(&parse(&s).unwrap()), s);
    }
}
#[test]
fn inositol_keeps_all_nine_stereoisomers() {
    let v = keys("OC1C(O)C(O)C(O)C(O)C1O");
    assert_eq!(v.len(), 9);
    for s in v {
        assert_eq!(canonical_smiles(&parse(&s).unwrap()), s);
    }
}
#[test]
fn methylene_rich_alkanes_do_not_consume_the_budget() {
    for n in [8, 16, 32] {
        assert_eq!(keys(&"C".repeat(n)).len(), 1);
    }
    assert_eq!(keys("CCCCCCCCCCCCCCCC(O)F").len(), 2);
}
#[test]
fn explicit_and_mixed_identical_hydrogens_do_not_consume_the_budget() {
    let implicit = parse("NC12CC3CC(CC(C3)C1)C2").unwrap();
    let explicit = add_hydrogens(&implicit);
    assert_eq!(enumerate_stereoisomers(&explicit).len(), 4);
    // Two identical explicit H plus a non-node H cannot be stereogenic.
    let mixed = "[H][CH]([H])C";
    assert_eq!(keys(mixed).len(), 1);
}
#[test]
fn isotope_distinguished_terminal_hydrogens_remain_stereogenic() {
    assert_eq!(keys("[H]C([2H])(F)Cl").len(), 2);
    assert_eq!(keys("[2H]C(F)Cl").len(), 2);
}
#[test]
fn identical_terminal_halogen_and_methyl_groups_are_excluded() {
    assert_eq!(keys("FC(F)(Br)Cl").len(), 1);
    assert_eq!(keys("CC(C)(Br)Cl").len(), 1);
    let repeated = "FC(F)(Br)Cl.".repeat(8) + "C";
    assert_eq!(keys(&repeated).len(), 1);
}
#[test]
fn atom_maps_and_nonchemical_tags_do_not_make_identical_ligands_different() {
    let mut m = parse("[H:1]C([H:2])(F)Cl").unwrap();
    m.set_tag(AtomIdx(0), Some(1));
    m.set_tag(AtomIdx(2), Some(2));
    let v = enumerate_stereoisomers(&m);
    assert_eq!(v.len(), 1);
    assert_eq!(
        v[0].atom(AtomIdx(1)).chirality,
        chematic_core::Chirality::None
    );
    assert_eq!(v[0].atom_tag(AtomIdx(0)).unwrap().get(), 1);
}
#[test]
fn seven_real_candidate_centers_still_respect_the_existing_cap() {
    let s = std::iter::repeat_n("FC(Cl)(Br)I", 7)
        .collect::<Vec<_>>()
        .join(".");
    assert!(keys(&s).is_empty());
}
#[test]
fn six_real_candidate_centers_are_enumerated_completely() {
    // The distinct atom maps label the disconnected centers for this test,
    // making all 64 indexed assignments observable in the canonical key.
    let s = (1..=6)
        .map(|i| format!("F[CH:{i}](Cl)Br"))
        .collect::<Vec<_>>()
        .join(".");
    assert_eq!(keys(&s).len(), 64);
}
#[test]
fn preserves_full_graph_metadata_in_both_copy_and_assignment_paths() {
    for s in ["FC(Cl)(Br)I.*.C1CC1", "N[C@@H](C)C(=O)O.*.C1CC1"] {
        let mut m = parse(s).unwrap();
        let star = m.atoms().find(|(_, a)| a.wildcard).unwrap().0;
        m.set_r_group(star, RGroupLabel::numbered(3).unwrap());
        m.set_tag(AtomIdx(0), Some(71));
        m.set_stereo_groups_with_bonds(vec![(
            StereoGroup::new(StereoGroupKind::Or(2), vec![AtomIdx(1)]),
            vec![BondIdx(0)],
        )]);
        let expected = canonical_smiles(&m);
        for out in enumerate_stereoisomers(&m) {
            assert_eq!(out.atom_count(), m.atom_count());
            assert_eq!(out.bond_count(), m.bond_count());
            assert_eq!(out.atom_tag(AtomIdx(0)).unwrap().get(), 71);
            assert_eq!(out.r_group_label(star).unwrap().number(), Some(3));
            assert_eq!(out.stereo_groups(), m.stereo_groups());
            assert_eq!(
                out.stereo_group_bond_indices(0),
                m.stereo_group_bond_indices(0)
            );
            for (b, bond) in m.bonds() {
                assert_eq!(out.bond(b), bond);
                assert_eq!(
                    out.smiles_ring_closure_metadata(b),
                    m.smiles_ring_closure_metadata(b)
                );
            }
            if s.starts_with('N') {
                assert_eq!(canonical_smiles(&out), expected);
            }
        }
    }
}
#[test]
fn results_do_not_mutate_original_chiral_state() {
    let m = parse("FC(Cl)(Br)I").unwrap();
    let before = canonical_smiles(&m);
    assert_eq!(enumerate_stereoisomers(&m).len(), 2);
    assert_eq!(canonical_smiles(&m), before);
    assert_eq!(m.atom(AtomIdx(1)).chirality, chematic_core::Chirality::None);
}
#[test]
fn aromatic_multiple_bonds_and_wildcards_do_not_get_tetrahedral_tags() {
    for s in [
        "c1ccccc1",
        "C=C(F)Cl",
        "C#CC",
        "FC(~Cl)(Br)I",
        "F*(Cl)(Br)I",
    ] {
        let m = parse(s).unwrap();
        for out in enumerate_stereoisomers(&m) {
            for (i, a) in m.atoms() {
                if a.aromatic
                    || a.wildcard
                    || m.neighbors(i).any(|(_, b)| {
                        !matches!(
                            m.bond(b).order,
                            chematic_core::BondOrder::Single
                                | chematic_core::BondOrder::Up
                                | chematic_core::BondOrder::Down
                        )
                    })
                {
                    assert_eq!(out.atom(i).chirality, a.chirality);
                }
            }
        }
    }
}

#[test]
fn new_assignments_invalidate_derived_cip_annotations_but_copies_keep_them() {
    let mut source = parse("N[C@@H](C)C(O)F").unwrap();
    source.set_cip_code(AtomIdx(1), Some(chematic_core::CipCode::R));
    source.set_cip_code(AtomIdx(3), Some(chematic_core::CipCode::S));
    let forms = enumerate_stereoisomers(&source);
    assert_eq!(forms.len(), 2);
    for form in forms {
        assert!(form.atoms().all(|(_, atom)| atom.cip_code.is_none()));
        assert_eq!(
            form.atom(AtomIdx(1)).chirality,
            source.atom(AtomIdx(1)).chirality
        );
    }
    assert_eq!(
        source.atom(AtomIdx(1)).cip_code,
        Some(chematic_core::CipCode::R)
    );
    let mut fixed = parse("N[C@@H](C)C(=O)O").unwrap();
    fixed.set_cip_code(AtomIdx(1), Some(chematic_core::CipCode::S));
    assert_eq!(
        enumerate_stereoisomers(&fixed)[0].atom(AtomIdx(1)).cip_code,
        Some(chematic_core::CipCode::S)
    );
}
