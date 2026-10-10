use crate::{
    InchiOutputAtom, InchiOutputStereo0D, rdkit_canonical_smiles, rdkit_molecule_from_inchi_output,
};
use chematic_core::{AtomIdx, BondOrder};
fn atoms(source: &str) -> Vec<InchiOutputAtom> {
    let mol = crate::parse(source).unwrap();
    mol.atoms()
        .map(|(idx, a)| InchiOutputAtom {
            element: a.element.symbol().into(),
            charge: a.charge,
            num_iso_h: [chematic_core::implicit_hcount(&mol, idx) as i8, 0, 0, 0],
            bonds: mol
                .neighbors(idx)
                .map(|(nb, bi)| {
                    (
                        nb.0 as usize,
                        match mol.bond(bi).order {
                            BondOrder::Single => 1,
                            BondOrder::Double => 2,
                            BondOrder::Triple => 3,
                            _ => panic!("fixture bond"),
                        },
                        0,
                    )
                })
                .collect(),
            ..Default::default()
        })
        .collect()
}
fn canonical(atoms: &[InchiOutputAtom], stereo: &[InchiOutputStereo0D]) -> String {
    rdkit_canonical_smiles(&rdkit_molecule_from_inchi_output(atoms, stereo).unwrap()).unwrap()
}
#[test]
fn inchi_raw_output_rejects_unknown_elements_bond_types_and_unbonded_stereo_ligands() {
    let mut bad = atoms("CC");
    bad[0].element = "Xx".into();
    assert!(
        rdkit_molecule_from_inchi_output(&bad, &[])
            .err()
            .unwrap()
            .to_string()
            .contains("element Xx")
    );
    for (bond_type, bond_stereo, context) in [
        (0, 0, "bond type 0"),
        (5, 0, "bond type 5"),
        (9, 0, "bond type 9"),
        (1, 1, "2D bond stereo"),
    ] {
        let mut bad = atoms("CC");
        bad[0].bonds[0] = (1, bond_type, bond_stereo);
        assert!(
            rdkit_molecule_from_inchi_output(&bad, &[])
                .err()
                .unwrap()
                .to_string()
                .contains(context)
        );
    }
    let raw = atoms("C(F)(Cl)(Br)I.C");
    for (neighbors, context) in [
        ([1, 2, 3, 5], "neighbour not bonded"),
        ([0, 1, 2, 3], "non-bonded ligand"),
    ] {
        let stereo = InchiOutputStereo0D {
            neighbor: neighbors,
            central_atom: 0,
            stereo_type: 2,
            parity: 1,
        };
        assert!(
            rdkit_molecule_from_inchi_output(&raw, &[stereo])
                .err()
                .unwrap()
                .to_string()
                .contains(context)
        );
    }
}
#[test]
fn inchi_raw_tetrahedral_parity_tracks_ligand_permutations_and_implicit_h() {
    for (source, neighbors) in [
        ("C(F)(Cl)(Br)I", [1, 2, 3, 4]),
        ("C(F)(Cl)Br", [0, 1, 2, 3]),
    ] {
        let raw = atoms(source);
        let make = |neighbor, parity| InchiOutputStereo0D {
            neighbor,
            central_atom: 0,
            stereo_type: 2,
            parity,
        };
        let odd = canonical(&raw, &[make(neighbors, 1)]);
        let even = canonical(&raw, &[make(neighbors, 2)]);
        assert_ne!(odd, even);
        let mut flipped = neighbors;
        flipped.swap(2, 3);
        assert_eq!(canonical(&raw, &[make(flipped, 2)]), odd);
        for parity in [0, 4] {
            assert_eq!(
                canonical(&raw, &[make(neighbors, parity)]),
                canonical(&raw, &[])
            );
        }
        assert_eq!(
            rdkit_molecule_from_inchi_output(&raw, &[make(neighbors, 1)])
                .unwrap()
                .atom_count(),
            raw.len()
        );
    }
}
#[test]
fn inchi_raw_double_bond_parity_and_high_priority_substituent_swaps_are_consistent() {
    let raw = atoms("FC(Cl)=C(Br)I");
    let make = |neighbor, parity| InchiOutputStereo0D {
        neighbor,
        central_atom: -1,
        stereo_type: 1,
        parity,
    };
    let odd = canonical(&raw, &[make([0, 1, 3, 4], 1)]);
    let even = canonical(&raw, &[make([0, 1, 3, 4], 2)]);
    assert_ne!(odd, even);
    assert_eq!(canonical(&raw, &[make([2, 1, 3, 4], 2)]), odd);
    assert_eq!(canonical(&raw, &[make([2, 1, 3, 5], 1)]), odd);
    let unrelated = make([0, 0, 5, 3], 1);
    assert_eq!(canonical(&raw, &[unrelated]), canonical(&raw, &[]));
}
#[test]
fn inchi_raw_isotope_hydrogen_and_radical_metadata_survive_conversion() {
    for (slot, isotope) in [(1, 1), (2, 2), (3, 3)] {
        let mut raw = atoms("C");
        raw[0].num_iso_h = [3, 0, 0, 0];
        raw[0].num_iso_h[slot] = 1;
        let result = rdkit_molecule_from_inchi_output(&raw, &[]).unwrap();
        assert_eq!(result.atom_count(), 2);
        assert_eq!(result.atom(AtomIdx(1)).isotope, Some(isotope));
        assert_eq!(result.bond_count(), 1);
    }
    for (shift, isotope) in [(10001, 13), (9999, 11)] {
        let mut raw = atoms("C");
        raw[0].isotopic_mass = shift;
        let result = rdkit_molecule_from_inchi_output(&raw, &[]).unwrap();
        assert_eq!(result.atom(AtomIdx(0)).isotope, Some(isotope));
    }
    for (source, radical, hydrogens, expected) in [("C", 2, 3, "[CH3]"), ("O", 3, 0, "[O]")] {
        let mut raw = atoms(source);
        raw[0].radical = radical;
        raw[0].num_iso_h[0] = hydrogens;
        assert_eq!(canonical(&raw, &[]), expected);
    }
}

#[test]
fn unusual_inchi_valence_cleanup_preserves_connectivity_and_localizes_charge() {
    // These are deliberately pre-sanitization InChI output graphs. The
    // reader's cleanup transfers bond order and charge without deleting atoms
    // or changing connectivity. Compare to explicit, valence-correct products.
    for (source, expected) in [
        ("CN(=C)=N", "C=[N+](C)[NH-]"),
        ("CN(=C)=[NH2+]", "C=[N+](C)N"),
        ("CN(=C)=O", "C=[N+](C)[O-]"),
        ("CN(=C)=S", "C=[N+](C)[S-]"),
        ("[S-](=N)(=C)#N", "C=[S](=[N])=N"),
        ("Cl#S", "ClS"),
    ] {
        let raw = atoms(source);
        let actual = rdkit_molecule_from_inchi_output(&raw, &[]).unwrap();
        let original = crate::parse(source).unwrap();
        assert_eq!(actual.atom_count(), original.atom_count(), "{source}");
        assert_eq!(actual.bond_count(), original.bond_count(), "{source}");
        for (idx, atom) in original.atoms() {
            assert_eq!(actual.atom(idx).element, atom.element, "{source}");
            let mut before: Vec<_> = original.neighbors(idx).map(|(n, _)| n).collect();
            let mut after: Vec<_> = actual.neighbors(idx).map(|(n, _)| n).collect();
            before.sort_unstable();
            after.sort_unstable();
            assert_eq!(after, before, "{source}");
        }
        assert_eq!(
            rdkit_canonical_smiles(&actual).unwrap(),
            rdkit_canonical_smiles(&crate::parse(expected).unwrap()).unwrap(),
            "{source}"
        );
    }
}
