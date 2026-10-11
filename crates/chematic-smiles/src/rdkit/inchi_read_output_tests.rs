use super::*;

#[test]
fn inchi_output_preserves_isotopes_radicals_and_stereo() {
    let alkene = [
        InchiOutputAtom {
            element: "C".into(),
            bonds: vec![(1, 1, 0)],
            num_iso_h: [3, 0, 0, 0],
            ..Default::default()
        },
        InchiOutputAtom {
            element: "C".into(),
            bonds: vec![(0, 1, 0), (2, 2, 0)],
            num_iso_h: [1, 0, 0, 0],
            ..Default::default()
        },
        InchiOutputAtom {
            element: "C".into(),
            bonds: vec![(1, 2, 0), (3, 1, 0)],
            num_iso_h: [1, 0, 0, 0],
            ..Default::default()
        },
        InchiOutputAtom {
            element: "C".into(),
            bonds: vec![(2, 1, 0)],
            num_iso_h: [3, 0, 0, 0],
            ..Default::default()
        },
    ];
    let stereo = [InchiOutputStereo0D {
        neighbor: [0, 1, 2, 3],
        central_atom: -1,
        stereo_type: INCHI_STEREO_TYPE_DOUBLE_BOND,
        parity: INCHI_PARITY_ODD,
    }];
    let mol = rdkit_molecule_from_inchi_output(&alkene, &stereo).unwrap();
    let written = crate::canonical_smiles(&mol);
    assert!(written.contains('/') || written.contains('\\'));

    let isotopic = [InchiOutputAtom {
        element: "C".into(),
        isotopic_mass: (ISOTOPIC_SHIFT_FLAG + 1) as i16,
        radical: 3,
        num_iso_h: [1, 0, 1, 0],
        ..Default::default()
    }];
    let mol = rdkit_molecule_from_inchi_output(&isotopic, &[]).unwrap();
    assert_eq!(mol.atom_count(), 2);
    assert!(mol.atoms().any(|(_, a)| a.isotope == Some(2)));
}

#[test]
fn inchi_output_rejects_bond_stereo_and_bad_tetrahedral_neighbors() {
    let bad_bond_stereo = [
        InchiOutputAtom {
            element: "C".into(),
            bonds: vec![(1, 1, 1)],
            ..Default::default()
        },
        InchiOutputAtom {
            element: "C".into(),
            ..Default::default()
        },
    ];
    assert!(rdkit_molecule_from_inchi_output(&bad_bond_stereo, &[]).is_err());

    let tetra = [
        InchiOutputAtom {
            element: "C".into(),
            bonds: vec![(1, 1, 0)],
            ..Default::default()
        },
        InchiOutputAtom {
            element: "F".into(),
            bonds: vec![(0, 1, 0)],
            ..Default::default()
        },
    ];
    let stereo = [InchiOutputStereo0D {
        neighbor: [0, 1, 2, 3],
        central_atom: 0,
        stereo_type: INCHI_STEREO_TYPE_TETRAHEDRAL,
        parity: INCHI_PARITY_EVEN,
    }];
    assert!(rdkit_molecule_from_inchi_output(&tetra, &stereo).is_err());
}
