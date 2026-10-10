use super::*;

fn mol2(atoms: &str, bonds: &str, atom_count: usize, bond_count: usize) -> String {
    format!(
        "@<TRIPOS>MOLECULE\ntest\n{atom_count} {bond_count} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n{atoms}@<TRIPOS>BOND\n{bonds}"
    )
}

#[test]
fn mol2_missing_and_malformed_records_fail_with_specific_diagnostics() {
    let valid_atom = "1 C1 0 0 0 C.3 1 LIG 0\n";
    for (input, diagnostic) in [
        ("".to_string(), "MOLECULE"),
        ("@<TRIPOS>MOLECULE\nx\n1 0\n".to_string(), "ATOM"),
        (mol2(valid_atom, "", 0, 0), "no atoms"),
        (mol2(valid_atom, "", 10000, 0), "EOF"),
        (mol2("\n", "", 1, 0), "no info"),
        (mol2("1 C\n", "", 1, 0), "premature"),
        (mol2("1 C 0\n", "", 1, 0), "premature"),
        (mol2("1 C 0 0\n", "", 1, 0), "premature"),
        (mol2("1 C 0 0 0\n", "", 1, 0), "premature"),
        (mol2("1 C x 0 0 C.3\n", "", 1, 0), "coordinates"),
        (mol2("1 C 0 0 0 Qq\n", "", 1, 0), "Element"),
        (mol2("1 C 0 0 0 ANY\n", "", 1, 0), "query"),
        (mol2(valid_atom, "1 1\n", 1, 1), "bond line"),
        (mol2(valid_atom, "1 x 1 1\n", 1, 1), "unsigned"),
        (mol2(valid_atom, "1 1 2 1\n", 1, 1), "index"),
        (mol2(valid_atom, "1 1 1 1\n", 1, 1), "self-bond"),
    ] {
        let error = rdkit_mol_from_mol2_block(&input, true, true, true)
            .err()
            .unwrap();
        assert!(error.to_string().contains(diagnostic), "{input}: {error}");
    }
    for atom in ["Du", "HEV", "HET", "HAL"] {
        let input = mol2(&format!("1 X 0 0 0 {atom}\n"), "", 1, 0);
        assert!(matches!(
            rdkit_mol_from_mol2_block(&input, true, true, true),
            Err(RdkitSmilesError::Unsupported(_))
        ));
    }
    for bond in ["du", "un"] {
        let input = mol2(
            "1 C1 0 0 0 C.3\n2 C2 1.5 0 0 C.3\n",
            &format!("1 1 2 {bond}\n"),
            2,
            1,
        );
        assert!(matches!(
            rdkit_mol_from_mol2_block(&input, true, true, true),
            Err(RdkitSmilesError::Unsupported(_))
        ));
    }
}

#[test]
fn mol2_unity_charge_errors_do_not_silently_guess_charges() {
    let base = mol2("1 N 0 0 0 N.4\n", "", 1, 0);
    for (attr, diagnostic) in [
        ("", "EOF"),
        ("x 1\n", "UnityAtomAttr"),
        ("1 x\n", "UnityAtomAttr"),
        ("1 1\n", "EOF"),
        ("1 1\nAtomExpr nope\n", "formal charge"),
        ("2 1\nAtomExpr 1\n", "index"),
    ] {
        let input = format!("{base}@<TRIPOS>UNITY_ATOM_ATTR\n{attr}");
        let error = rdkit_mol_from_mol2_block(&input, false, false, true)
            .err()
            .unwrap();
        assert!(error.to_string().contains(diagnostic), "{error}");
    }
    let input = format!("{base}@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr 1\n\n");
    let read = rdkit_mol_from_mol2_block(&input, false, false, true).unwrap();
    assert_eq!(read.molecule.atom(chematic_core::AtomIdx(0)).charge, 1);
}

#[test]
fn mol2_lone_pair_records_and_unknown_bonds_are_skipped_explicitly() {
    let atoms = "1 C 0 0 0 C.3\n2 LP 0 0 1 LP\n";
    let input = mol2(atoms, "1 1 2 1\n", 2, 1);
    let read = rdkit_mol_from_mol2_block(&input, false, false, false).unwrap();
    assert_eq!(read.molecule.atom_count(), 1);
    assert_eq!(read.coords, vec![[0., 0., 0.]]);
    assert_eq!(read.molecule.bond_count(), 0);
    let atoms = "1 C 0 0 0 C.3\n2 C 1.5 0 0 C.3\n";
    let input = mol2(atoms, "1 1 2 nc\n", 2, 1);
    assert_eq!(
        rdkit_mol_from_mol2_block(&input, false, false, false)
            .unwrap()
            .molecule
            .bond_count(),
        0
    );
    let input = mol2(atoms, "1 1 2 1\n2 1 2 1\n", 2, 2);
    assert!(
        rdkit_mol_from_mol2_block(&input, false, false, false)
            .err()
            .unwrap()
            .to_string()
            .contains("already exists")
    );
}

#[test]
fn inchi_output_rejects_unmodeled_bonds_without_partial_molecules() {
    let carbon = InchiOutputAtom {
        element: "C".into(),
        num_iso_h: [4, 0, 0, 0],
        ..Default::default()
    };
    let unknown = InchiOutputAtom {
        element: "Qq".into(),
        ..Default::default()
    };
    assert!(
        rdkit_molecule_from_inchi_output(&[unknown], &[])
            .err()
            .unwrap()
            .to_string()
            .contains("element")
    );
    for (bond_type, stereo, diagnostic) in [(9, 0, "bond type"), (1, 1, "2D bond stereo")] {
        let mut left = carbon.clone();
        left.num_iso_h[0] = 3;
        left.bonds = vec![(1, bond_type, stereo)];
        let mut right = carbon.clone();
        right.num_iso_h[0] = 3;
        let error = rdkit_molecule_from_inchi_output(&[left, right], &[])
            .err()
            .unwrap();
        assert!(error.to_string().contains(diagnostic));
    }
}
