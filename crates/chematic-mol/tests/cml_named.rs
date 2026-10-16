use chematic_mol::{CmlNamedError as Error, parse_cml_strict, write_cml_named_checked};

fn ids(names: &[&str]) -> Vec<String> {
    names.iter().map(|s| s.to_string()).collect()
}

#[test]
fn names_and_references_are_written_without_replacement() {
    let mol = chematic_smiles::parse("CO").unwrap();
    let text = write_cml_named_checked(
        &mol,
        &ids(&["C9", "O_alpha-2"]),
        Some(&[[1., 2., 3.], [4., 5., 6.]]),
    )
    .unwrap();
    assert!(text.contains("id=\"C9\""));
    assert!(text.contains("atomRefs2=\"C9 O_alpha-2\""));
    assert!(!text.contains("id=\"a1\""));
    assert_eq!(
        chematic_chem::calc_mol_formula(&parse_cml_strict(&text).unwrap().0),
        "CH4O"
    );
}

#[test]
fn invalid_and_duplicate_ids_are_errors() {
    let mol = chematic_smiles::parse("CC").unwrap();
    assert_eq!(
        write_cml_named_checked(&mol, &ids(&["C1"]), None),
        Err(Error::IdCount)
    );
    for bad in ["", "1C", "C C", "C:1", "C.1", "C\"/>", "炭素"] {
        assert_eq!(
            write_cml_named_checked(&mol, &ids(&["C1", bad]), None),
            Err(Error::InvalidId { index: 1 })
        );
    }
    assert_eq!(
        write_cml_named_checked(&mol, &ids(&["C1", "C1"]), None),
        Err(Error::DuplicateId { index: 1 })
    );
}

#[test]
fn coordinate_contract() {
    let mol = chematic_smiles::parse("C").unwrap();
    assert_eq!(
        write_cml_named_checked(&mol, &ids(&["C1"]), Some(&[])),
        Err(Error::CoordinateCount)
    );
    for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            write_cml_named_checked(&mol, &ids(&["C1"]), Some(&[[0., 0., v]])),
            Err(Error::NonFiniteCoordinate { index: 0 })
        );
    }
    assert!(
        !write_cml_named_checked(&mol, &ids(&["C1"]), None)
            .unwrap()
            .contains("x3=")
    );
}

#[test]
fn charges_isotopes_aromaticity_and_explicit_hydrogens_round_trip() {
    for smiles in [
        "[NH4+]",
        "[2H]O[H]",
        "[13CH3]O",
        "c1cc[nH]c1",
        "C#N",
        "CC(=O)O",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let ids = (0..mol.atom_count())
            .map(|i| format!("Atom{i}"))
            .collect::<Vec<_>>();
        let out = write_cml_named_checked(&mol, &ids, None).unwrap();
        let parsed = parse_cml_strict(&out).unwrap().0;
        assert_eq!(
            chematic_chem::calc_mol_formula(&mol),
            chematic_chem::calc_mol_formula(&parsed),
            "{smiles}"
        );
        assert_eq!(mol.atom_count(), parsed.atom_count());
        for ((_, a), (_, b)) in mol.atoms().zip(parsed.atoms()) {
            assert_eq!(
                (a.element, a.isotope, a.charge, a.aromatic),
                (b.element, b.isotope, b.charge, b.aromatic)
            );
        }
    }
}

#[test]
fn stereo_maps_and_nonmolecular_bonds_are_refused() {
    for smiles in [
        "C[C@H](F)Cl",
        "C/C=C/C",
        "[CH3:7]O",
        "*",
        "N->[Cu+2]",
        "C~C",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let ids = (0..mol.atom_count())
            .map(|i| format!("Atom{i}"))
            .collect::<Vec<_>>();
        assert!(
            matches!(
                write_cml_named_checked(&mol, &ids, None),
                Err(Error::Unsupported(_))
            ),
            "{smiles}"
        );
    }
}

#[test]
fn graph_and_id_size_limits() {
    let empty = chematic_core::MoleculeBuilder::new().build();
    assert_eq!(
        write_cml_named_checked(&empty, &[], None),
        Err(Error::SizeLimit)
    );
    let mol = chematic_smiles::parse("C").unwrap();
    assert!(write_cml_named_checked(&mol, &["C".repeat(128)], None).is_ok());
    assert_eq!(
        write_cml_named_checked(&mol, &["C".repeat(129)], None),
        Err(Error::InvalidId { index: 0 })
    );
}
