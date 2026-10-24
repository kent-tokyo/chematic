//! Open Babel #3005: explicit isotope loss, without changing rigid atom data.
use chematic_core::{Atom, Element, MoleculeBuilder};
use chematic_mol::{PdbqtWriteError, project_pdbqt_isotopes, write_pdbqt, write_pdbqt_checked};
use serde_json::Value;

#[test]
fn labelled_corpus_requires_explicit_projection_and_records_every_label() {
    let cases: Vec<Value> = serde_json::from_str(include_str!(
        "../../../validation/openbabel_issues/pdbqt-isotopes/cases.json"
    ))
    .unwrap();
    for case in cases {
        let mol = chematic_smiles::parse(case["smiles"].as_str().unwrap()).unwrap();
        let before = chematic_smiles::write(&mol);
        let n = mol.atom_count();
        let coords: Vec<_> = (0..n).map(|i| (i as f64, -2.345, 6.789)).collect();
        let charges = vec![-0.123; n];
        let projection = project_pdbqt_isotopes(&mol, &coords, &charges, "labels").unwrap();
        assert_eq!(
            serde_json::to_value(&projection.isotopes_omitted).unwrap(),
            case["isotopes_omitted"]
        );
        assert_eq!(
            projection.pdbqt,
            write_pdbqt(&mol, &coords, &charges, "labels")
        );
        let checked = write_pdbqt_checked(&mol, &coords, &charges, "labels");
        if let Some(loss) = projection.isotopes_omitted.first() {
            assert_eq!(
                checked,
                Err(PdbqtWriteError::IsotopeLoss {
                    atom: loss.source_atom_index,
                    mass_number: loss.mass_number,
                })
            );
        } else {
            assert_eq!(checked.unwrap(), projection.pdbqt);
        }
        assert_eq!(
            chematic_smiles::write(&mol),
            before,
            "source graph must not be mutated"
        );
        let lines: Vec<_> = projection
            .pdbqt
            .lines()
            .filter(|s| s.starts_with("ATOM"))
            .collect();
        assert_eq!(lines.len(), n);
        let (back, back_xyz, back_q) = chematic_mol::parse_pdbqt(&projection.pdbqt).unwrap();
        assert_eq!(back_xyz, coords);
        assert_eq!(back_q, charges);
        for ((idx, atom), line) in back.atoms().zip(lines) {
            assert_eq!(line.len(), 79);
            assert_eq!(line[6..11].trim().parse::<u32>().unwrap(), idx.0 + 1);
            assert_eq!(
                atom.element.atomic_number(),
                case["atomic_numbers"][idx.0 as usize].as_u64().unwrap() as u8
            );
            assert_eq!(atom.isotope, None);
        }
    }
}

#[test]
fn isotope_labels_never_change_hydrogen_types_atom_count_or_supplied_charges() {
    for (labelled, ordinary, types) in [
        (
            "[2H]O[13CH2][3H]",
            "[H]O[CH2][H]",
            vec!["HD", "OA", "C", "H"],
        ),
        ("[3H]N([2H])C", "[H]N([H])C", vec!["HD", "N", "HD", "C"]),
    ] {
        let a = chematic_smiles::parse(labelled).unwrap();
        let b = chematic_smiles::parse(ordinary).unwrap();
        let coords = vec![(1.25, -2.5, 3.75); a.atom_count()];
        let charges = vec![-0.2, 0.3, -0.4, 0.1];
        let report = project_pdbqt_isotopes(&a, &coords, &charges, "same").unwrap();
        assert_eq!(
            report.pdbqt,
            write_pdbqt_checked(&b, &coords, &charges, "same").unwrap()
        );
        let observed: Vec<_> = report
            .pdbqt
            .lines()
            .filter(|s| s.starts_with("ATOM"))
            .map(|s| s[77..79].trim())
            .collect();
        assert_eq!(observed, types);
    }
}

#[test]
fn acetophenone_d8_keeps_all_atoms_without_creating_donor_hydrogens() {
    let labelled = "[2H]C([2H])([2H])C(=O)c1c([2H])c([2H])c([2H])c([2H])c1[2H]";
    let plain = labelled.replace("[2H]", "[H]");
    let a = chematic_smiles::parse(labelled).unwrap();
    let b = chematic_smiles::parse(&plain).unwrap();
    let coords = vec![(1., 2., 3.); 17];
    let charges: Vec<_> = (0..17).map(|i| i as f64 * 0.01).collect();
    let report = project_pdbqt_isotopes(&a, &coords, &charges, "same").unwrap();
    assert_eq!(report.isotopes_omitted.len(), 8);
    assert_eq!(
        report.pdbqt,
        write_pdbqt_checked(&b, &coords, &charges, "same").unwrap()
    );
    let types: Vec<_> = report
        .pdbqt
        .lines()
        .filter(|s| s.starts_with("ATOM"))
        .map(|s| s[77..79].trim())
        .collect();
    assert_eq!(types.len(), 17);
    assert_eq!(types.iter().filter(|&&s| s == "H").count(), 8);
    assert!(!types.contains(&"HD"));
    assert!(!report.pdbqt.contains("TORSDOF"));
}

#[test]
fn explicit_projection_keeps_all_other_checks() {
    let mol = chematic_smiles::parse("[2H]O").unwrap();
    assert!(matches!(
        project_pdbqt_isotopes(&mol, &[], &[0.; 2], ""),
        Err(PdbqtWriteError::CountMismatch { .. })
    ));
    for bad in [f64::NAN, f64::INFINITY, 10000.] {
        assert!(project_pdbqt_isotopes(&mol, &[(bad, 0., 0.); 2], &[0.; 2], "").is_err());
    }
    assert!(project_pdbqt_isotopes(&mol, &[(0., 0., 0.); 2], &[100.; 2], "").is_err());
    assert!(project_pdbqt_isotopes(&mol, &[(0., 0., 0.); 2], &[0.; 2], "\nATOM").is_err());
    let unsupported = chematic_smiles::parse("[10B]").unwrap();
    assert!(matches!(
        project_pdbqt_isotopes(&unsupported, &[(0., 0., 0.)], &[0.], ""),
        Err(PdbqtWriteError::UnsupportedAtom { .. })
    ));
}

#[test]
fn even_explicit_zero_isotope_metadata_requires_projection() {
    let mut b = MoleculeBuilder::new();
    let mut atom = Atom::new(Element::C);
    atom.isotope = Some(0);
    b.add_atom(atom);
    let mol = b.build();
    let error = write_pdbqt_checked(&mol, &[(0., 0., 0.)], &[0.], "").unwrap_err();
    assert_eq!(
        error,
        PdbqtWriteError::IsotopeLoss {
            atom: 0,
            mass_number: 0
        }
    );
    assert!(error.to_string().contains("project_pdbqt_isotopes"));
    let report = project_pdbqt_isotopes(&mol, &[(0., 0., 0.)], &[0.], "").unwrap();
    assert_eq!(report.isotopes_omitted[0].mass_number, 0);
}
