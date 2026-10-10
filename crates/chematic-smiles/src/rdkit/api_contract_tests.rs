//! Public API contracts: atom numbering, geometry and chemical identity.
use super::*;
use crate::{parse, parse_template};

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-6, "{actual} != {expected}");
}

#[test]
fn topological_and_cartesian_distance_contracts() {
    let mol = parse("C=CO").unwrap();
    assert_eq!(
        rdkit_distance_matrix(&mol, false, false).unwrap(),
        vec![vec![0., 1., 2.], vec![1., 0., 1.], vec![2., 1., 0.]]
    );
    assert_eq!(
        rdkit_distance_matrix(&mol, true, false).unwrap(),
        vec![vec![0., 0.5, 1.5], vec![0.5, 0., 1.], vec![1.5, 1., 0.]]
    );
    let weighted = rdkit_distance_matrix(&mol, false, true).unwrap();
    close(weighted[0][0], 1.);
    close(weighted[2][2], 0.75);
    let coords = [[0., 0., 0.], [3., 0., 0.], [3., 4., 0.]];
    let spatial = rdkit_distance_matrix_3d(&mol, &coords, true).unwrap();
    close(spatial[0][1], 3.);
    close(spatial[1][2], 4.);
    close(spatial[0][2], 5.);
    close(spatial[2][2], 0.75);
    assert!(
        rdkit_distance_matrix_3d(&mol, &coords[..2], false)
            .unwrap_err()
            .to_string()
            .contains("coordinates")
    );
    let fragments = parse("C.O").unwrap();
    close(
        rdkit_distance_matrix(&fragments, true, false).unwrap()[0][1],
        1e8,
    );
    let aromatic = parse("c1ccccc1").unwrap();
    close(
        rdkit_distance_matrix(&aromatic, true, false).unwrap()[0][1],
        2. / 3.,
    );
}

#[test]
fn chemistry_diagnostics_preserve_problem_atom_indices() {
    let valid = parse("CCO").unwrap();
    assert!(rdkit_detect_chemistry_problems(&valid).unwrap().is_empty());
    assert!(
        rdkit_detect_chemistry_problems_sanitized(&valid)
            .unwrap()
            .is_empty()
    );
    let ring = parse("c1cccc1").unwrap();
    let problems = rdkit_detect_chemistry_problems(&ring).unwrap();
    assert_eq!(problems.len(), 1);
    assert_eq!(problems[0].kind, "KekulizeException");
    assert_eq!(problems[0].atoms, vec![0, 1, 2, 3, 4]);
    assert!(!problems[0].message.is_empty());
    assert!(rdkit_detect_chemistry_problems_sanitized(&ring).is_err());
    let invalid_valence = parse_template("CO(C)C").unwrap();
    let problems = rdkit_detect_chemistry_problems(&invalid_valence).unwrap();
    assert!(
        problems
            .iter()
            .any(|p| p.kind == "AtomValenceException" && p.atoms == [1])
    );
}

#[test]
fn random_smiles_preserve_stereo_and_seed_reproducibility() {
    let mol = parse("N[C@@H](C)C(=O)Oc1ccccc1").unwrap();
    let params = RdkitSmilesParams::default();
    let expected = rdkit_canonical_smiles(&mol).unwrap();
    let first = rdkit_random_smiles(&mol, 24, 42, &params).unwrap();
    assert_eq!(first, rdkit_random_smiles(&mol, 24, 42, &params).unwrap());
    assert_eq!(first.len(), 24);
    for smiles in first {
        assert_eq!(
            rdkit_canonical_smiles(&parse(&smiles).unwrap()).unwrap(),
            expected,
            "{smiles}"
        );
    }
    assert!(
        rdkit_random_smiles(&mol, 0, 42, &params)
            .unwrap()
            .is_empty()
    );
    let rooted = RdkitSmilesParams {
        canonical: false,
        rooted_at_atom: Some(0),
        ..params
    };
    for smiles in rdkit_random_smiles_with(&mol, 8, 19, &rooted).unwrap() {
        assert!(smiles.starts_with('N'));
        assert_eq!(
            rdkit_canonical_smiles(&parse(&smiles).unwrap()).unwrap(),
            expected
        );
    }
}

#[test]
fn hybridization_and_pi_counts_use_original_atom_order() {
    for (smiles, hybrid, pi) in [
        ("CC#N", vec![3, 1, 1], vec![0, 2, 2]),
        ("C=O", vec![2, 2], vec![1, 1]),
        ("c1ccccc1", vec![2; 6], vec![1; 6]),
    ] {
        let mol = parse(smiles).unwrap();
        assert_eq!(rdkit_hybridizations(&mol).unwrap(), hybrid, "{smiles}");
        assert_eq!(rdkit_num_pi_electrons(&mol).unwrap(), pi, "{smiles}");
    }
    let explicit_h = parse("[H]C").unwrap();
    assert!(rdkit_hybridizations(&explicit_h).is_none());
    assert!(rdkit_num_pi_electrons(&explicit_h).is_none());
}

#[test]
fn hydrogen_addition_and_first_match_typing_contracts() {
    let methane = parse("C").unwrap();
    let (types, heavy) = rdkit_addhs_first_pattern(&methane, &["[#6]", "[#1]", "*"]).unwrap();
    assert_eq!(heavy, 1);
    assert_eq!(types, vec![Some(0), Some(1), Some(1), Some(1), Some(1)]);
    assert_eq!(
        rdkit_addhs_first_pattern(&methane, &["*", "[#6]"])
            .unwrap()
            .0,
        vec![Some(0); 5]
    );
    assert_eq!(
        rdkit_addhs_first_pattern(&methane, &["[#8]"]).unwrap().0,
        vec![None; 5]
    );
    let expanded = rdkit_added_hs_molecule(&methane).unwrap();
    assert_eq!(expanded.atom_count(), 5);
    assert_eq!(expanded.bond_count(), 4);
    assert_eq!(rdkit_num_atoms(&methane).unwrap(), 1);
    let explicit = parse("[H]C([H])([H])[H]").unwrap();
    let (suppressed, map) = rdkit_hydrogen_suppressed_with_map(&explicit).unwrap();
    assert_eq!(map, vec![1]);
    assert_eq!(suppressed.atom_count(), 1);
    assert_eq!(rdkit_canonical_smiles(&suppressed).unwrap(), "C");
    assert!(rdkit_hydrogen_suppressed_with_map(&methane).is_none());
}

#[test]
fn alignment_recovers_a_rigid_transform_and_checks_dimensions() {
    let mol = parse("CCO").unwrap();
    let reference = [[0., 0., 0.], [1.5, 0., 0.], [1.5, 1.2, 0.3]];
    let probe = reference.map(|[x, y, z]| [-y + 2., x - 1., z + 0.5]);
    for map in [None, Some(&[(0, 0), (1, 1), (2, 2)][..])] {
        let fit = rdkit_align_mol(
            &mol,
            &probe,
            &reference,
            map,
            Some(&[1., 2., 3.]),
            false,
            100,
        )
        .unwrap();
        close(fit.rmsd, 0.);
        for (p, r) in probe.iter().zip(reference) {
            for (actual, expected) in fit.transform.transform_point(*p).into_iter().zip(r) {
                close(actual, expected);
            }
        }
        assert_eq!(fit.transform.rows()[3], [0., 0., 0., 1.]);
    }
    close(
        rdkit_best_rms(&mol, &probe, &reference, 100, true, None)
            .unwrap()
            .rmsd,
        0.,
    );
    close(
        rdkit_align_points(&reference, &probe, None, false, 100)
            .unwrap()
            .0,
        0.,
    );
    assert!(rdkit_align_mol(&mol, &probe[..2], &reference, None, None, false, 100).is_err());
    assert!(rdkit_align_points(&reference, &probe, Some(&[1.]), false, 100).is_err());
    let translated = reference.map(|[x, y, z]| [x + 2., y - 1., z + 0.5]);
    close(
        rdkit_calc_rms(&mol, &translated, &reference, 100, false, None).unwrap(),
        5.25_f64.sqrt(),
    );
}

#[test]
fn fragment_writer_respects_selected_atoms_and_bonds() {
    let mol = parse("OC(=O)c1ccccc1").unwrap();
    let params = RdkitSmilesParams::default();
    assert_eq!(
        rdkit_fragment_smiles(&mol, &[0, 1, 2], None, &params).unwrap(),
        "O=CO"
    );
    assert_eq!(
        rdkit_fragment_smiles(&mol, &[0, 1, 2], Some(&[]), &params).unwrap(),
        "C.O.O"
    );
    assert!(rdkit_fragment_smiles(&mol, &[99], None, &params).is_err());
    assert_eq!(
        rdkit_cx_smiles(&parse("CCO").unwrap(), &params).unwrap(),
        "CCO"
    );
}
