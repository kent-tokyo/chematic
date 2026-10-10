use crate::rdkit_uff::*;
use chematic_smiles::{parse, rdkit_added_hs_molecule};

#[test]
fn uff_gradients_match_energy_and_conserve_total_force() {
    for smiles in [
        "CCCC",
        "CC(=O)N",
        "c1ccccc1",
        "C1CC1",
        "C1CCC1",
        "CS(=O)C",
        "CS(=O)(=O)C",
        "CP(F)(F)(F)F",
        "C[Si](C)(C)C",
    ] {
        let mol = rdkit_added_hs_molecule(&parse(smiles).unwrap()).unwrap();
        let coords: Vec<_> = (0..mol.atom_count())
            .map(|i| {
                let t = i as f64 * 1.7;
                [t.cos() * 2.3, t.sin() * 2.3, i as f64 * 0.7]
            })
            .collect();
        let field = RdkitUffField::new(&mol, &coords, 10., true);
        assert_eq!(field.num_atoms(), coords.len());
        let mut pos: Vec<f64> = coords.iter().flatten().copied().collect();
        let e = field.energy(&pos);
        assert!(e.is_finite(), "{smiles}");
        let terms = field.debug_terms();
        assert!(terms.iter().any(|(kind, _, _)| *kind == "bond"), "{smiles}");
        for (_, atoms, params) in terms {
            assert!(atoms.iter().all(|&i| i < coords.len()));
            assert!(params.iter().all(|v| v.is_finite()));
        }
        let mut grad = vec![0.; pos.len()];
        field.gradient(&pos, &mut grad);
        assert!(grad.iter().all(|v| v.is_finite()), "{smiles}");
        for axis in 0..3 {
            let force: f64 = grad.iter().skip(axis).step_by(3).sum();
            assert!(
                force.abs() < 1e-6 * (1. + grad.iter().map(|v| v.abs()).sum::<f64>()),
                "{smiles}: {force}"
            );
        }
        for i in 0..pos.len() {
            let old = pos[i];
            pos[i] = old + 1e-5;
            let plus = field.energy(&pos);
            pos[i] = old - 1e-5;
            let minus = field.energy(&pos);
            pos[i] = old;
            let numeric = (plus - minus) / 2e-5;
            assert!(
                (numeric - grad[i]).abs() < 1e-3 * (1. + numeric.abs()),
                "{smiles} axis {i}: numeric {numeric}, analytic {}",
                grad[i]
            );
        }
    }
}

#[test]
fn uff_minimization_decreases_energy_and_handles_unparameterized_atoms() {
    let mol = rdkit_added_hs_molecule(&parse("CCO").unwrap()).unwrap();
    let coords: Vec<_> = (0..mol.atom_count())
        .map(|i| [i as f64, (i as f64 * 1.7).sin(), (i as f64 * 1.3).cos()])
        .collect();
    let field = RdkitUffField::new(&mol, &coords, 10., false);
    let mut pos: Vec<_> = coords.iter().flatten().copied().collect();
    let before = field.energy(&pos);
    assert!([0, 1].contains(&field.minimize(&mut pos, 20, 1e-4)));
    assert!(field.energy(&pos) <= before);
    let unparameterized = parse("[*]").unwrap();
    assert!(!rdkit_uff_has_all_params(&unparameterized));
    let field = RdkitUffField::new(&unparameterized, &[[0.; 3]], 10., true);
    assert_eq!(field.energy(&[0.; 3]), 0.);
    assert_eq!(field.minimize(&mut [0.; 3], 20, 1e-4), 0);
}

#[test]
fn pyramidal_inversion_retains_rdkit_reference_energy_and_gradient() {
    // Independently measured using rdkit==2026.3.1, Chem.AddHs and
    // UFFGetMoleculeForceField with the coordinates below. This API promises
    // RDKit's expressions. For P/As inversion its analytic gradient differs
    // from a finite difference of its energy; do not replace the parity
    // contract with the finite-difference contract used for organic cases.
    // Upstream: Release_2026_03_1/Code/ForceField/UFF/Inversion.cpp:113.
    for (smiles, expected_energy, expected_y_gradient) in [
        ("CP(C)C", 67330.09078592298, -1224.619113639887),
        ("C[As](C)C", 66325.53503120068, -1087.4359414743303),
    ] {
        let mol = rdkit_added_hs_molecule(&parse(smiles).unwrap()).unwrap();
        let coords: Vec<_> = (0..mol.atom_count())
            .map(|i| {
                let t = i as f64 * 1.7;
                [t.cos() * 2.3, t.sin() * 2.3, i as f64 * 0.7]
            })
            .collect();
        let field = RdkitUffField::new(&mol, &coords, 10., true);
        let pos: Vec<_> = coords.iter().flatten().copied().collect();
        let mut grad = vec![0.; pos.len()];
        field.gradient(&pos, &mut grad);
        assert!(
            (field.energy(&pos) - expected_energy).abs() < 1e-7,
            "{smiles}"
        );
        assert!(
            (grad[1] - expected_y_gradient).abs() < 1e-8,
            "{smiles}: {}",
            grad[1]
        );
        assert!(
            field
                .debug_terms()
                .iter()
                .any(|(kind, _, _)| *kind == "inversion")
        );
    }
}
