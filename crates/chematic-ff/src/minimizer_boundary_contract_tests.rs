use crate::*;
use chematic_smiles::parse;

fn geometry(n: usize) -> Vec<[f64; 3]> {
    (0..n)
        .map(|i| {
            [
                i as f64 * 1.25,
                if i % 2 == 0 { 0. } else { 0.8 },
                (i as f64 * 0.7).sin() * 0.35,
            ]
        })
        .collect()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 1e-8 * (1. + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

#[test]
fn torsion_scan_is_rigid_transform_invariant_and_keeps_input_coordinates() {
    for smiles in ["CCCC", "CCCCC", "CCCO", "CCCN"] {
        let mol = parse(smiles).unwrap();
        let xyz = geometry(mol.atom_count());
        let original = xyz.clone();
        let scan = mmff94_torsion_scan(&mol, &xyz, 0, 1, 2, 3, 12).unwrap();
        assert_eq!(scan.len(), 12);
        near(scan[0].1, mmff94_total_energy(&mol, &xyz).unwrap());
        let moved = xyz
            .iter()
            .map(|p| [-p[1] + 7., p[0] - 3., p[2] + 5.])
            .collect::<Vec<_>>();
        let translated = mmff94_torsion_scan(&mol, &moved, 0, 1, 2, 3, 12).unwrap();
        for (i, ((angle, energy), (_, reference))) in scan.iter().zip(translated.iter()).enumerate()
        {
            near(*angle, i as f64 * 30.);
            near(*energy, *reference);
        }
        let min = scan.iter().map(|x| x.1).fold(f64::INFINITY, f64::min);
        let max = scan.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max - min > 1e-3,
            "{smiles}: scan must change internal geometry"
        );
        assert_eq!(xyz, original);
        for requested in [0, 1, 2] {
            assert_eq!(
                mmff94_torsion_scan(&mol, &xyz, 0, 1, 2, 3, requested)
                    .unwrap()
                    .len(),
                2
            );
        }
    }
}

#[test]
fn minimizers_respect_iteration_limits_and_prepared_energy_decreases() {
    for smiles in ["CCCC", "CCCCCCCCCCCCCCCCCC"] {
        let mol = parse(smiles).unwrap();
        let xyz = geometry(mol.atom_count());
        let model = Mmff94EnergyModel::new(&mol).unwrap();
        let initial = model.energy(&xyz);
        let mut zero = xyz.clone();
        let result = model.minimize_lbfgs(&mut zero, 0).unwrap();
        assert_eq!(zero, xyz);
        assert_eq!(result.iterations, 0);
        assert_eq!(result.termination, Mmff94TerminationReason::IterationLimit);
        near(result.energy, initial);
        let mut numeric = xyz.clone();
        let result = model.minimize_lbfgs(&mut numeric, 3).unwrap();
        assert!(result.iterations <= 3);
        assert!(result.energy < initial);
        near(result.energy, model.energy(&numeric));
        let mut full = xyz.clone();
        let result = minimize_mmff94_full(&mol, &mut full, 3).unwrap();
        assert!(result.iterations <= 3);
        assert!(result.energy.is_finite());
        near(result.energy, mmff94_total_energy(&mol, &full).unwrap());
    }
}

#[test]
fn cutoff_minimizer_reduces_its_consistent_objective() {
    let mol = parse("CCCCCC").unwrap();
    let mut xyz = geometry(mol.atom_count());
    let model = Mmff94EnergyModel::new(&mol).unwrap();
    let initial = model.energy(&xyz);
    let result = model
        .minimize_lbfgs_bounded_analytic_cutoff(&mut xyz, 5)
        .unwrap();
    assert!(result.energy < initial);
    near(result.energy, model.energy(&xyz));
    let components = [
        model.bond_angle_gradient(&xyz),
        model.stretch_bend_gradient(&xyz),
        model.torsion_gradient(&xyz),
        model.oop_gradient(&xyz),
        model.cutoff_nonbonded_gradient(&xyz),
    ];
    let total = model.bounded_analytic_gradient_cutoff(&xyz);
    for i in 0..xyz.len() {
        for axis in 0..3 {
            near(total[i][axis], components.iter().map(|g| g[i][axis]).sum());
        }
    }
}

#[test]
fn mmff_boundary_types_and_charges_match_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-mmff-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let text = row["smiles"].as_str().unwrap();
        let mol = chematic_chem::add_hydrogens(&parse(text).unwrap());
        let expected: Vec<u8> = row["types"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as u8)
            .collect();
        let actual = assign_mmff94_numeric_types(&mol);
        if actual.as_ref().ok() != Some(&expected) {
            failures.push(format!("types {text}: {actual:?} != {expected:?}"));
        }
        let expected: Vec<f64> = row["charges"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();
        match mmff94_charges_numeric(&mol) {
            Ok(actual)
                if actual.len() == expected.len()
                    && actual
                        .iter()
                        .zip(&expected)
                        .all(|(a, e)| (a - e).abs() < 1e-8) => {}
            other => failures.push(format!("charges {text}: {other:?} != {expected:?}")),
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn mmff_unsupported_metal_typing_returns_an_actionable_error() {
    for text in ["[Na+]", "[Mg+2]", "[Zn+2]", "[U]"] {
        let mol = parse(text).unwrap();
        let error = assign_mmff94_numeric_types(&mol).unwrap_err();
        assert!(error.to_string().contains("unsupported element"));
        let error = Mmff94EnergyModel::new(&mol).err().unwrap();
        assert!(error.to_string().contains("type assignment failed"));
    }
}

#[test]
fn empty_minimizers_converge_without_allocating_a_search() {
    let mol = chematic_core::MoleculeBuilder::new().build();
    let model = Mmff94EnergyModel::new(&mol).unwrap();
    for result in [
        minimize_mmff94_full(&mol, &mut [], 10).unwrap(),
        model.minimize_lbfgs(&mut [], 10).unwrap(),
        model
            .minimize_bfgs(&mut [], 10, Mmff94Convergence::default())
            .unwrap(),
    ] {
        assert_eq!(result.iterations, 0);
        assert!(result.converged);
        assert_eq!(result.energy, 0.);
        assert_eq!(
            result.termination,
            Mmff94TerminationReason::GradientConverged
        );
    }
}

#[test]
fn cutoff_neighbor_list_energy_gradient_and_large_bfgs_fallback_are_consistent() {
    let mol = parse(&"CO".repeat(40)).unwrap();
    let model = Mmff94EnergyModel::new(&mol).unwrap();
    let coords = geometry(mol.atom_count());
    let gradients = model.cutoff_nonbonded_gradient(&coords);
    let step = 1e-5;
    for atom in [0, 19, 40, 79] {
        for axis in 0..3 {
            let mut plus = coords.clone();
            let mut minus = coords.clone();
            plus[atom][axis] += step;
            minus[atom][axis] -= step;
            let (v1, e1) = model.cutoff_nonbonded_energy(&plus);
            let (v2, e2) = model.cutoff_nonbonded_energy(&minus);
            let reference = (v1 + e1 - v2 - e2) / (2. * step);
            assert!((reference - gradients[atom][axis]).abs() < 1e-5 * (1. + reference.abs()));
        }
    }
    let mut work = coords.clone();
    let initial = model.energy(&coords);
    let result = model
        .minimize_lbfgs_bounded_analytic_cutoff(&mut work, 2)
        .unwrap();
    assert!(result.energy.is_finite());
    assert!(result.iterations <= 2);
    assert!(initial.is_finite());
    let mol = parse(&"C".repeat(crate::mmff94_minimizer::BFGS_DENSE_MAX_ATOMS + 1)).unwrap();
    let model = Mmff94EnergyModel::new(&mol).unwrap();
    let coords = geometry(mol.atom_count());
    let mut bfgs = coords.clone();
    let mut lbfgs = coords.clone();
    let a = model
        .minimize_bfgs(&mut bfgs, 0, Mmff94Convergence::default())
        .unwrap();
    let b = model
        .minimize_lbfgs_bounded_analytic(&mut lbfgs, 0)
        .unwrap();
    assert_eq!(bfgs, coords);
    assert_eq!(lbfgs, coords);
    near(a.energy, b.energy);
    assert_eq!(a.termination, b.termination);
}
