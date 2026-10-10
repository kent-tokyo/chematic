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

#[test]
fn heteroatom_force_field_gradients_match_numerical_energy_derivatives() {
    // These less common central bonds exercise empirical MMFF torsion rules.
    // Numerical differentiation checks the independent energy/force contract.
    for source in [
        "COOC",
        "CSSC",
        "COSC",
        "CSOC",
        "C[SiH2][SiH2]C",
        "CPPC",
        "CSPC",
        "CPSC",
        "CNP(C)C",
        "C=NN=C",
        "CC=NNC",
        "CN=NC",
        "C=CNC",
        "C=CPC",
        "C=CSC",
        "C=COC",
        "CC=CC",
        "C[NH+]=NC",
        "C=NSC",
        "C=NPC",
        "CS(=O)SC",
        "CP(=O)(C)PC",
        "C[SiH2]OC",
        "C[SiH2]NC",
    ] {
        let mol = chematic_chem::add_hydrogens(&parse(source).unwrap());
        let model = Mmff94EnergyModel::new(&mol).unwrap();
        let xyz = geometry(mol.atom_count());
        let gradient = model.bounded_analytic_gradient(&xyz);
        assert!(model.energy(&xyz).is_finite(), "{source}");
        let mut work = xyz.clone();
        let h = 1e-5;
        for i in 0..xyz.len() {
            for axis in 0..3 {
                work[i][axis] = xyz[i][axis] + h;
                let plus = model.energy(&work);
                work[i][axis] = xyz[i][axis] - h;
                let minus = model.energy(&work);
                work[i][axis] = xyz[i][axis];
                let expected = (plus - minus) / (2.0 * h);
                assert!(
                    (gradient[i][axis] - expected).abs() < 2e-4 * (1.0 + expected.abs()),
                    "{source} atom {i} axis {axis}: {} != {expected}",
                    gradient[i][axis]
                );
            }
        }
        for axis in 0..3 {
            assert!(
                gradient.iter().map(|g| g[axis]).sum::<f64>().abs() < 1e-6,
                "net internal force: {source}"
            );
        }
    }
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
    for text in ["[Na]", "[Mg]", "[Zn]", "[Fe+]", "[Cu+3]", "C[Na]", "[U]"] {
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

#[test]
fn legacy_bond_charge_increments_are_antisymmetric_for_every_public_type() {
    use crate::mmff94::MMFF94Type::*;
    let types = [
        C_sp3,
        C_sp2_Alkene,
        C_sp_Alkyne,
        C_Aromatic,
        C_Carbonyl,
        C_Carboxylic,
        C_Carbamate,
        C_Ester,
        C_Amide,
        C_Imide,
        C_CarbamideN,
        N_sp3_Amine,
        N_sp3_AmineAromatic,
        N_sp2_Imine,
        N_sp2_Aromatic,
        N_sp2_Carbonyl,
        N_sp_Nitrile,
        N_Amide,
        N_Carbamate,
        N_Ester,
        N_Imide,
        N_Aromatic_5ring,
        N_Aromatic_6ring,
        N_Aromatic_Pyridine,
        N_Aromatic_Pyrrole,
        N_Aromatic_Imidazole,
        N_Aromatic_Triazole,
        N_Aromatic_Tetrazole,
        N_Aromatic_Pyrimidine,
        N_Aromatic_Pyrazine,
        O_Alcohol,
        O_Phenol,
        O_Ether,
        O_Carbonyl,
        O_Carboxylic,
        O_Carbamate,
        O_Ester,
        O_Amide,
        O_Imide,
        O_CarbamideN,
        O_Sulfoxide,
        O_Sulfone,
        S_Thiol,
        S_Thioether,
        S_Disulfide,
        S_Sulfoxide,
        S_Sulfone,
        S_Aromatic,
        P_sp3,
        P_Oxide,
        Si_sp3,
        Si_sp2,
        F,
        Cl,
        Br,
        I,
        H_Carbon,
        H_Nitrogen,
        H_Oxygen,
        H_Sulfur,
        H_Halogen,
        H_Aromatic,
        Generic,
    ];
    for a in types {
        for b in types {
            let forward = crate::mmff94_bci::bci(a, b);
            let reverse = crate::mmff94_bci::bci(b, a);
            assert!(forward.is_finite() && reverse.is_finite());
            assert!(
                (forward + reverse).abs() < 1e-12,
                "{a:?} -> {b:?}: {forward} vs {reverse}"
            );
        }
    }
}

#[test]
fn legacy_charge_assignment_preserves_formal_charge_across_functional_groups() {
    for source in [
        "C",
        "C=C",
        "C#C",
        "c1ccccc1",
        "CN",
        "CNC",
        "CC(=O)N",
        "CC(=O)NC(=O)C",
        "CO",
        "COC",
        "CC=O",
        "CC(=O)O",
        "CC(=O)OC",
        "CC(=O)[O-]",
        "CS",
        "CSC",
        "CSSC",
        "CS(=O)C",
        "CS(=O)(=O)C",
        "c1ccsc1",
        "CP",
        "CP(=O)(O)O",
        "C[SiH3]",
        "C[SiH]=C",
        "CF",
        "CCl",
        "CBr",
        "CI",
        "[NH4+]",
        "[H]F",
        "[H]Cl",
        "[H]Br",
    ] {
        let mol = chematic_chem::add_hydrogens(&parse(source).unwrap());
        let expected = mol.atoms().map(|(_, a)| f64::from(a.charge)).sum::<f64>();
        let charges = crate::mmff94_bci::mmff94_charges_bci(&mol).unwrap();
        assert_eq!(charges.len(), mol.atom_count());
        assert!(charges.iter().all(|x| x.is_finite()));
        near(charges.iter().sum(), expected);
    }
    let mol = parse("[Xe]").unwrap();
    let e = crate::mmff94_bci::mmff94_charges_bci(&mol).unwrap_err();
    assert!(matches!(
        e,
        crate::mmff94::AssignError::UnsupportedElement(_)
    ));
    assert!(e.to_string().contains("54"));
}

#[test]
fn dreiding_parameter_families_have_finite_physical_scales() {
    use crate::dreiding::DREIDINGType::*;
    use crate::params::*;
    use chematic_core::BondOrder;
    let types = [
        H_, C_3, C_2, C_1, C_R, N_3, N_2, N_1, N_R, O_3, O_2, O_R, S_3, S_R, P_3, F_, Cl, Br, I_,
        X_,
    ];
    for &a in &types {
        let (radius, depth) = dreiding_vdw(a);
        assert!((2.0..5.0).contains(&radius));
        assert!(depth.is_finite() && depth > 0.0);
        assert!((0.0..=std::f64::consts::PI).contains(&dreiding_angle(a)));
        assert!(!a.symbol().is_empty());
        for &b in &types {
            for order in [
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Triple,
                BondOrder::Quadruple,
                BondOrder::Aromatic,
                BondOrder::Up,
                BondOrder::Down,
                BondOrder::QueryAny,
            ] {
                let distance = dreiding_bond_len(a, b, order);
                assert!(distance.is_finite() && (0.5..3.0).contains(&distance));
                assert!(dreiding_torsion_barrier(a, b) >= 0.0);
            }
        }
    }
    assert!(
        dreiding_bond_len(C_3, C_3, BondOrder::Single)
            > dreiding_bond_len(C_2, C_2, BondOrder::Double)
    );
    assert!(
        dreiding_bond_len(C_2, C_2, BondOrder::Double)
            > dreiding_bond_len(C_1, C_1, BondOrder::Triple)
    );
}

#[test]
fn uff_parameter_families_have_positive_radii_and_well_depths() {
    use crate::uff::UffType::*;
    for atom_type in [
        C_3, C_2, C_1, C_R, N_3, N_2, N_1, N_R, O_3, O_2, O_1, O_R, S_3, S_2, S_R, P_3, P_R, H_,
        F_, Cl, Br, I_, Li, Na, K, Ca, Mg, Fe, Co, Ni, Cu, Zn, Mn, Cr, V_, Mo, W_, Pd, Pt, Au, Ag,
        Hg, Al, Si, Unknown,
    ] {
        assert!((0.2..2.5).contains(&atom_type.r1()));
        assert!((0.0..=180.0).contains(&atom_type.theta0()));
        assert!((2.0..5.0).contains(&atom_type.x1()));
        assert!(atom_type.d1().is_finite() && atom_type.d1() > 0.0);
    }
    assert!(C_3.r1() > C_2.r1() && C_2.r1() > C_1.r1());
    assert!(F_.r1() < Cl.r1() && Cl.r1() < Br.r1() && Br.r1() < I_.r1());
    assert_eq!(C_1.theta0(), 180.0);
    assert_eq!(C_2.theta0(), 120.0);
}
