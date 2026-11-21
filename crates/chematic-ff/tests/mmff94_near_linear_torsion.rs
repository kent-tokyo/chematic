use chematic_ff::rdkit_mmff::{MmffVariant, RdkitMmffField};

#[test]
fn torsion_gradient_tracks_energy_at_near_linear_angles() {
    let mol = chematic_smiles::parse("CCCC").unwrap();
    let coords = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [3.0, 0.000002, 0.0],
        [4.0, 0.4, 0.8],
    ];
    let field = RdkitMmffField::with_terms(
        &mol,
        &coords,
        MmffVariant::Mmff94,
        100.0,
        true,
        [false, false, false, false, true, false, false],
    )
    .unwrap();
    let mut position = coords.as_flattened().to_vec();
    let mut grad = vec![0.0; position.len()];
    field.gradient_with_consistent_torsions(&position, &mut grad);
    for axis in 0..position.len() {
        let step = 1e-10;
        position[axis] += step;
        let plus = field.energy(&position);
        position[axis] -= 2.0 * step;
        let minus = field.energy(&position);
        position[axis] += step;
        let numeric = (plus - minus) / (2.0 * step);
        assert!(grad[axis].is_finite() && numeric.is_finite());
        assert!(
            (grad[axis] - numeric).abs() <= 0.01 * (1.0 + numeric.abs()),
            "axis {axis}: analytic {} numeric {numeric}",
            grad[axis]
        );
    }
}

#[test]
fn exactly_linear_torsion_has_zero_gradient_inside_energy_cutoff() {
    let mol = chematic_smiles::parse("CCCC").unwrap();
    let coords = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [4.0, 0.4, 0.8],
    ];
    let field = RdkitMmffField::with_terms(
        &mol,
        &coords,
        MmffVariant::Mmff94,
        100.0,
        true,
        [false, false, false, false, true, false, false],
    )
    .unwrap();
    let mut grad = vec![0.0; coords.len() * 3];
    field.gradient_with_consistent_torsions(coords.as_flattened(), &mut grad);
    assert!(grad.iter().all(|value| *value == 0.0));
}

#[test]
fn consistent_and_compatible_torsions_agree_away_from_degeneracy() {
    let mol = chematic_smiles::parse("CCCC").unwrap();
    let coords = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.5, 1.0, 0.0],
        [3.0, 1.5, 1.2],
    ];
    let field = RdkitMmffField::with_terms(
        &mol,
        &coords,
        MmffVariant::Mmff94,
        100.0,
        true,
        [false, false, false, false, true, false, false],
    )
    .unwrap();
    let mut compatible = vec![0.0; coords.len() * 3];
    let mut consistent = compatible.clone();
    field.gradient(coords.as_flattened(), &mut compatible);
    field.gradient_with_consistent_torsions(coords.as_flattened(), &mut consistent);
    for (a, b) in compatible.iter().zip(consistent) {
        assert!((a - b).abs() < 1e-10);
    }
}

#[test]
fn near_linear_amine_plane_terms_match_finite_differences() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("mmff94_near_linear_amine.json")).unwrap();
    let smiles = fixture["smiles"].as_str().unwrap();
    let coords: Vec<[f64; 3]> = serde_json::from_value(fixture["coords"].clone()).unwrap();
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
    assert_eq!(coords.len(), mol.atom_count());
    for (terms, step) in [
        ([false, false, false, false, true, false, false], 1e-10),
        ([false, false, false, true, false, false, false], 1e-8),
    ] {
        let field =
            RdkitMmffField::with_terms(&mol, &coords, MmffVariant::Mmff94, 100.0, true, terms)
                .unwrap();
        let mut flat = coords.as_flattened().to_vec();
        let mut gradient = vec![0.0; flat.len()];
        field.gradient_with_consistent_planes(&flat, &mut gradient);
        for axis in 0..flat.len() {
            flat[axis] += step;
            let plus = field.energy(&flat);
            flat[axis] -= 2.0 * step;
            let minus = field.energy(&flat);
            flat[axis] += step;
            let numeric = (plus - minus) / (2.0 * step);
            assert!(gradient[axis].is_finite() && numeric.is_finite());
            assert!(
                (gradient[axis] - numeric).abs() <= 0.02 * (1.0 + numeric.abs()),
                "axis {axis}: analytic {} numeric {numeric}",
                gradient[axis]
            );
        }
    }
}

#[test]
fn ordinary_out_of_plane_gradient_matches_energy_and_compatible_gradient() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse("c1ccccc1").unwrap());
    let mut coords = Vec::new();
    for radius in [1.4, 2.4] {
        for i in 0..6 {
            let angle = i as f64 * std::f64::consts::PI / 3.0;
            coords.push([
                radius * angle.cos(),
                radius * angle.sin(),
                0.05 * (2.0 * angle).sin(),
            ]);
        }
    }
    let field = RdkitMmffField::with_terms(
        &mol,
        &coords,
        MmffVariant::Mmff94,
        100.0,
        true,
        [false, false, false, true, false, false, false],
    )
    .unwrap();
    let mut flat = coords.as_flattened().to_vec();
    let mut gradient = vec![0.0; flat.len()];
    let mut compatible = gradient.clone();
    field.gradient_with_consistent_planes(&flat, &mut gradient);
    field.gradient(&flat, &mut compatible);
    for axis in 0..flat.len() {
        let step = 1e-6;
        flat[axis] += step;
        let plus = field.energy(&flat);
        flat[axis] -= 2.0 * step;
        let minus = field.energy(&flat);
        flat[axis] += step;
        let numeric = (plus - minus) / (2.0 * step);
        assert!((gradient[axis] - numeric).abs() < 1e-5 * (1.0 + numeric.abs()));
        assert!((gradient[axis] - compatible[axis]).abs() < 1e-10);
    }
}
