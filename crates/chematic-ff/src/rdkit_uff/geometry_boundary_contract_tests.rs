use super::*;

#[test]
fn bridging_hydrogen_matches_rdkit_energy_and_gradient() {
    for json in [
        include_str!("../../../../validation/rdkit-2026.09.2-uff-bridging-hydrogen.json"),
        include_str!("../../../../validation/rdkit-2026.09.2-uff-bridging-hydrogen-exposed.json"),
    ] {
        let fixture: serde_json::Value = serde_json::from_str(json).unwrap();
        let mol = chematic_chem::add_hydrogens(
            &chematic_smiles::parse(fixture["smiles"].as_str().unwrap()).unwrap(),
        );
        let coords: Vec<[f64; 3]> = fixture["coords"]
            .as_array()
            .unwrap()
            .iter()
            .map(|point| {
                [
                    point[0].as_f64().unwrap(),
                    point[1].as_f64().unwrap(),
                    point[2].as_f64().unwrap(),
                ]
            })
            .collect();
        assert_eq!(mol.atom_count(), coords.len());
        let field = RdkitUffField::new(&mol, &coords, 10.0, true);
        let pos = coords.as_flattened();
        assert!((field.energy(pos) - fixture["energy"].as_f64().unwrap()).abs() < 1e-8);
        let mut gradient = vec![0.0; pos.len()];
        field.gradient(pos, &mut gradient);
        for (actual, expected) in gradient.iter().zip(fixture["gradient"].as_array().unwrap()) {
            assert!((actual - expected.as_f64().unwrap()).abs() < 1e-7);
        }
    }
}

fn matches_native(actual: f64, expected: &serde_json::Value) -> bool {
    match expected.as_str() {
        Some("nan") => actual.is_nan(),
        Some("+inf") => actual == f64::INFINITY,
        Some("-inf") => actual == f64::NEG_INFINITY,
        None => {
            let expected = expected.as_f64().unwrap();
            actual.is_finite() && (actual - expected).abs() < 1e-7 + expected.abs() * 1e-9
        }
        other => panic!("unrecognized reference {other:?}"),
    }
}

#[test]
fn small_ring_parameters_and_degenerate_coordinate_guards_match_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../../validation/rdkit-2026.03.1-uff-geometries-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let mode = row["mode"].as_str().unwrap();
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(source).unwrap());
        let coords: Vec<[f64; 3]> = row["coords"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| {
                [
                    v[0].as_f64().unwrap(),
                    v[1].as_f64().unwrap(),
                    v[2].as_f64().unwrap(),
                ]
            })
            .collect();
        assert_eq!(coords.len(), mol.atom_count());
        let field = RdkitUffField::new(&mol, &coords, 100., false);
        assert_eq!(field.num_atoms(), mol.atom_count());
        let pos = coords.as_flattened().to_vec();
        let energy = field.energy(&pos);
        assert!(
            matches_native(energy, &row["energy"]),
            "{source},{mode}: E={energy} != {}",
            row["energy"]
        );
        let mut gradient = vec![0.; pos.len()];
        field.gradient(&pos, &mut gradient);
        for (axis, expected) in row["gradient"].as_array().unwrap().iter().enumerate() {
            assert!(
                matches_native(gradient[axis], expected),
                "{source},{mode},axis={axis}: {} != {expected}",
                gradient[axis]
            );
        }
        // Degenerate inputs compare the native guard/nonfinite behavior;
        // they do not represent usable or optimized molecular geometries.
        if mode == "distorted" {
            assert!(energy.is_finite() && gradient.iter().all(|v| v.is_finite()));
            for axis in 0..pos.len() {
                let mut plus = pos.clone();
                let mut minus = pos.clone();
                plus[axis] += 1e-6;
                minus[axis] -= 1e-6;
                let numerical = (field.energy(&plus) - field.energy(&minus)) / 2e-6;
                assert!(
                    (gradient[axis] - numerical).abs() < 1e-3 + gradient[axis].abs() * 2e-5,
                    "{source},axis={axis}: {} != {numerical}",
                    gradient[axis]
                );
            }
        }
    }
    assert_eq!(fixture["rows"].as_array().unwrap().len(), 45);
}
