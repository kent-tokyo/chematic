use super::*;

#[test]
fn public_default_field_matches_native_ionic_pairs_and_collision_guards() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../../validation/rdkit-2026.03.1-mmff-nonbonded-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let mol = chematic_smiles::parse(source).unwrap();
        let types = crate::mmff94_numeric::assign_mmff94_numeric_types(&mol).unwrap();
        let charges = crate::mmff94_numeric::mmff94_charges_numeric(&mol).unwrap();
        for i in 0..2 {
            assert_eq!(
                u64::from(types[i]),
                row["atom_types"][i].as_u64().unwrap(),
                "{source}"
            );
            assert!(
                (charges[i] - row["charges"][i].as_f64().unwrap()).abs() < 1e-12,
                "{source}"
            );
        }
        let distance = row["distance"].as_f64().unwrap();
        let ignore = row["ignore_interfrag"].as_bool().unwrap();
        let coords = [[0., 0., 0.], [distance, 0., 0.]];
        let field = RdkitMmffField::new(&mol, &coords, 100., ignore).unwrap();
        assert_eq!(field.num_atoms(), 2);
        let pos = [0., 0., 0., distance, 0., 0.];
        let energy = field.energy(&pos);
        let expected = row["energy"].as_f64().unwrap();
        assert!(
            energy.is_finite() && (energy - expected).abs() < 1e-8 + expected.abs() * 1e-10,
            "{source},d={distance},ignore={ignore}: {energy} != {expected}"
        );
        let mut gradient = [0.; 6];
        field.gradient(&pos, &mut gradient);
        for (axis, value) in row["gradient"].as_array().unwrap().iter().enumerate() {
            let expected = value.as_f64().unwrap();
            assert!(
                gradient[axis].is_finite()
                    && (gradient[axis] - expected).abs() < 1e-9 + expected.abs() * 1e-10,
                "{source},d={distance},ignore={ignore},axis={axis}: {} != {expected}",
                gradient[axis]
            );
        }
        for axis in 0..3 {
            assert!((gradient[axis] + gradient[axis + 3]).abs() < 1e-12);
        }
        if ignore {
            assert_eq!(energy, 0.);
            assert_eq!(gradient, [0.; 6]);
        }
        // At exactly coincident centres RDKit returns a deterministic
        // guard vector, not the derivative of a differentiable geometry.
        if distance >= 1. && !ignore {
            for axis in 0..6 {
                let mut plus = pos;
                let mut minus = pos;
                plus[axis] += 1e-5;
                minus[axis] -= 1e-5;
                let finite_difference = (field.energy(&plus) - field.energy(&minus)) / (2e-5);
                assert!(
                    (gradient[axis] - finite_difference).abs() < 1e-4 + gradient[axis].abs() * 2e-6,
                    "{source},d={distance},axis={axis}"
                );
            }
        }
    }
    assert_eq!(fixture["rows"].as_array().unwrap().len(), 104);
}
