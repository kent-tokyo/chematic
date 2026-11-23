//! Fixed-coordinate RDKit 2026.09.2 UFF checks, including exposed development outliers.
use chematic_ff::rdkit_uff::RdkitUffField;

#[test]
fn non_tetrahedral_hybridization_matches_rdkit_uff() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("uff_nontetrahedral_fixture.json")).unwrap();
    for case in fixture["results"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let mol = chematic_chem::add_hydrogens(
            &chematic_smiles::parse(case["smiles"].as_str().unwrap()).unwrap(),
        );
        let coords: Vec<[f64; 3]> = serde_json::from_value(case["coords"].clone()).unwrap();
        let pos: Vec<f64> = coords.iter().flatten().copied().collect();
        let expected_gradient: Vec<f64> = serde_json::from_value(case["gradient"].clone()).unwrap();
        assert_eq!(pos.len(), expected_gradient.len(), "{name} atom count");
        let field = RdkitUffField::new(&mol, &coords, 10.0, true);
        let energy = field.energy(&pos);
        let expected = case["energy"].as_f64().unwrap();
        assert!(
            (energy - expected).abs() < 1e-6,
            "{name} energy {energy} vs {expected}"
        );
        let mut gradient = vec![0.0; pos.len()];
        field.gradient(&pos, &mut gradient);
        for (i, (actual, expected)) in gradient.iter().zip(&expected_gradient).enumerate() {
            assert!(
                (actual - expected).abs() < 1e-5,
                "{name} gradient {i}: {actual} vs {expected}"
            );
        }
    }
}
