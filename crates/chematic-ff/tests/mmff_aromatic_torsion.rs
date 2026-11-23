//! Exposed development case; RDKit 2026.09.2 on fixed coordinates.
use chematic_ff::rdkit_mmff::RdkitMmffField;

#[test]
fn mmff_aromatic_type_set_controls_empirical_torsion() {
    let case: serde_json::Value =
        serde_json::from_str(include_str!("mmff_aromatic_torsion_fixture.json")).unwrap();
    let mol = chematic_chem::add_hydrogens(
        &chematic_smiles::parse(case["smiles"].as_str().unwrap()).unwrap(),
    );
    let coords: Vec<[f64; 3]> = serde_json::from_value(case["coords"].clone()).unwrap();
    let pos: Vec<f64> = coords.iter().flatten().copied().collect();
    let expected_gradient: Vec<f64> = serde_json::from_value(case["gradient"].clone()).unwrap();
    let field = RdkitMmffField::new(&mol, &coords, 100.0, true).unwrap();
    let actual = field.energy(&pos);
    let expected = case["energy"].as_f64().unwrap();
    assert!(
        (actual - expected).abs() < 1e-6,
        "energy {actual} vs {expected}"
    );
    let mut gradient = vec![0.0; pos.len()];
    field.gradient(&pos, &mut gradient);
    assert_eq!(gradient.len(), expected_gradient.len());
    for (i, (actual, expected)) in gradient.iter().zip(&expected_gradient).enumerate() {
        assert!(
            (actual - expected).abs() < 1e-5,
            "gradient {i}: {actual} vs {expected}"
        );
    }
}
