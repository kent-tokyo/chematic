//! Exposed development outlier with sp carbons in a six-membered ring.
use chematic_ff::rdkit_mmff::RdkitMmffField;

#[test]
fn cumulene_ring_is_not_mmff_aromatic() {
    let case: serde_json::Value =
        serde_json::from_str(include_str!("mmff_sp_ring_fixture.json")).unwrap();
    let mol = chematic_chem::add_hydrogens(
        &chematic_smiles::parse(case["smiles"].as_str().unwrap()).unwrap(),
    );
    let expected_types: Vec<u8> = serde_json::from_value(case["atom_types"].clone()).unwrap();
    let expected_charges: Vec<f64> = serde_json::from_value(case["charges"].clone()).unwrap();
    let types = chematic_ff::assign_mmff94_numeric_types(&mol).unwrap();
    let charges = chematic_ff::mmff94_charges_numeric(&mol).unwrap();
    assert_eq!(types, expected_types);
    for (index, (actual, expected)) in charges.iter().zip(&expected_charges).enumerate() {
        assert!(
            (actual - expected).abs() < 1e-8,
            "charge {index}: {actual} vs {expected}"
        );
    }
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
    for (index, (actual, expected)) in gradient.iter().zip(&expected_gradient).enumerate() {
        assert!(
            (actual - expected).abs() < 1e-5,
            "gradient {index}: {actual} vs {expected}"
        );
    }
}
