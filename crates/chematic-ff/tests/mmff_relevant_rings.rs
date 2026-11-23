//! Exposed fused-ring development cases against RDKit 2026.09.2.
use chematic_ff::rdkit_mmff::RdkitMmffField;

#[test]
fn relevant_ring_order_matches_rdkit_mmff_development_cases() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("mmff_relevant_rings_fixture.json")).unwrap();
    for case in fixture["results"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let mol = chematic_chem::add_hydrogens(
            &chematic_smiles::parse(case["smiles"].as_str().unwrap()).unwrap(),
        );
        let expected_types: Vec<u8> = serde_json::from_value(case["atom_types"].clone()).unwrap();
        let expected_charges: Vec<f64> = serde_json::from_value(case["charges"].clone()).unwrap();
        let types = chematic_ff::assign_mmff94_numeric_types(&mol).unwrap();
        let charges = chematic_ff::mmff94_charges_numeric(&mol).unwrap();
        assert_eq!(types, expected_types, "{name} atom types");
        assert_eq!(charges.len(), expected_charges.len(), "{name} charge count");
        for (i, (actual, expected)) in charges.iter().zip(&expected_charges).enumerate() {
            assert!(
                (actual - expected).abs() < 1e-8,
                "{name} charge {i}: {actual} vs {expected}"
            );
        }
        let coords: Vec<[f64; 3]> = serde_json::from_value(case["coords"].clone()).unwrap();
        let pos: Vec<f64> = coords.iter().flatten().copied().collect();
        let expected_gradient: Vec<f64> = serde_json::from_value(case["gradient"].clone()).unwrap();
        assert_eq!(pos.len(), expected_gradient.len(), "{name} gradient count");
        let field = RdkitMmffField::new(&mol, &coords, 100.0, true).unwrap();
        let actual = field.energy(&pos);
        let expected = case["energy"].as_f64().unwrap();
        assert!(
            (actual - expected).abs() < 1e-6,
            "{name} energy {actual} vs {expected}"
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
