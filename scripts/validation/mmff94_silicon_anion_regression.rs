use chematic_ff::rdkit_mmff::{MmffVariant, RdkitMmffField};
use chematic_ff::{assign_mmff94_numeric_types, mmff94_charges_numeric};

#[test]
fn aromatic_silicon_anions_match_rdkit_force_field_oracle() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("mmff94_silicon_anion_fixture.json")).unwrap();
    let molecules = fixture["molecules"].as_array().unwrap();
    assert_eq!(molecules.len(), 3);
    for row in molecules {
        let name = row["name"].as_str().unwrap();
        let smiles = row["smiles"].as_str().unwrap();
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let charges = mmff94_charges_numeric(&mol).unwrap();
        let expected_types: Vec<u8> = serde_json::from_value(row["types"].clone()).unwrap();
        let expected_charges: Vec<f64> = serde_json::from_value(row["charges"].clone()).unwrap();
        assert_eq!(types, expected_types, "{name}: types");
        for (i, (actual, expected)) in charges.iter().zip(expected_charges).enumerate() {
            assert!(
                (actual - expected).abs() < 1e-9,
                "{name} atom {i}: {actual} != {expected}"
            );
        }
        let coords: Vec<[f64; 3]> = serde_json::from_value(row["coords"].clone()).unwrap();
        for (label, variant) in [
            ("MMFF94", MmffVariant::Mmff94),
            ("MMFF94s", MmffVariant::Mmff94s),
        ] {
            let field =
                RdkitMmffField::with_terms(&mol, &coords, variant, 100.0, true, [true; 7]).unwrap();
            let reference = row["fields"][label]["energy"].as_f64().unwrap();
            let expected: Vec<f64> =
                serde_json::from_value(row["fields"][label]["gradient"].clone()).unwrap();
            let energy = field.energy(coords.as_flattened());
            assert!(
                (energy - reference).abs() < 1e-6,
                "{name} {label}: {energy} != {reference}"
            );
            let mut gradient = vec![0.0; coords.len() * 3];
            field.gradient(coords.as_flattened(), &mut gradient);
            assert_eq!(gradient.len(), expected.len());
            for (i, (actual, expected)) in gradient.iter().zip(expected).enumerate() {
                assert!(
                    (actual - expected).abs() < 1e-6,
                    "{name} {label} axis {i}: {actual} != {expected}"
                );
            }
        }
    }
}
