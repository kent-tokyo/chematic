use crate::rdkit_uff::{RdkitUffField, rdkit_uff_has_all_params};
#[test]
fn uff_boundary_elements_energies_and_gradients_match_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-uff-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    // The earlier RDKit release had an inversion-gradient sign error for
    // P/As/Sb/Bi. Preserve its energy fixture, but compare the affected Sb/Bi
    // gradients with the same coordinates and settings on 2026.09.2.
    let inversion: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.09.2-uff-inversion-boundary.json"
    ))
    .unwrap();
    assert_eq!(inversion["rdkit_version"], "2026.09.2");
    assert_eq!(inversion["rows"].as_array().unwrap().len(), 6);
    let mut failures = Vec::new();
    for (index, row) in fixture["rows"].as_array().unwrap().iter().enumerate() {
        let source = row["smiles"].as_str().unwrap();
        let parent = chematic_smiles::parse(source).unwrap();
        let mol = if row["coords"].as_array().unwrap().len() == 1 {
            parent
        } else {
            chematic_chem::add_hydrogens(&parent)
        };
        let coords: Vec<[f64; 3]> = serde_json::from_value(row["coords"].clone()).unwrap();
        let field = RdkitUffField::new(
            &mol,
            &coords,
            row["threshold"].as_f64().unwrap(),
            row["ignore"].as_bool().unwrap(),
        );
        let pos: Vec<f64> = coords.iter().flatten().copied().collect();
        let energy = field.energy(&pos);
        let expected = row["energy"].as_f64().unwrap();
        if (energy - expected).abs() > 1e-6 * (1. + expected.abs()) {
            failures.push(format!("{source}: energy {energy} != {expected}"));
        }
        let mut gradient = vec![0.; pos.len()];
        field.gradient(&pos, &mut gradient);
        let corrected = inversion["rows"]
            .as_array()
            .unwrap()
            .iter()
            .find(|item| item["original_row"].as_u64() == Some(index as u64));
        if let Some(corrected) = corrected {
            assert_eq!(corrected["smiles"], row["smiles"]);
        }
        let expected: Vec<f64> =
            serde_json::from_value(corrected.unwrap_or(row)["gradient"].clone()).unwrap();
        for (i, (a, b)) in gradient.iter().zip(expected.iter()).enumerate() {
            if !a.is_finite() || (a - b).abs() > 1e-5 * (1. + b.abs()) {
                failures.push(format!("{source}: gradient {i}: {a} != {b}"));
                break;
            }
        }
        if rdkit_uff_has_all_params(&mol) != row["all_params"].as_bool().unwrap() {
            failures.push(format!("{source}: parameter availability"));
        }
        assert_eq!(field.num_atoms(), coords.len());
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
