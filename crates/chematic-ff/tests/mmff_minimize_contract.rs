use chematic_ff::rdkit_mmff::{MmffVariant, RdkitMmffField};

fn model(terms: [bool; 7]) -> (RdkitMmffField, Vec<f64>) {
    let data: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/openbabel_issues/mmff-dielectric/cases.json"
    ))
    .unwrap();
    let row = &data["cases"][0];
    let mol = chematic_chem::add_hydrogens(
        &chematic_smiles::parse(row["smiles"].as_str().unwrap()).unwrap(),
    );
    let coords: Vec<[f64; 3]> = serde_json::from_value(row["coords"].clone()).unwrap();
    let field =
        RdkitMmffField::with_terms(&mol, &coords, MmffVariant::Mmff94, 100., false, terms).unwrap();
    (field, coords.into_iter().flatten().collect())
}

#[test]
fn minimize_agrees_with_detailed_optimization_and_reduces_energy() {
    let (field, mut direct) = model([true; 7]);
    let mut detailed = direct.clone();
    let initial = field.energy(&direct);
    let status = field.minimize(&mut direct, 200, 1e-4);
    let result = field.optimize_detailed(&mut detailed, 200);
    assert_eq!(status, result.status);
    assert_eq!(direct, detailed);
    assert!(result.energy < initial);
    assert_eq!(field.energy(&direct), result.energy);
}

#[test]
fn minimize_with_no_terms_leaves_coordinates_unchanged() {
    let (field, mut coords) = model([false; 7]);
    let original = coords.clone();
    assert_eq!(field.minimize(&mut coords, 100, 1e-4), 0);
    assert_eq!(coords, original);
    assert_eq!(field.energy(&coords), 0.);
    let detailed = field.optimize_detailed(&mut coords, 100);
    assert_eq!(detailed.status, 0);
    assert_eq!(detailed.energy, 0.);
    assert_eq!(detailed.iterations, 0);
    assert_eq!(detailed.final_gradient_test, Some(0.));
    assert_eq!(coords, original);
    assert_eq!(field.optimize(&mut coords, 100), (0, 0.));
    assert_eq!(coords, original);
}
