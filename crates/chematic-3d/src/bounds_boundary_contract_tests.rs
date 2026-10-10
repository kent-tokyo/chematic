use crate::rdkit_embed::rdkit_bounds_matrix;
use serde_json::Value;
#[test]
fn distance_bound_flags_and_boundary_molecules_match_pinned_rdkit() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-bounds-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let text = row["smiles"].as_str().unwrap();
        let mol = chematic_smiles::parse(text).unwrap();
        let mol = if row["add_hs"].as_bool().unwrap() {
            chematic_chem::add_hydrogens(&mol)
        } else {
            mol
        };
        let actual = rdkit_bounds_matrix(
            &mol,
            row["set15"].as_bool().unwrap(),
            row["smooth"].as_bool().unwrap(),
            row["macrocycle14"].as_bool().unwrap(),
        )
        .unwrap();
        let expected = row["bounds"].as_array().unwrap();
        assert_eq!(actual.len(), expected.len());
        let mut mismatch = None;
        for (i, line) in actual.iter().enumerate() {
            assert_eq!(line.len(), expected[i].as_array().unwrap().len());
            for (j, value) in line.iter().enumerate() {
                let expected = expected[i][j].as_f64().unwrap();
                if !value.is_finite() || (value - expected).abs() > 1e-8 {
                    mismatch = Some((i, j, *value, expected));
                    break;
                }
            }
            if mismatch.is_some() {
                break;
            }
        }
        if let Some(m) = mismatch {
            failures.push(format!(
                "{text}, H={},set15={},smooth={},macro14={}: {m:?}",
                row["add_hs"], row["set15"], row["smooth"], row["macrocycle14"]
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
