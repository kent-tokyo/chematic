//! Same-coordinate USR comparison with RDKit 2026.09.1.

use chematic_3d::{RdkitUsrError, rdkit_usr_descriptors, usr_descriptors};
use serde_json::Value;

const COORDINATES: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/rdkit_3d_g1_fixtures.json"
));
const ORACLE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/rdkit_usr_3d_2026_09_1.json"
));

#[test]
fn rdkit_usr_matches_frozen_explicit_h_coordinates() {
    let fixtures: Value = serde_json::from_str(COORDINATES).unwrap();
    let oracle: Value = serde_json::from_str(ORACLE).unwrap();
    assert_eq!(oracle["rdkit_version"], "2026.09.1");
    let cases = fixtures["fixtures"].as_array().unwrap();
    let rows = oracle["rows"].as_array().unwrap();
    assert_eq!(cases.len(), 24);
    assert_eq!(rows.len(), cases.len());
    let mut max_error = 0.0f64;
    let mut different_from_legacy = 0;
    for (case, row) in cases.iter().zip(rows) {
        assert_eq!(case["name"], row["name"]);
        assert_eq!(case["conf_id"], row["conf_id"]);
        let coords: Vec<[f64; 3]> = case["atoms"]
            .as_array()
            .unwrap()
            .iter()
            .map(|atom| {
                [
                    atom["x"].as_f64().unwrap(),
                    atom["y"].as_f64().unwrap(),
                    atom["z"].as_f64().unwrap(),
                ]
            })
            .collect();
        let actual = rdkit_usr_descriptors(&coords).unwrap();
        let legacy = usr_descriptors(&coords);
        for (index, expected) in row["expected"].as_array().unwrap().iter().enumerate() {
            let error = (actual[index] - expected.as_f64().unwrap()).abs();
            max_error = max_error.max(error);
            assert!(
                error < 1e-9,
                "{}[{}] {index}: {error}",
                row["name"],
                row["conf_id"]
            );
        }
        if actual
            .iter()
            .zip(legacy)
            .any(|(left, right)| (left - right).abs() > 1e-6)
        {
            different_from_legacy += 1;
        }
    }
    assert_eq!(different_from_legacy, rows.len());
    eprintln!("RDKit USR maximum absolute error across 24 conformers: {max_error:.3e}");
}

#[test]
fn rdkit_usr_rejects_missing_or_nonfinite_geometry() {
    assert_eq!(
        rdkit_usr_descriptors(&[[0.0; 3]; 2]),
        Err(RdkitUsrError::TooFewAtoms { count: 2 })
    );
    assert_eq!(
        rdkit_usr_descriptors(&[[0.0; 3], [1.0, f64::NAN, 0.0], [2.0; 3]]),
        Err(RdkitUsrError::NonFiniteCoordinate { atom_index: 1 })
    );
}
