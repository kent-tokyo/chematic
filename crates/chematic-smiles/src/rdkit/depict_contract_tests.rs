//! Coordinate fixtures measured with RDKit 2026.03.1, native RDDepictor,
//! Compute2DCoords defaults and forceRDKit=True (no CoordGen or templates).
use super::rdkit_2d_coords;
#[test]
fn ring_and_nontetrahedral_depictions_match_rdkit_reference() {
    let fixtures: serde_json::Value = serde_json::from_str(include_str!(
        "../../../../validation/rdkit-2026.03.1-depict-contract.json"
    ))
    .unwrap();
    for row in fixtures.as_array().unwrap() {
        let smiles = row["smiles"].as_str().unwrap();
        let mol = crate::parse(smiles).unwrap();
        let actual = rdkit_2d_coords(&mol).unwrap();
        let expected = row["coords"].as_array().unwrap();
        assert_eq!(actual.len(), expected.len());
        for (i, point) in actual.iter().enumerate() {
            for axis in 0..2 {
                let reference = expected[i][axis].as_f64().unwrap();
                assert!(
                    (point[axis] - reference).abs() < 1e-9,
                    "{smiles} atom {i} axis {axis}: {} != {reference}",
                    point[axis]
                );
            }
        }
    }
}
