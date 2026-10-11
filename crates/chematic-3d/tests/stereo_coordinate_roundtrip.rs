//! Serialized geometry must preserve the exact coordinates used for stereo checks.
use chematic_3d::{
    coords::{Coords3D, Point3},
    stereo_constraints::verify_stereo,
};
use chematic_core::AtomIdx;

#[test]
fn json_coordinates_preserve_boundary_stereo() {
    // The four alkene reference atoms from chembl_tier_b_0126. Without
    // float_roundtrip, JSON decoding changes coordinate bits and flips the
    // E/Z report at the 90-degree boundary. Do not soften the stereo verifier.
    let mol = chematic_smiles::parse("C/C=C/C").unwrap();
    let literal = [
        [-1.6994382329723745, -3.412312097402026, -2.204046205120537],
        [
            -1.8533642073648884,
            -2.8083132813418286,
            -0.8321184001674955,
        ],
        [-2.9297120835404935, -1.994864722790872, -0.7343006938254443],
        [
            -3.1094875313017245,
            -0.5127242579901362,
            -0.9733856679237677,
        ],
    ];
    // Literal decimal input mirrors Python json.dumps, avoiding any Rust
    // serializer/parser coupling in the regression.
    let json = r#"[
        [-1.6994382329723745,-3.412312097402026,-2.204046205120537],
        [-1.8533642073648884,-2.8083132813418286,-0.8321184001674955],
        [-2.9297120835404935,-1.994864722790872,-0.7343006938254443],
        [-3.1094875313017245,-0.5127242579901362,-0.9733856679237677]
    ]"#;
    let parsed: Vec<[f64; 3]> = serde_json::from_str(json).unwrap();
    let value: serde_json::Value = serde_json::from_str(json).unwrap();
    let make_coords = |points: &[[f64; 3]]| {
        let mut coords = Coords3D::new_zeroed(points.len());
        for (i, p) in points.iter().enumerate() {
            coords.set(AtomIdx(i as u32), Point3::new(p[0], p[1], p[2]));
        }
        coords
    };
    for (i, point) in literal.iter().enumerate() {
        for (j, component) in point.iter().enumerate() {
            let expected = f64::to_bits(*component);
            assert_eq!(parsed[i][j].to_bits(), expected);
            // The common scorer parses Value before reading coordinates.
            assert_eq!(value[i][j].as_f64().unwrap().to_bits(), expected);
        }
    }
    assert!(verify_stereo(&mol, &make_coords(&literal)).is_fully_satisfied());
    assert!(verify_stereo(&mol, &make_coords(&parsed)).is_fully_satisfied());
}
