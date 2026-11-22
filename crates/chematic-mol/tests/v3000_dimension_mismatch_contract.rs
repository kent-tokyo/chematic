use chematic_mol::{
    CoordinateDimension, GeometryRank, Stereo3DDiagnostic, read_mol_v3000_with_diagnostics,
};

fn block(dimension: &str, fourth_z: i32) -> String {
    format!(
        "mismatch\n  test              {dimension}\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 4 3 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 F 1 0 0 0\nM  V30 3 Cl 0 1 0 0\nM  V30 4 Br 0 0 {fourth_z} 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 1 3\nM  V30 3 1 1 4\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n"
    )
}

#[test]
fn v3000_reports_declared_2d_with_nonzero_3d_coordinates() {
    let report = read_mol_v3000_with_diagnostics(&block("2D", 1)).unwrap();
    assert_eq!(report.coordinate_dimension, CoordinateDimension::TwoD);
    assert_eq!(report.geometry_rank, GeometryRank::ThreeD);
    assert!(report.stereo3d_diagnostics.iter().any(|d| matches!(
        d,
        Stereo3DDiagnostic::DeclaredTwoDButNonzeroZ {
            observed: GeometryRank::ThreeD
        }
    )));
    assert!(report.conformer.is_some());
}
