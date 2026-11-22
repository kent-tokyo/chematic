use chematic_crystal::{CrystalError, FractionalCoord, Lattice};

#[test]
fn finite_vectors_with_overflowing_volume_are_rejected() {
    let error = Lattice::from_matrix([[1e150, 0.0, 0.0], [0.0, 1e150, 0.0], [0.0, 0.0, 1e150]])
        .unwrap_err();
    assert!(matches!(error, CrystalError::NonFiniteVolume));
}

#[test]
fn fractional_coordinates_report_nonfinite_input_before_wrapping() {
    let coord = FractionalCoord::new([0.0, f64::INFINITY, 0.0]);
    assert!(!coord.is_finite());
    assert!(!coord.wrapped().is_finite());
}

#[test]
fn orthorhombic_cell_reports_the_nonpositive_axis() {
    for (lengths, expected_axis) in [
        ([0.0, 3.0, 4.0], "a"),
        ([2.0, -1.0, 4.0], "b"),
        ([2.0, 3.0, 0.0], "c"),
    ] {
        let error = Lattice::orthorhombic(lengths[0], lengths[1], lengths[2]).unwrap_err();
        assert!(
            matches!(error, CrystalError::NonPositiveLength { axis, .. } if axis == expected_axis)
        );
    }
}
