use chematic_mol::gaussian::{
    GaussianError, GaussianParseLimits, parse_gaussian_log, parse_gaussian_log_with_limits,
};

fn log(rows: &str) -> String {
    format!(
        " Standard orientation:\n -----\n Center Atomic Coordinates\n -----\n{rows}\n -----\n SCF Done: E(RHF) = -40.25 A.U.\n"
    )
}

#[test]
fn gaussian_03_five_column_geometry_is_read_in_angstroms() {
    let result = parse_gaussian_log(&log(" 1 6 0.0 1.25 -2.5\n 2 8 2.0 0.0 0.5")).unwrap();
    assert_eq!(result.mol.atom_count(), 2);
    assert_eq!(result.coords, [(0.0, 1.25, -2.5), (2.0, 0.0, 0.5)]);
    assert_eq!(result.scf_energy, Some(-40.25));
}

#[test]
fn gaussian_log_rejects_invalid_elements_coordinates_and_atom_limits() {
    for (row, expected) in [
        (
            "1 999 0 0.0 0.0 0.0",
            GaussianError::UnknownElement("999".into()),
        ),
        (
            "1 xx 0 0.0 0.0 0.0",
            GaussianError::UnknownElement("xx".into()),
        ),
        (
            "1 6 0 x 0.0 0.0",
            GaussianError::InvalidCoordinate("x".into()),
        ),
        (
            "1 6 0 0.0 y 0.0",
            GaussianError::InvalidCoordinate("y".into()),
        ),
        (
            "1 6 0 0.0 0.0 z",
            GaussianError::InvalidCoordinate("z".into()),
        ),
        (
            "1 6 0 NaN 0.0 0.0",
            GaussianError::NonFiniteCoordinate("NaN 0.0 0.0".into()),
        ),
    ] {
        assert!(
            matches!(parse_gaussian_log(&log(row)), Err(e) if e == expected),
            "{row}"
        );
    }
    let limits = GaussianParseLimits {
        max_atoms: 1,
        ..Default::default()
    };
    assert!(matches!(
        parse_gaussian_log_with_limits(&log("1 6 0 0 0 0\n2 8 0 1 0 0"), &limits),
        Err(GaussianError::ResourceLimit {
            resource: "atoms",
            actual: 2,
            limit: 1
        })
    ));
}
