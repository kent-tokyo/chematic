use chematic_mol::{MoldenLimits, parse_molden, parse_molden_with_limits};

#[test]
fn malformed_geometry_and_vibration_fields_identify_the_failed_boundary() {
    let cases = [
        (
            "[Molden Format]\n[Atoms\nC 1 6 0 0 0\n",
            "unclosed section header",
        ),
        (
            "[Molden Format]\n[Atoms] Angs\nC 1 6 0 0\n",
            "Atoms row needs symbol",
        ),
        (
            "[Molden Format]\n[Atoms] Angs\nC x 6 0 0 0\n",
            "invalid atom index",
        ),
        (
            "[Molden Format]\n[Atoms] Angs\nC 1 x 0 0 0\n",
            "invalid atomic number",
        ),
        (
            "[Molden Format]\n[Atoms] Angs\nC 1 6 bad 0 0\n",
            "invalid number",
        ),
        (
            "[Molden Format]\n[FR-COORD]\nC 0 0\n",
            "FR-COORD row needs symbol",
        ),
        (
            "[Molden Format]\n[FREQ]\n100 200\n",
            "frequency row needs exactly one value",
        ),
        (
            "[Molden Format]\n[FR-NORM-COORD]\nvibration 1\n0 0 0\n",
            "normal modes require FR-COORD",
        ),
        (
            "[Molden Format]\n[FR-COORD]\nC 0 0 0\n[FR-NORM-COORD]\nmode 1\n0 0 0\n",
            "expected vibration and mode index",
        ),
        (
            "[Molden Format]\n[FR-COORD]\nC 0 0 0\n[FR-NORM-COORD]\nvibration x\n0 0 0\n",
            "invalid mode index",
        ),
    ];
    for (source, expected) in cases {
        let error = match parse_molden(source) {
            Ok(_) => panic!("malformed input accepted: {source}"),
            Err(error) => error,
        };
        assert!(error.detail.contains(expected), "{source}: {error}");
    }
}

#[test]
fn atom_intensity_and_total_row_limits_are_enforced_before_output() {
    let limits = MoldenLimits {
        max_atoms: 0,
        max_modes: 0,
        max_displacements: 0,
        ..Default::default()
    };
    for (source, message) in [
        (
            "[Molden Format]\n[Atoms] Angs\nC 1 6 0 0 0\n".to_string(),
            "Atoms count exceeds limit",
        ),
        (
            "[Molden Format]\n[INT]\n1.0\n".to_string(),
            "intensity count exceeds limit",
        ),
        (
            format!("[Molden Format]\n[FREQ]\n{}", "1.0\n".repeat(33)),
            "record count exceeds limits",
        ),
    ] {
        let error = match parse_molden_with_limits(&source, &limits) {
            Ok(_) => panic!("limit should be enforced"),
            Err(error) => error,
        };
        assert!(error.detail.contains(message), "{error}");
    }
}
