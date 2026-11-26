use chematic_mol::parse_orca_output;

#[test]
fn truncated_geometry_rows_do_not_hide_later_energy() {
    for row in ["H 1.0 2.0", "Xx 1.0 2.0 3.0", "H 1.0 missing 3.0"] {
        let text = format!(
            "CARTESIAN COORDINATES (ANGSTROEM)\n--------------------\nC 0.0 0.0 0.0\n{row}\nFINAL SINGLE POINT ENERGY -12.5\n"
        );
        let parsed = parse_orca_output(&text).unwrap();
        assert_eq!(parsed.trajectory.len(), 1, "{row}");
        assert_eq!(parsed.trajectory[0].coords, [(0.0, 0.0, 0.0)], "{row}");
        assert_eq!(parsed.final_energy_hartree, Some(-12.5), "{row}");
    }
}

#[test]
fn malformed_energy_does_not_replace_the_last_valid_value() {
    let text = "FINAL SINGLE POINT ENERGY -12.5\nFINAL SINGLE POINT ENERGY pending\n";
    let parsed = parse_orca_output(text).unwrap();
    assert_eq!(parsed.final_energy_hartree, Some(-12.5));
}

#[test]
fn frequency_table_does_not_accept_invalid_numbers_or_units() {
    for row in ["1: 123.0 Hz", "1: unknown cm**-1"] {
        let text = format!(
            "VIBRATIONAL FREQUENCIES\n-----------------------\n0: 100.0 cm**-1\n{row}\nFINAL SINGLE POINT ENERGY -12.5\n"
        );
        let parsed = parse_orca_output(&text).unwrap();
        assert_eq!(parsed.frequencies_cm1, [100.0], "{row}");
        assert_eq!(parsed.final_energy_hartree, Some(-12.5), "{row}");
    }
}
