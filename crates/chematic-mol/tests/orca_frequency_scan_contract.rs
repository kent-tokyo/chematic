use chematic_mol::parse_orca_output;

#[test]
fn frequency_table_stops_after_bounded_nonmatching_headers() {
    // A frequency table may have headings before its first numeric entry.
    let headers = (0..10).map(|i| format!("header {i}\n")).collect::<String>();
    let prefix = format!("VIBRATIONAL FREQUENCIES\n-----------------------\n{headers}");
    let accepted = parse_orca_output(&format!(
        "{prefix}0: -123.45 cm**-1 imaginary\n\nFINAL SINGLE POINT ENERGY -12.5\n"
    ))
    .unwrap();
    assert_eq!(accepted.frequencies_cm1, [-123.45]);
    assert_eq!(accepted.final_energy_hartree, Some(-12.5));

    // An eleventh nonmatching line ends the search; later numerical text
    // must not be misread as a frequency table entry.
    let stopped = parse_orca_output(&format!(
        "{prefix}extra heading\n0: -123.45 cm**-1 imaginary\nFINAL SINGLE POINT ENERGY -12.5\n"
    ))
    .unwrap();
    assert!(stopped.frequencies_cm1.is_empty());
    assert_eq!(stopped.final_energy_hartree, Some(-12.5));
}
