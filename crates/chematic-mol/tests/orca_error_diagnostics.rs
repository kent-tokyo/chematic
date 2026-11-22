use chematic_mol::parse_orca_cartesian;

#[test]
fn orca_cartesian_refusals_explain_units_and_geometry_requirements() {
    let cases = [
        (
            "! Angs Bohrs\n* xyz 0 1\nH 0 0 0\n*\n",
            "conflicting Angs/Bohrs",
        ),
        (
            "! UnitsAU\n* xyz 0 1\nH 0 0 0\n*\n",
            "unsupported ORCA coordinate unit keyword: UnitsAU",
        ),
        (
            "%coords\nend\n* xyz 0 1\nH 0 0 0\n*\n",
            "does not interpret %coords blocks",
        ),
        ("! HF\n", "requires one nonempty embedded xyz block"),
    ];
    for (source, expected) in cases {
        let error = parse_orca_cartesian(source).unwrap_err();
        assert!(error.to_string().contains(expected), "{}", error);
    }
}
