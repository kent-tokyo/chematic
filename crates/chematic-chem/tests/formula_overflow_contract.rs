use chematic_chem::{FormulaParseError, parse_formula};

#[test]
fn atom_counts_are_checked_when_elements_or_groups_are_combined() {
    for formula in ["C4294967295C", "C4294967295(C)"] {
        let error = parse_formula(formula).unwrap_err();
        assert_eq!(error, FormulaParseError::InvalidCount("1".into()));
        assert_eq!(error.to_string(), "invalid atom count: 1");
    }
    assert_eq!(parse_formula("C4294967295").unwrap()["C"], u32::MAX);
}
