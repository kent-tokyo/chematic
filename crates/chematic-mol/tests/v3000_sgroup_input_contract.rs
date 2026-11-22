use chematic_mol::{V3000SGroupKind, parse_v3000_sgroup_line};

#[test]
fn generic_sgroup_preserves_atom_ids_parent_and_attributes() {
    let parsed = parse_v3000_sgroup_line("7 GEN 3 ATOMS=(2 1 9) LABEL=group").unwrap();
    assert_eq!(parsed.id, 7);
    assert_eq!(parsed.kind, V3000SGroupKind::Gen);
    assert_eq!(parsed.parent_id, Some(3));
    assert_eq!(parsed.atom_ids, [1, 9]);
    assert!(
        parsed
            .attributes
            .iter()
            .any(|(name, value)| name == "LABEL" && value == "group")
    );
}

#[test]
fn malformed_sgroup_lists_and_identifiers_are_rejected() {
    for (line, expected) in [
        (
            "x SUP 0 ATOMS=(1 1)",
            "SGROUP id is not a non-negative integer",
        ),
        (
            "1 SUP x ATOMS=(1 1)",
            "SGROUP parent is not a non-negative integer",
        ),
        ("1 SUP 0 ATOMS=(1 1", "SGROUP ATOMS value is not closed"),
        (
            "1 SUP 0 ATOMS=()",
            "SGROUP ATOMS value is missing its count",
        ),
        (
            "1 SUP 0 ATOMS=(x 1)",
            "SGROUP ATOMS count is not an integer",
        ),
        (
            "1 SUP 0 ATOMS=(1 x)",
            "SGROUP ATOMS contains a non-integer atom id",
        ),
        ("1 SUP 0 ATOMS=(2 1)", "declares 2 ids but contains 1"),
    ] {
        let error = parse_v3000_sgroup_line(line).unwrap_err();
        assert!(error.to_string().contains(expected), "{line}: {error}");
    }
}
