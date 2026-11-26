use chematic_mol::{MolJsonError, MolJsonParseLimits, parse_moljson_with_limits};

#[test]
fn string_budget_counts_utf8_bytes_even_in_unknown_fields() {
    let input = r#"{"atoms":[],"bonds":[],"extra":"ééé"}"#;
    let mut limits = MolJsonParseLimits {
        max_string_bytes: 5,
        ..Default::default()
    };
    assert!(matches!(
        parse_moljson_with_limits(input, &limits),
        Err(MolJsonError::ResourceLimit {
            resource: "string bytes",
            actual: 6,
            limit: 5
        })
    ));
    limits.max_string_bytes = 6;
    assert!(parse_moljson_with_limits(input, &limits).is_ok());
}
