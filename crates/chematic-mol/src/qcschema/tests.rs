use super::*;

fn minimal_molecule() -> &'static str {
    r#"{"symbols":["H"],"geometry":[0.0,0.0,0.0]}"#
}

#[test]
fn bounded_molecule_parser_rejects_oversized_input() {
    let limits = QcSchemaParseLimits {
        max_input_bytes: 8,
        ..Default::default()
    };
    let err = parse_qcschema_molecule_with_limits(minimal_molecule(), &limits).unwrap_err();
    assert!(matches!(
        err,
        QcSchemaError::ResourceLimit {
            resource: "input bytes",
            ..
        }
    ));
}

#[test]
fn bounded_molecule_parser_rejects_oversized_arrays() {
    let limits = QcSchemaParseLimits {
        max_array_items: 2,
        ..Default::default()
    };
    let err = parse_qcschema_molecule_with_limits(minimal_molecule(), &limits).unwrap_err();
    assert!(matches!(
        err,
        QcSchemaError::ResourceLimit {
            resource: "array items",
            ..
        }
    ));
}

#[test]
fn bounded_atomic_parsers_apply_json_limits_before_schema_validation() {
    let limits = QcSchemaParseLimits {
        max_string_bytes: 8,
        ..Default::default()
    };
    let input = r#"{"x":"longvalue","symbols":["H"],"geometry":[0.0,0.0,0.0]}"#;
    assert!(matches!(
        parse_atomic_input_with_limits(input, &limits),
        Err(QcSchemaError::ResourceLimit {
            resource: "string bytes",
            ..
        })
    ));
    assert!(matches!(
        parse_atomic_result_with_limits(input, &limits),
        Err(QcSchemaError::ResourceLimit {
            resource: "string bytes",
            ..
        })
    ));
}
