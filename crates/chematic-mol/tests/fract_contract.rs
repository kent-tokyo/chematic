use chematic_mol::{
    FractNumber, FractParseLimits, FractRecord, parse_fract, parse_fract_json,
    parse_fract_with_limits, write_fract, write_fract_with_limits,
};
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/fract_preservation/cases.json"
    ))
    .unwrap()
}
#[test]
fn authored_fields_and_nominal_bits_survive_three_cycles() {
    for case in packet()["valid"].as_array().unwrap() {
        let mut record = parse_fract(case["text"].as_str().unwrap()).unwrap();
        let expected: FractRecord = serde_json::from_value(case["expected"].clone()).unwrap();
        assert_eq!(record, expected, "{}", case["id"]);
        for _ in 0..3 {
            record = parse_fract(&write_fract(&record).unwrap()).unwrap();
            assert_eq!(record, expected, "{}", case["id"]);
            for (a, b) in record.atoms.iter().zip(&expected.atoms) {
                for (x, y) in a.fractional.iter().zip(&b.fractional) {
                    assert_eq!(x.value().to_bits(), y.value().to_bits());
                }
            }
        }
    }
    let zero = FractNumber::new("-0.000(02)").unwrap();
    assert_eq!(zero.value().to_bits(), (-0.0f64).to_bits());
    assert_eq!(zero.uncertainty_digits(), Some("02"));
}
#[test]
fn unsupported_records_are_errors_and_limits_cover_both_directions() {
    for case in packet()["refused"].as_array().unwrap() {
        assert!(
            parse_fract(case["text"].as_str().unwrap()).is_err(),
            "{}",
            case["id"]
        );
    }
    assert_eq!(
        parse_fract("line\n4 5 6 90 90 90\n\nXx 0 0 0\n")
            .unwrap_err()
            .line,
        4
    );
    let text = "limit\n4 5 6 90 90 90\nC .1 .2 .3\n";
    let record = parse_fract(text).unwrap();
    for limits in [
        FractParseLimits {
            max_atoms: 0,
            ..Default::default()
        },
        FractParseLimits {
            max_input_bytes: 10,
            ..Default::default()
        },
        FractParseLimits {
            max_line_len: 5,
            ..Default::default()
        },
    ] {
        assert!(parse_fract_with_limits(text, &limits).is_err());
        assert!(write_fract_with_limits(&record, &limits).is_err());
    }
}
#[test]
fn writer_uses_edits_and_json_refuses_unrepresented_meaning() {
    let mut record = parse_fract("old\n4 5 6 90 90 90\nC .1 .2 .3\n").unwrap();
    record.title = "edited".into();
    record.atoms[0].fractional[0] = FractNumber::new("-2.12345678901234567(004)").unwrap();
    record.cell[0] = FractNumber::new("7.25(2)").unwrap();
    assert_eq!(parse_fract(&write_fract(&record).unwrap()).unwrap(), record);
    let limits = FractParseLimits::default();
    let json = serde_json::to_string(&record).unwrap();
    assert_eq!(parse_fract_json(&json, &limits).unwrap(), record);
    for invalid in [
        json.replacen("\"title\":", "\"title\":\"a\",\"title\":", 1),
        json.replacen("\"title\":", "\"t\\u0069tle\":\"a\",\"title\":", 1),
        json.replacen("\"element\":", "\"element\":\"O\",\"element\":", 1),
        json.replacen("\"element\":", "\"charge\":1,\"element\":", 1),
        json.replacen("\"title\":", "\"bonds\":[],\"title\":", 1),
        json.replace("\"7.25(2)\"", "7.25"),
        json.replace("\"7.25(2)\"", "null"),
        json.replace("\"7.25(2)\"", "\"NaN\""),
    ] {
        assert!(parse_fract_json(&invalid, &limits).is_err(), "{invalid}");
    }
    record.atoms[0].element = "C1".into();
    assert!(write_fract(&record).is_err());
    record.atoms[0].element = "C".into();
    record.title.push('\n');
    assert!(write_fract(&record).is_err());
}
