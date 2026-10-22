use chematic_mol::{
    CjsonParseLimits, parse_cjson, parse_cjson_checked, parse_cjson_checked_with_limits,
    parse_cjson_record, write_cjson, write_cjson_record,
};
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!("../../../validation/cjson_checked/cases.json")).unwrap()
}
fn layouts(input: &str) -> [String; 3] {
    // Compact without parsing: negative fixtures deliberately contain duplicate keys.
    let mut compact = String::new();
    let mut quoted = false;
    let mut escaped = false;
    for c in input.chars() {
        if quoted {
            compact.push(c);
            if escaped {
                escaped = false;
            } else if c == '\\' {
                escaped = true;
            } else if c == '"' {
                quoted = false;
            }
        } else if c == '"' {
            quoted = true;
            compact.push(c);
        } else if !c.is_ascii_whitespace() {
            compact.push(c);
        }
    }
    [input.to_owned(), input.replace('\n', "\r\n"), compact]
}
#[test]
fn complete_and_compatibility_graphs_keep_source_indexed_values() {
    let packet = packet();
    for case in packet["valid"]
        .as_array()
        .unwrap()
        .iter()
        .chain(packet["compatibility"].as_array().unwrap())
    {
        for source in layouts(case["input"].as_str().unwrap()) {
            let (mol, coords) = parse_cjson_checked(&source).unwrap();
            let output: Value = serde_json::from_str(&write_cjson(&mol, &coords)).unwrap();
            let expected = &case["expected"];
            assert_eq!(
                output["atoms"]["elements"]["number"], expected["numbers"],
                "{}",
                case["id"]
            );
            assert_eq!(
                output["atoms"]
                    .get("formalCharges")
                    .cloned()
                    .unwrap_or(serde_json::json!(vec![0; mol.atom_count()])),
                expected["charges"]
            );
            assert_eq!(
                output["atoms"]
                    .get("isotopes")
                    .cloned()
                    .unwrap_or(serde_json::json!(vec![0; mol.atom_count()])),
                expected["isotopes"]
            );
            assert_eq!(output["bonds"]["connections"]["index"], expected["indices"]);
            let orders = |v: &Value| {
                v.as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_f64().unwrap())
                    .collect::<Vec<_>>()
            };
            assert_eq!(
                orders(&output["bonds"]["order"]),
                orders(&expected["orders"])
            );
            let record = parse_cjson_record(&source).unwrap();
            assert_eq!(write_cjson_record(&record), source);
            assert_eq!(record.molecule().atom_count(), mol.atom_count());
            assert_eq!(record.coordinates(), coords);
        }
    }
}
#[test]
fn invalid_declared_fields_never_return_a_candidate() {
    for case in packet()["invalid"].as_array().unwrap() {
        for source in layouts(case["input"].as_str().unwrap()) {
            assert!(parse_cjson_checked(&source).is_err(), "{}", case["id"]);
            assert!(parse_cjson_record(&source).is_err(), "{}", case["id"]);
        }
    }
}
#[test]
fn legacy_reader_and_opaque_retention_remain_available() {
    let source = packet()["invalid"][0]["input"].as_str().unwrap().to_owned();
    assert!(parse_cjson(&source).is_ok());
    assert!(parse_cjson_checked(&source).is_err());
    let source = r#"{"atoms":{"elements":{"number":[6]}},"vendor":{"note":"not interpreted"}}"#;
    let record = parse_cjson_record(source).unwrap();
    assert_eq!(write_cjson_record(&record), source);
    assert_eq!(parse_cjson_checked(source).unwrap().0.atom_count(), 1);
}
#[test]
fn bounded_projection_and_numeric_coordinates_match_record_reader() {
    let source = packet()["valid"][0]["input"].as_str().unwrap().to_owned();
    let exact = CjsonParseLimits {
        max_input_bytes: source.len(),
        max_atoms: 2,
        max_bonds: 1,
        ..Default::default()
    };
    assert!(parse_cjson_checked_with_limits(&source, &exact).is_ok());
    for limits in [
        CjsonParseLimits {
            max_input_bytes: source.len() - 1,
            ..exact
        },
        CjsonParseLimits {
            max_atoms: 1,
            ..exact
        },
        CjsonParseLimits {
            max_bonds: 0,
            ..exact
        },
        CjsonParseLimits {
            max_array_items: 2,
            ..exact
        },
        CjsonParseLimits {
            max_json_depth: 1,
            ..exact
        },
        CjsonParseLimits {
            max_string_bytes: 2,
            ..exact
        },
    ] {
        assert!(matches!(
            parse_cjson_checked_with_limits(&source, &limits),
            Err(chematic_mol::CjsonError::ResourceLimit { .. })
        ));
    }
    let input = r#"{"atoms":{"elements":{"number":[6]},"coords":{"3d":[-1.4175225513503917,-0.0,5e-324]}}}"#;
    let coords = parse_cjson_checked(input).unwrap().1;
    let expected = [-1.4175225513503917f64, -0.0, 5e-324];
    for (a, b) in [coords[0].0, coords[0].1, coords[0].2].iter().zip(expected) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}
