use chematic_mol::{
    UnixyzParseLimits, UnixyzRecord, parse_unixyz, parse_unixyz_json, parse_unixyz_with_limits,
    unixyz_nominal_value, write_unixyz, write_unixyz_with_limits,
};
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/unixyz_preservation/cases.json"
    ))
    .unwrap()
}
#[test]
fn title_atomic_numbers_and_decimal_tokens_survive_three_cycles() {
    for c in packet()["valid"].as_array().unwrap() {
        let expected: UnixyzRecord = serde_json::from_value(c["expected"].clone()).unwrap();
        let mut r = parse_unixyz(c["text"].as_str().unwrap()).unwrap();
        assert_eq!(r, expected, "{}", c["id"]);
        for _ in 0..3 {
            r = parse_unixyz(&write_unixyz(&r).unwrap()).unwrap();
            assert_eq!(r, expected, "{}", c["id"]);
        }
    }
    assert_eq!(
        unixyz_nominal_value("-0.000").unwrap().to_bits(),
        (-0.0f64).to_bits()
    );
    // This is a different header layout, not an alias for the existing XYZ reader.
    assert!(chematic_mol::xyz::parse_xyz("title\n1\n6 .1 .2 .3\n").is_err());
}
#[test]
fn malformed_whole_input_and_both_direction_limits_refuse() {
    for c in packet()["refused"].as_array().unwrap() {
        assert!(
            parse_unixyz(c["text"].as_str().unwrap()).is_err(),
            "{}",
            c["id"]
        );
    }
    let text = packet()["valid"][0]["text"].as_str().unwrap().to_string();
    let r = parse_unixyz(&text).unwrap();
    for l in [
        UnixyzParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        UnixyzParseLimits {
            max_input_bytes: 10,
            ..Default::default()
        },
        UnixyzParseLimits {
            max_line_len: 3,
            ..Default::default()
        },
    ] {
        assert!(parse_unixyz_with_limits(&text, &l).is_err());
        assert!(write_unixyz_with_limits(&r, &l).is_err());
    }
}
#[test]
fn edited_fields_and_strict_json_do_not_infer_metadata() {
    let mut r = parse_unixyz("old\n1\n6 .1 .2 .3\n").unwrap();
    r.title = "edited".into();
    r.atoms[0].atomic_number = 8;
    r.atoms[0].position_angstrom[0] = "-.12345678901234567".into();
    assert_eq!(parse_unixyz(&write_unixyz(&r).unwrap()).unwrap(), r);
    let json = serde_json::to_string(&r).unwrap();
    let l = UnixyzParseLimits::default();
    assert_eq!(parse_unixyz_json(&json, &l).unwrap(), r);
    for bad in [
        json.replacen("\"title\":", "\"title\":\"a\",\"title\":", 1),
        json.replacen("\"title\":", "\"t\\u0069tle\":\"a\",\"title\":", 1),
        json.replacen(
            "\"atomic_number\":",
            "\"atomic_number\":6,\"atomic_number\":",
            1,
        ),
        json.replacen("\"atomic_number\":", "\"charge\":1,\"atomic_number\":", 1),
        json.replacen("\"title\":", "\"bonds\":[],\"title\":", 1),
        json.replace("\"-.12345678901234567\"", "-0.125"),
        json.replace("\"-.12345678901234567\"", "\"1e-999\""),
        json.replace("\"atomic_number\":8", "\"atomic_number\":8.5"),
    ] {
        assert!(parse_unixyz_json(&bad, &l).is_err(), "{bad}");
    }
    r.title.push('\n');
    assert!(write_unixyz(&r).is_err());
    r.title.clear();
    r.atoms[0].position_angstrom[0] = "NaN".into();
    assert!(write_unixyz(&r).is_err());
}
