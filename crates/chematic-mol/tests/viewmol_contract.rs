use chematic_mol::{
    ViewmolParseLimits, ViewmolRecord, parse_viewmol, parse_viewmol_json,
    parse_viewmol_with_limits, write_viewmol, write_viewmol_with_limits,
};
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/viewmol_preservation/cases.json"
    ))
    .unwrap()
}
#[test]
fn raw_units_titles_flags_and_projected_bits_survive_three_cycles() {
    for c in packet()["valid"].as_array().unwrap() {
        let expected: ViewmolRecord = serde_json::from_value(c["expected"].clone()).unwrap();
        let mut r = parse_viewmol(c["text"].as_str().unwrap()).unwrap();
        assert_eq!(r, expected, "{}", c["id"]);
        for _ in 0..3 {
            let p = r.positions_angstrom().unwrap();
            for (actual, reference) in p.iter().zip(c["expected_angstrom"].as_array().unwrap()) {
                for (a, e) in actual.iter().zip(reference.as_array().unwrap()) {
                    assert_eq!(a.to_bits(), e.as_f64().unwrap().to_bits(), "{}", c["id"]);
                }
            }
            r = parse_viewmol(&write_viewmol(&r).unwrap()).unwrap();
            assert_eq!(r, expected, "{}", c["id"]);
        }
    }
}
#[test]
fn invalid_input_and_read_write_projection_limits_refuse() {
    for c in packet()["refused"].as_array().unwrap() {
        assert!(
            parse_viewmol(c["text"].as_str().unwrap()).is_err(),
            "{}",
            c["id"]
        );
    }
    let text = packet()["valid"][0]["text"].as_str().unwrap().to_string();
    let r = parse_viewmol(&text).unwrap();
    for l in [
        ViewmolParseLimits {
            max_atoms: 0,
            ..Default::default()
        },
        ViewmolParseLimits {
            max_input_bytes: 10,
            ..Default::default()
        },
        ViewmolParseLimits {
            max_line_len: 2,
            ..Default::default()
        },
    ] {
        assert!(parse_viewmol_with_limits(&text, &l).is_err());
        assert!(write_viewmol_with_limits(&r, &l).is_err());
    }
    for c in packet()["projection_refused"].as_array().unwrap() {
        let r = parse_viewmol(c["text"].as_str().unwrap()).unwrap();
        assert_eq!(serde_json::to_value(&r).unwrap(), c["expected"]);
        assert_eq!(parse_viewmol(&write_viewmol(&r).unwrap()).unwrap(), r);
        assert!(r.positions_angstrom().is_err());
    }
}
#[test]
fn edited_fields_and_strict_json_preserve_unknown_state_as_refusal() {
    let mut r = parse_viewmol("$coord 2\n.1 .2 .3 C x\n$end\n").unwrap();
    r.title = Some("edited".into());
    r.factor_to_angstrom = ".5".into();
    r.atoms[0].position[0] = "-.12345678901234567".into();
    r.atoms[0].fixed_axes = Some("zy".into());
    assert_eq!(parse_viewmol(&write_viewmol(&r).unwrap()).unwrap(), r);
    let json = serde_json::to_string(&r).unwrap();
    let l = ViewmolParseLimits::default();
    assert_eq!(parse_viewmol_json(&json, &l).unwrap(), r);
    for bad in [
        json.replacen("\"title\":", "\"title\":null,\"title\":", 1),
        json.replacen("\"title\":", "\"t\\u0069tle\":null,\"title\":", 1),
        json.replacen("\"element\":", "\"element\":\"N\",\"element\":", 1),
        json.replacen("\"element\":", "\"charge\":1,\"element\":", 1),
        json.replacen("\"title\":", "\"cell\":[],\"title\":", 1),
        json.replace("\"-.12345678901234567\"", "-.125"),
        json.replace("\".5\"", "0.5"),
        json.replace("\"zy\"", "\"\""),
        json.replace("\"zy\"", "\"xx\""),
    ] {
        assert!(parse_viewmol_json(&bad, &l).is_err(), "{bad}");
    }
    r.title = Some("$coord".into());
    assert!(write_viewmol(&r).is_err());
}
