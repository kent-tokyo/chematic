use chematic_mol::{
    GprParseLimits, GprRecord, gpr_nominal_value, parse_gpr, parse_gpr_json, parse_gpr_with_limits,
    write_gpr, write_gpr_with_limits,
};
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/gpr_preservation/cases.json"
    ))
    .unwrap()
}
#[test]
fn fields_tokens_bond_codes_and_absent_charges_survive_three_cycles() {
    for c in packet()["valid"].as_array().unwrap() {
        let expected: GprRecord = serde_json::from_value(c["expected"].clone()).unwrap();
        let mut record = parse_gpr(c["text"].as_str().unwrap()).unwrap();
        assert_eq!(record, expected, "{}", c["id"]);
        for _ in 0..3 {
            record = parse_gpr(&write_gpr(&record).unwrap()).unwrap();
            assert_eq!(record, expected, "{}", c["id"]);
        }
    }
    assert_eq!(
        gpr_nominal_value("-0.000").unwrap().to_bits(),
        (-0.0f64).to_bits()
    );
    assert_eq!(gpr_nominal_value(".125").unwrap(), 0.125);
}
#[test]
fn unsupported_input_and_both_direction_limits_refuse() {
    for c in packet()["refused"].as_array().unwrap() {
        assert!(
            parse_gpr(c["text"].as_str().unwrap()).is_err(),
            "{}",
            c["id"]
        );
    }
    let text = packet()["valid"][2]["text"].as_str().unwrap().to_string();
    let record = parse_gpr(&text).unwrap();
    for limits in [
        GprParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        GprParseLimits {
            max_bonds: 0,
            ..Default::default()
        },
        GprParseLimits {
            max_input_bytes: 10,
            ..Default::default()
        },
        GprParseLimits {
            max_line_len: 5,
            ..Default::default()
        },
    ] {
        assert!(parse_gpr_with_limits(&text, &limits).is_err());
        assert!(write_gpr_with_limits(&record, &limits).is_err());
    }
}
#[test]
fn edited_writer_and_strict_json_never_infer_unsupported_fields() {
    let mut record = parse_gpr(packet()["valid"][2]["text"].as_str().unwrap()).unwrap();
    record.bonds[0].order = chematic_mol::GprBondOrder::D;
    record.coordinates[0].position_nm[0] = ".12345678901234567".into();
    record.partial_charges.as_mut().unwrap()[0].value = "-1.2345678901234567".into();
    assert_eq!(parse_gpr(&write_gpr(&record).unwrap()).unwrap(), record);
    let json = serde_json::to_string(&record).unwrap();
    let limits = GprParseLimits::default();
    assert_eq!(parse_gpr_json(&json, &limits).unwrap(), record);
    for bad in [
        json.replacen("\"atoms\":", "\"atoms\":[],\"atoms\":", 1),
        json.replacen("\"atoms\":", "\"at\\u006fms\":[],\"atoms\":", 1),
        json.replacen("\"id\":", "\"id\":0,\"id\":", 1),
        json.replacen("\"id\":", "\"formal_charge\":1,\"id\":", 1),
        json.replacen("\"atoms\":", "\"stereo\":null,\"atoms\":", 1),
        json.replace("\".12345678901234567\"", "0.125"),
        json.replace("\".12345678901234567\"", "\"1e-999\""),
        json.replace("\"atomic_number\":6", "\"atomic_number\":6.5"),
    ] {
        assert!(parse_gpr_json(&bad, &limits).is_err(), "{bad}");
    }
    record.partial_charges = None;
    assert!(!write_gpr(&record).unwrap().contains("!Charges"));
    record.coordinates[0].position_nm[0] = "NaN".into();
    assert!(write_gpr(&record).is_err());
}
