use chematic_mol::{
    DlpolyConfigLimits, dlpoly_config_from_json, parse_dlpoly_config,
    parse_dlpoly_config_with_limits, write_dlpoly_config,
};
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!("../../../validation/dlpoly_config/cases.json")).unwrap()
}

#[test]
fn shared_records_preserve_all_declared_fields_and_source_order() {
    for case in packet()["cases"].as_array().unwrap() {
        let record = parse_dlpoly_config(case["input"].as_str().unwrap()).unwrap();
        assert_eq!(
            serde_json::to_value(&record).unwrap(),
            case["expected"],
            "{}",
            case["id"]
        );
        let written = write_dlpoly_config(&record).unwrap();
        assert!(written.lines().all(|line| line.len() <= 72));
        assert_eq!(
            parse_dlpoly_config(&written).unwrap(),
            record,
            "{}",
            case["id"]
        );
        assert_eq!(
            dlpoly_config_from_json(&serde_json::to_string(&record).unwrap()).unwrap(),
            record
        );
    }
}

#[test]
fn invalid_and_unsupported_records_are_refused() {
    for case in packet()["refusals"].as_array().unwrap() {
        assert!(
            parse_dlpoly_config(case["input"].as_str().unwrap()).is_err(),
            "{}",
            case["id"]
        );
    }
}

#[test]
fn writer_serializes_edits_without_replaying_source_or_zeroing_vectors() {
    let mut record = parse_dlpoly_config(include_str!(
        "../../../validation/dlpoly_config/vectors.CONFIG"
    ))
    .unwrap();
    record.title = "edited source units".into();
    record.atoms[0].name = "DPD-A".into();
    record.atoms[0].velocity = Some([1.2345678901234567, 0.0, -2.5]);
    record.atoms[0].force = Some([-123.45678901234567, 1e-100, 2.5]);
    record.atoms[0].header_extra = vec!["opaque".into(), "-0.1250".into()];
    record.atoms.reverse();
    record.imcon = 3;
    record.cell_vectors = Some([[10., 1., 0.], [2., 11., 0.], [0., 0., -12.]]);
    assert_eq!(
        parse_dlpoly_config(&write_dlpoly_config(&record).unwrap()).unwrap(),
        record
    );
}

#[test]
fn mutable_records_cannot_bypass_presence_index_and_finiteness_checks() {
    let original = parse_dlpoly_config(include_str!(
        "../../../validation/dlpoly_config/vectors.CONFIG"
    ))
    .unwrap();
    let mut bad = original.clone();
    bad.atoms[0].force = None;
    assert!(write_dlpoly_config(&bad).is_err());
    let mut bad = original.clone();
    bad.levcfg = 0;
    assert!(write_dlpoly_config(&bad).is_err());
    let mut bad = original.clone();
    bad.declared_atom_count = Some(1);
    assert!(write_dlpoly_config(&bad).is_err());
    let mut bad = original.clone();
    bad.atoms[1].index = Some(1);
    assert!(write_dlpoly_config(&bad).is_err());
    let mut bad = original.clone();
    bad.atoms[0].header_extra = vec!["\n0 0".into()];
    assert!(write_dlpoly_config(&bad).is_err());
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut bad = original.clone();
        bad.atoms[0].velocity = Some([value; 3]);
        assert!(write_dlpoly_config(&bad).is_err());
        bad = original.clone();
        bad.atoms[0].force = Some([value; 3]);
        assert!(write_dlpoly_config(&bad).is_err());
    }
    let mut json = serde_json::to_value(original).unwrap();
    json["atoms"][0]["charge"] = 1.into();
    assert!(dlpoly_config_from_json(&json.to_string()).is_err());
    json["atoms"][0].as_object_mut().unwrap().remove("charge");
    json["units"] = "angstrom".into();
    assert!(dlpoly_config_from_json(&json.to_string()).is_err());
}

#[test]
fn byte_atom_and_fixed_record_limits_fail_before_unbounded_growth() {
    let input = include_str!("../../../validation/dlpoly_config/vectors.CONFIG");
    assert!(
        parse_dlpoly_config_with_limits(
            input,
            &DlpolyConfigLimits {
                max_input_bytes: input.len() - 1,
                max_atoms: 2
            }
        )
        .is_err()
    );
    assert!(
        parse_dlpoly_config_with_limits(
            input,
            &DlpolyConfigLimits {
                max_input_bytes: input.len(),
                max_atoms: 1
            }
        )
        .is_err()
    );
    assert!(
        parse_dlpoly_config_with_limits(
            input,
            &DlpolyConfigLimits {
                max_input_bytes: input.len(),
                max_atoms: 2
            }
        )
        .is_ok()
    );
    assert!(parse_dlpoly_config("title\n0 0 18446744073709551615\n").is_err());
    assert!(
        parse_dlpoly_config("title\n0 0\nX\n00000000000000000000000000000000000000000 0 0\n")
            .is_err()
    );
}

#[test]
fn numeric_output_is_exact_or_explicitly_refused_at_the_format_line_limit() {
    let mut record = parse_dlpoly_config("title\n0 0\nX\n0 0 0\n").unwrap();
    for n in [
        f64::MAX,
        -f64::MAX,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        -0.0,
    ] {
        record.atoms[0].position = [n, 0.0, 0.0];
        let back = parse_dlpoly_config(&write_dlpoly_config(&record).unwrap()).unwrap();
        assert_eq!(back.atoms[0].position[0].to_bits(), n.to_bits());
    }
    record.atoms[0].position = [-1.2345678901234567e-300; 3];
    assert!(
        write_dlpoly_config(&record)
            .unwrap_err()
            .detail
            .contains("precision loss")
    );
}

#[test]
fn empty_title_and_crlf_are_kept_without_accepting_blank_body_records() {
    let record = parse_dlpoly_config("\r\n0 0\r\nX\r\n0 0 0\r\n").unwrap();
    assert_eq!(record.title, "");
    assert!(parse_dlpoly_config("title\n0 0\n\nX\n0 0 0\n").is_err());
}
