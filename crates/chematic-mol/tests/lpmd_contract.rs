use chematic_mol::*;
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/lpmd_preservation/cases.json"
    ))
    .unwrap()
}
#[test]
fn all_frames_columns_and_decimal_tokens_survive_three_cycles() {
    for case in packet()["valid"].as_array().unwrap() {
        let expected: LpmdDocument = serde_json::from_value(case["record"].clone()).unwrap();
        let mut record = parse_lpmd(case["text"].as_str().unwrap()).unwrap();
        for _ in 0..3 {
            assert_eq!(record, expected, "{}", case["id"]);
            record = parse_lpmd(&write_lpmd(&record).unwrap()).unwrap();
        }
        assert_eq!(record, expected);
    }
    for case in packet()["refused"].as_array().unwrap() {
        let err =
            parse_lpmd(case["text"].as_str().unwrap()).expect_err(case["id"].as_str().unwrap());
        assert!(err.to_string().starts_with("LPMD line "));
    }
}
#[test]
fn projection_uses_all_oriented_vectors_without_wrapping_and_keeps_unsafe_raw_records() {
    let p = packet();
    let tilted = p["valid"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "triclinic")
        .unwrap();
    let record = parse_lpmd(tilted["text"].as_str().unwrap()).unwrap();
    assert_eq!(
        record.positions_angstrom().unwrap()[0][0],
        [0.75, 2.125, 8.5]
    );
    for case in p["projection_refused"].as_array().unwrap() {
        let record = parse_lpmd(case["text"].as_str().unwrap()).unwrap();
        assert!(record.positions_angstrom().is_err());
        assert_eq!(parse_lpmd(&write_lpmd(&record).unwrap()).unwrap(), record);
    }
}
#[test]
fn edited_fields_json_strictness_and_whole_document_limits() {
    let p = packet();
    let text = p["valid"][0]["text"].as_str().unwrap();
    let mut record = parse_lpmd(text).unwrap();
    record.frames[0].atoms[0][1] = ".5".into();
    record.frames[0].cell[0][0] = "10".into();
    assert_eq!(
        parse_lpmd(&write_lpmd(&record).unwrap())
            .unwrap()
            .positions_angstrom()
            .unwrap()[0][0][0],
        5.
    );
    let json = serde_json::to_string(&record).unwrap();
    assert!(
        parse_lpmd_json(
            &json.replacen("\"columns\":", "\"columns\":[],\"columns\":", 1),
            &Default::default()
        )
        .is_err()
    );
    assert!(
        parse_lpmd_json(
            &json.replacen("\"cell\":", "\"c\\u0065ll\":[],\"cell\":", 1),
            &Default::default()
        )
        .is_err()
    );
    assert!(
        parse_lpmd_json(
            &json.replacen('{', "{\"unknown\":0,", 1),
            &Default::default()
        )
        .is_err()
    );
    for limits in [
        LpmdParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        LpmdParseLimits {
            max_frames: 0,
            ..Default::default()
        },
        LpmdParseLimits {
            max_line_len: 10,
            ..Default::default()
        },
        LpmdParseLimits {
            max_input_bytes: 10,
            ..Default::default()
        },
    ] {
        assert!(parse_lpmd_with_limits(text, &limits).is_err());
        assert!(write_lpmd_with_limits(&record, &limits).is_err());
    }
    let multi = p["valid"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "three-frames-variable-count")
        .unwrap();
    assert!(
        parse_lpmd_with_limits(
            multi["text"].as_str().unwrap(),
            &LpmdParseLimits {
                max_atoms: 4,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert!(
        parse_lpmd_with_limits(
            multi["text"].as_str().unwrap(),
            &LpmdParseLimits {
                max_frames: 2,
                ..Default::default()
            }
        )
        .is_err()
    );
}
