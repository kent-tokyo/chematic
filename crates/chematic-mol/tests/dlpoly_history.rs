use chematic_mol::*;
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/dlpoly_history/cases.json"
    ))
    .unwrap()
}
#[test]
fn all_original_frames_and_fields_round_trip() {
    for c in packet()["valid"].as_array().unwrap() {
        let g = parse_dlpoly_history(c["input"].as_str().unwrap())
            .unwrap_or_else(|e| panic!("{}: {e}", c["id"]));
        assert_eq!(
            serde_json::to_value(&g).unwrap(),
            c["expected"],
            "{}",
            c["id"]
        );
        assert_eq!(
            parse_dlpoly_history(&write_dlpoly_history(&g).unwrap()).unwrap(),
            g,
            "{}",
            c["id"]
        );
    }
}
#[test]
fn unsupported_and_truncated_records_refuse() {
    for c in packet()["invalid"].as_array().unwrap() {
        assert!(
            parse_dlpoly_history(c["input"].as_str().unwrap()).is_err(),
            "{}",
            c["id"]
        );
    }
}
#[test]
fn edits_and_strict_json_do_not_replay_old_source() {
    let mut g = parse_dlpoly_history(include_str!(
        "../../../validation/dlpoly_history/two-frames.HISTORY"
    ))
    .unwrap();
    g.frames[1].atoms[0].charge = -0.7654321;
    g.frames[1].atoms[0].force = Some([3.0, 4.0, 5.0]);
    assert_eq!(
        parse_dlpoly_history(&write_dlpoly_history(&g).unwrap()).unwrap(),
        g
    );
    let value = serde_json::to_value(&g).unwrap();
    for path in ["root", "frame", "atom"] {
        let mut json = value.clone();
        match path {
            "root" => json["unknown"] = 1.into(),
            "frame" => json["frames"][0]["time"] = 1.into(),
            _ => json["frames"][0]["atoms"][0]["formal_charge"] = 1.into(),
        }
        assert!(dlpoly_history_from_json(&json.to_string()).is_err());
    }
}
#[test]
fn limits_precision_and_level_contracts_refuse() {
    let s = include_str!("../../../validation/dlpoly_history/two-frames.HISTORY");
    for limits in [
        DlpolyHistoryLimits {
            max_atoms: 1,
            ..Default::default()
        },
        DlpolyHistoryLimits {
            max_frames: 1,
            ..Default::default()
        },
        DlpolyHistoryLimits {
            max_atom_records: 3,
            ..Default::default()
        },
        DlpolyHistoryLimits {
            max_input_bytes: s.len() - 1,
            ..Default::default()
        },
    ] {
        assert!(parse_dlpoly_history_with_limits(s, &limits).is_err());
    }
    let mut g = parse_dlpoly_history(s).unwrap();
    g.frames[0].atoms[0].position[0] = 0.123456789012345;
    assert!(write_dlpoly_history(&g).is_err());
    g.frames[0].atoms[0].position[0] = f64::NAN;
    assert!(write_dlpoly_history(&g).is_err());
    g.frames[0].atoms[0].position[0] = 0.0;
    g.frames[0].atoms[0].velocity = None;
    assert!(write_dlpoly_history(&g).is_err());
}
#[test]
fn fixed_record_lengths_and_signed_zero() {
    let mut g = parse_dlpoly_history(include_str!(
        "../../../validation/dlpoly_history/two-frames.HISTORY"
    ))
    .unwrap();
    g.frames[0].atoms[0].position[0] = -0.0;
    let text = write_dlpoly_history(&g).unwrap();
    let rows: Vec<_> = text.lines().collect();
    assert_eq!(rows[1].len(), 30);
    assert_eq!(rows[2].len(), 60);
    assert_eq!(rows[3].len(), 42);
    assert_eq!(rows[4].len(), 36);
    assert_eq!(
        parse_dlpoly_history(&text).unwrap().frames[0].atoms[0].position[0].to_bits(),
        (-0.0_f64).to_bits()
    );
}
