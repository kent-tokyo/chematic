use chematic_mol::*;
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/nwchem_geometry/cases.json"
    ))
    .unwrap()
}
#[test]
fn original_cases_preserve_every_structured_field() {
    for c in packet()["valid"].as_array().unwrap() {
        let g = parse_nwchem_geometry(c["input"].as_str().unwrap())
            .unwrap_or_else(|e| panic!("{}: {e}", c["id"]));
        assert_eq!(
            serde_json::to_value(&g).unwrap(),
            c["expected"],
            "{}",
            c["id"]
        );
        assert_eq!(
            parse_nwchem_geometry(&write_nwchem_geometry(&g).unwrap()).unwrap(),
            g
        );
    }
}
#[test]
fn unsupported_and_ambiguous_inputs_refuse() {
    for c in packet()["invalid"].as_array().unwrap() {
        assert!(
            parse_nwchem_geometry(c["input"].as_str().unwrap()).is_err(),
            "{}",
            c["id"]
        );
    }
}
#[test]
fn edits_and_unknown_json_fields() {
    let mut g = parse_nwchem_geometry(include_str!(
        "../../../validation/nwchem_geometry/centers.nw"
    ))
    .unwrap();
    g.centers[0].nuclear_charge = Some(-0.12345678901234566);
    g.centers[0].mass_amu = Some(18.01010101010101);
    g.units = Some("pm".into());
    assert_eq!(
        parse_nwchem_geometry(&write_nwchem_geometry(&g).unwrap()).unwrap(),
        g
    );
    let mut json = serde_json::to_value(&g).unwrap();
    json["centers"][0]["formal_charge"] = 1.into();
    assert!(nwchem_geometry_from_json(&json.to_string()).is_err());
    g.centers[0].mass_amu = Some(0.0);
    assert!(write_nwchem_geometry(&g).is_err());
}
#[test]
fn limits_and_invalid_edited_records_refuse() {
    let s = include_str!("../../../validation/nwchem_geometry/centers.nw");
    for limits in [
        NwchemLimits {
            max_centers: 0,
            ..Default::default()
        },
        NwchemLimits {
            max_line_len: 20,
            ..Default::default()
        },
        NwchemLimits {
            max_input_bytes: s.len() - 1,
            ..Default::default()
        },
    ] {
        assert!(parse_nwchem_geometry_with_limits(s, &limits).is_err());
    }
    let mut g = parse_nwchem_geometry(s).unwrap();
    g.autosym = Some(false);
    assert!(write_nwchem_geometry(&g).is_err());
    g.autosym = Some(true);
    g.centers[0].velocity_au = Some([f64::INFINITY, 0.0, 0.0]);
    assert!(write_nwchem_geometry(&g).is_err());
}
#[test]
fn extreme_finite_values_retain_bits_and_dummy_tags_stay_tags() {
    let mut g = parse_nwchem_geometry("geometry\nbqH 0 0 0 charge -0.5\nend\n").unwrap();
    for v in [f64::from_bits(1), f64::MIN_POSITIVE, f64::MAX, -0.0] {
        g.centers[0].position[0] = v;
        let got = parse_nwchem_geometry(&write_nwchem_geometry(&g).unwrap()).unwrap();
        assert_eq!(got.centers[0].position[0].to_bits(), v.to_bits());
        assert_eq!(got.centers[0].tag, "bqH");
        assert_eq!(got.centers[0].nuclear_charge, Some(-0.5));
    }
}
