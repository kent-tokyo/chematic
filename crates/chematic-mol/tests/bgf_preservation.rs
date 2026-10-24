use chematic_mol::*;
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!("../../../validation/bgf/cases.json")).unwrap()
}
#[test]
fn shared_bgf_packet_preserves_all_fields() {
    for case in packet()["valid"].as_array().unwrap() {
        let record = parse_bgf(case["input"].as_str().unwrap())
            .unwrap_or_else(|e| panic!("{}: {e}", case["id"]));
        assert_eq!(
            serde_json::to_value(&record).unwrap(),
            case["expected"],
            "{}",
            case["id"]
        );
        assert_eq!(parse_bgf(&write_bgf(&record).unwrap()).unwrap(), record);
    }
}
#[test]
fn shared_bgf_refusals_are_not_silently_projected() {
    for case in packet()["invalid"].as_array().unwrap() {
        assert!(
            parse_bgf(case["input"].as_str().unwrap()).is_err(),
            "{}",
            case["id"]
        );
    }
}
#[test]
fn edited_record_and_strict_json() {
    let mut r = parse_bgf(include_str!("../../../validation/bgf/carbonyl.bgf")).unwrap();
    r.atoms[0].partial_charge = -0.98765;
    r.atoms[0].lone_pairs = 3;
    assert_eq!(parse_bgf(&write_bgf(&r).unwrap()).unwrap(), r);
    let mut json = serde_json::to_value(&r).unwrap();
    json["atoms"][0]["formal_charge"] = 1.into();
    assert!(bgf_record_from_json(&json.to_string()).is_err());
    r.atoms[0].partial_charge = 0.123456;
    assert!(write_bgf(&r).is_err());
    r.atoms[0].partial_charge = f64::NAN;
    assert!(write_bgf(&r).is_err());
}
#[test]
fn field_and_resource_limits() {
    let s = include_str!("../../../validation/bgf/carbonyl.bgf");
    for limits in [
        BgfLimits {
            max_atoms: 1,
            ..Default::default()
        },
        BgfLimits {
            max_input_bytes: s.len() - 1,
            ..Default::default()
        },
        BgfLimits {
            max_line_len: 30,
            ..Default::default()
        },
    ] {
        assert!(parse_bgf_with_limits(s, &limits).is_err());
    }
    let mut r = parse_bgf(s).unwrap();
    r.atoms[0].name = "TOO-LONG".into();
    assert!(write_bgf(&r).is_err());
    r.atoms[0].name = "C".into();
    r.atoms[0].position_angstrom[0] = 10000.0;
    assert!(write_bgf(&r).is_err());
    r.atoms[0].position_angstrom[0] = 1.0;
    r.connectivity[0].orders = Some(vec![1]);
    assert!(write_bgf(&r).is_err());
}

#[test]
fn eleven_neighbors_fit_and_a_twelfth_refuses() {
    let mut r = parse_bgf(include_str!("../../../validation/bgf/carbonyl.bgf")).unwrap();
    let template = r.atoms[0].clone();
    r.atoms = (1..=13)
        .map(|id| {
            let mut a = template.clone();
            a.id = id;
            a
        })
        .collect();
    r.connectivity = vec![BgfConnectivity {
        atom: 1,
        neighbors: (2..=12).collect(),
        orders: Some(vec![1; 11]),
    }];
    assert_eq!(parse_bgf(&write_bgf(&r).unwrap()).unwrap(), r);
    r.connectivity[0].neighbors.push(13);
    r.connectivity[0].orders.as_mut().unwrap().push(1);
    assert!(write_bgf(&r).is_err());
}
