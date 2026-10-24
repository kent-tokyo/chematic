use chematic_mol::{
    TinkerParseLimits, parse_tinker_xyz, parse_tinker_xyz_with_limits, tinker_xyz_from_json,
    write_tinker_xyz,
};
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!("../../../validation/tinker_xyz/cases.json")).unwrap()
}

#[test]
fn authored_records_parse_and_roundtrip_without_perception() {
    for case in packet()["cases"].as_array().unwrap() {
        let record = parse_tinker_xyz(case["input"].as_str().unwrap()).unwrap();
        assert_eq!(
            serde_json::to_value(&record).unwrap(),
            case["expected"],
            "{}",
            case["id"]
        );
        let written = write_tinker_xyz(&record).unwrap();
        assert_eq!(
            parse_tinker_xyz(&written).unwrap(),
            record,
            "{}",
            case["id"]
        );
        assert_eq!(
            tinker_xyz_from_json(&serde_json::to_string(&record).unwrap()).unwrap(),
            record
        );
    }
}

#[test]
fn unsupported_or_invalid_input_is_refused() {
    for case in packet()["refusals"].as_array().unwrap() {
        assert!(
            parse_tinker_xyz(case["input"].as_str().unwrap()).is_err(),
            "{}",
            case["id"]
        );
    }
}

#[test]
fn writes_current_fields_including_custom_names_types_cell_and_connections() {
    let mut record =
        parse_tinker_xyz(include_str!("../../../validation/tinker_xyz/water.txyz")).unwrap();
    record.title = "edited".into();
    record.atoms[0].name = "OW".into();
    record.atoms[0].atom_type = 90001;
    record.atoms[0].position_angstrom[0] = 1.2345678901234567;
    record.atoms[0].neighbors = vec![2];
    record.atoms[2].neighbors.clear();
    record.cell = Some(chematic_mol::TinkerCell {
        lengths_angstrom: [12.345678901234567, 20.0, 30.0],
        angles_degrees: [70.0, 80.0, 90.0],
    });
    assert_eq!(
        parse_tinker_xyz(&write_tinker_xyz(&record).unwrap()).unwrap(),
        record
    );
}

#[test]
fn constructed_data_and_json_cannot_bypass_validation() {
    let original =
        parse_tinker_xyz(include_str!("../../../validation/tinker_xyz/water.txyz")).unwrap();
    for invalid in ["title\ninjected", " leading"] {
        let mut r = original.clone();
        r.title = invalid.into();
        assert!(write_tinker_xyz(&r).is_err());
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut r = original.clone();
        r.atoms[0].position_angstrom[0] = value;
        assert!(write_tinker_xyz(&r).is_err());
    }
    for neighbor in [0, 1, 4, usize::MAX] {
        let mut r = original.clone();
        r.atoms[0].neighbors.push(neighbor);
        assert!(write_tinker_xyz(&r).is_err());
    }
    let mut json = serde_json::to_value(&original).unwrap();
    json["formal_charge"] = 1.into();
    assert!(tinker_xyz_from_json(&json.to_string()).is_err());
    json.as_object_mut().unwrap().remove("formal_charge");
    json["atoms"][0]["isotope"] = 18.into();
    assert!(tinker_xyz_from_json(&json.to_string()).is_err());
    let duplicate = "{\"title\":\"a\",\"title\":\"b\",\"atoms\":[],\"cell\":null}";
    assert!(tinker_xyz_from_json(duplicate).is_err());
}

#[test]
fn parser_limits_apply_before_count_allocation_and_neighbor_growth() {
    let input = include_str!("../../../validation/tinker_xyz/water.txyz");
    let defaults = TinkerParseLimits::default();
    for limits in [
        TinkerParseLimits {
            max_input_bytes: input.len() - 1,
            ..defaults
        },
        TinkerParseLimits {
            max_atoms: 2,
            ..defaults
        },
        TinkerParseLimits {
            max_neighbors: 1,
            ..defaults
        },
        TinkerParseLimits {
            max_line_len: 5,
            ..defaults
        },
    ] {
        assert!(parse_tinker_xyz_with_limits(input, &limits).is_err());
    }
    assert!(parse_tinker_xyz("18446744073709551615\n").is_err());
    assert!(
        parse_tinker_xyz_with_limits(
            input,
            &TinkerParseLimits {
                max_input_bytes: input.len(),
                max_atoms: 3,
                max_neighbors: 2,
                ..defaults
            }
        )
        .is_ok()
    );
}

#[test]
fn finite_f64_extremes_survive_text_roundtrip() {
    let mut r = parse_tinker_xyz("1\n1 X 0 0 0 0\n").unwrap();
    for x in [
        f64::MAX,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        -0.0,
        -f64::MAX,
    ] {
        r.atoms[0].position_angstrom = [x; 3];
        let back = parse_tinker_xyz(&write_tinker_xyz(&r).unwrap()).unwrap();
        assert_eq!(back.atoms[0].position_angstrom[0].to_bits(), x.to_bits());
    }
}

#[test]
fn crlf_and_tabs_in_atom_rows_are_accepted() {
    let input = "  1  title with spaces\r\n1\tX\t0\t0\t0\t0\r\n";
    assert_eq!(parse_tinker_xyz(input).unwrap().title, "title with spaces");
}
