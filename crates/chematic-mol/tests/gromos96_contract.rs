use chematic_mol::*;
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/gromos96_preservation/cases.json"
    ))
    .unwrap()
}
fn read(text: &str, unit: Gromos96LengthUnit) -> Gromos96Document {
    parse_gromos96_with_limits(text, unit, &Default::default()).unwrap()
}
#[test]
fn retains_all_blocks_current_payload_and_caller_owned_units() {
    for c in packet()["valid"].as_array().unwrap() {
        let expected: Gromos96Document = serde_json::from_value(c["record"].clone()).unwrap();
        let mut r = read(c["text"].as_str().unwrap(), expected.length_unit);
        for _ in 0..3 {
            assert_eq!(r, expected, "{}", c["id"]);
            r = read(&write_gromos96(&r).unwrap(), expected.length_unit);
        }
        assert_eq!(r, expected);
    }
    for c in packet()["refused"].as_array().unwrap() {
        let error =
            parse_gromos96(c["text"].as_str().unwrap()).expect_err(c["id"].as_str().unwrap());
        assert!(error.to_string().starts_with("GROMOS96 line "));
    }
}
#[test]
fn units_scale_positions_and_box_together_without_element_or_velocity_inference() {
    let p = packet();
    let text = p["valid"][0]["text"].as_str().unwrap();
    let unknown = parse_gromos96(text).unwrap();
    assert_eq!(unknown.length_unit, Gromos96LengthUnit::Unspecified);
    assert!(unknown.geometry_angstrom().is_err());
    let r = read(text, Gromos96LengthUnit::Nm);
    let g = r.geometry_angstrom().unwrap();
    assert_eq!(g[0].positions_angstrom[0], [1.25, 2.5, 5.]);
    assert_eq!(
        g[0].box_vectors_angstrom,
        Some([[20., 0., 0.], [0., 40., 0.], [0., 0., 80.]])
    );
    assert_eq!(g[0].identities.as_ref().unwrap()[1].atom_name, "Cl");
    let direct = read(text, Gromos96LengthUnit::Angstrom)
        .geometry_angstrom()
        .unwrap();
    assert_eq!(direct[0].positions_angstrom[0], [0.125, 0.25, 0.5]);
    for c in p["projection_refused"].as_array().unwrap() {
        let expected: Gromos96Document = serde_json::from_value(c["record"].clone()).unwrap();
        let r = read(c["text"].as_str().unwrap(), expected.length_unit);
        assert!(r.geometry_angstrom().is_err());
        assert_eq!(
            read(&write_gromos96(&r).unwrap(), expected.length_unit),
            expected
        );
    }
}
#[test]
fn edits_json_ambiguity_limits_and_aggregate_frame_bounds() {
    let p = packet();
    let text = p["valid"][0]["text"].as_str().unwrap();
    let mut r = read(text, Gromos96LengthUnit::Nm);
    r.blocks[1].lines[0].replace_range(24..39, &format!("{:>15}", ".5"));
    assert_eq!(
        read(&write_gromos96(&r).unwrap(), r.length_unit)
            .geometry_angstrom()
            .unwrap()[0]
            .positions_angstrom[0][0],
        5.
    );
    let json = serde_json::to_string(&r).unwrap();
    let limits = Gromos96ParseLimits::default();
    for bad in [
        json.replacen("\"blocks\":", "\"blocks\":[],\"blocks\":", 1),
        json.replacen("\"lines\":", "\"l\\u0069nes\":[],\"lines\":", 1),
        json.replacen('{', "{\"unknown\":1,", 1),
    ] {
        assert!(parse_gromos96_json(&bad, &limits).is_err());
    }
    for l in [
        Gromos96ParseLimits {
            max_frames: 0,
            ..limits
        },
        Gromos96ParseLimits {
            max_atoms: 1,
            ..limits
        },
        Gromos96ParseLimits {
            max_input_bytes: 10,
            ..limits
        },
        Gromos96ParseLimits {
            max_line_len: 10,
            ..limits
        },
    ] {
        assert!(parse_gromos96_with_limits(text, r.length_unit, &l).is_err());
        assert!(write_gromos96_with_limits(&r, &l).is_err());
    }
    let multi = p["valid"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "two-frames-full")
        .unwrap();
    assert!(
        parse_gromos96_with_limits(
            multi["text"].as_str().unwrap(),
            r.length_unit,
            &Gromos96ParseLimits {
                max_atoms: 3,
                ..limits
            }
        )
        .is_err()
    );
    assert!(
        parse_gromos96_with_limits(
            multi["text"].as_str().unwrap(),
            r.length_unit,
            &Gromos96ParseLimits {
                max_frames: 1,
                ..limits
            }
        )
        .is_err()
    );
}
