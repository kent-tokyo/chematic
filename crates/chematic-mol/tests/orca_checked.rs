use chematic_core::Element;
use chematic_mol::{
    OrcaAtom, OrcaCoords, OrcaInput, OrcaWriteLimits, parse_orca_input, write_orca_input_checked,
    write_orca_input_checked_with_limits, write_orca_input_json_checked,
    write_orca_input_json_checked_with_limits,
};

#[test]
fn authored_packet_preserves_fields_or_refuses_without_output() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/orca-checked-2026-10-24/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let result = write_orca_input_json_checked(case["input_json"].as_str().unwrap());
        if case["expected"] == "refuse" {
            assert!(result.is_err(), "{}", case["id"]);
            continue;
        }
        let text = result.unwrap_or_else(|e| panic!("{}: {e}", case["id"]));
        let input = parse_orca_input(&text).unwrap();
        assert_eq!(
            serde_json::json!(input.comments),
            case["fields"]["comments"]
        );
        assert_eq!(
            serde_json::json!(input.keywords),
            case["fields"]["keywords"]
        );
        let blocks: Vec<_> = input
            .blocks
            .iter()
            .map(|b| {
                serde_json::json!({
                    "name": b.name, "raw":b.raw, "has_end":b.has_end
                })
            })
            .collect();
        assert_eq!(serde_json::json!(blocks), case["fields"]["blocks"]);
        let Some(OrcaCoords::Xyz {
            charge,
            multiplicity,
            atoms,
        }) = &input.coords
        else {
            panic!("xyz missing");
        };
        let expected = &case["fields"]["coords"];
        assert_eq!(*charge, expected["charge"].as_i64().unwrap() as i32);
        assert_eq!(
            *multiplicity,
            expected["multiplicity"].as_u64().unwrap() as u32
        );
        assert_eq!(atoms.len(), expected["atoms"].as_array().unwrap().len());
        for (i, a) in atoms.iter().enumerate() {
            assert_eq!(a.element.symbol(), expected["atoms"][i]["element"]);
            assert_eq!(serde_json::json!(a.frozen), expected["atoms"][i]["frozen"]);
            assert_eq!(serde_json::json!(a.extra), expected["atoms"][i]["extra"]);
            for (axis, v) in [a.x, a.y, a.z].into_iter().enumerate() {
                assert_eq!(
                    format!("{:016x}", v.to_bits()),
                    case["coordinate_bits"][i][axis]
                );
            }
        }
        assert_eq!(write_orca_input_checked(&input).unwrap(), text);
    }
}

fn input(value: f64) -> OrcaInput {
    OrcaInput {
        keywords: vec!["Bohrs".into()],
        coords: Some(OrcaCoords::Xyz {
            charge: -1,
            multiplicity: 2,
            atoms: vec![OrcaAtom {
                element: Element::C,
                x: value,
                y: -value,
                z: -0.0,
                frozen: [true, false, true],
                extra: Some("NewGTO \"opaque\" end".into()),
            }],
        }),
        ..Default::default()
    }
}

#[test]
fn finite_binary64_pattern_trials_keep_bits_and_nonfinite_values_refuse() {
    let mut bits = 0x1234_5678_9abc_def0u64;
    for _ in 0..4096 {
        bits = bits
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let value = f64::from_bits(bits);
        let result = write_orca_input_checked(&input(value));
        if !value.is_finite() {
            assert!(result.is_err());
            continue;
        }
        let reparsed = parse_orca_input(&result.unwrap()).unwrap();
        let Some(OrcaCoords::Xyz { atoms, .. }) = reparsed.coords else {
            panic!("missing xyz");
        };
        assert_eq!(atoms[0].x.to_bits(), value.to_bits());
        assert_eq!(atoms[0].y.to_bits(), (-value).to_bits());
        assert_eq!(atoms[0].z.to_bits(), (-0.0f64).to_bits());
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(write_orca_input_checked(&input(value)).is_err());
    }
}

#[test]
fn atom_output_and_transport_byte_budgets_cover_exact_boundaries() {
    let data = input(1.2345678901234567);
    let text = write_orca_input_checked(&data).unwrap();
    let limits = OrcaWriteLimits {
        max_output_bytes: text.len(),
        max_atoms: 1,
        ..Default::default()
    };
    assert_eq!(
        write_orca_input_checked_with_limits(&data, &limits).unwrap(),
        text
    );
    for limits in [
        OrcaWriteLimits {
            max_output_bytes: text.len() - 1,
            ..limits
        },
        OrcaWriteLimits {
            max_atoms: 0,
            ..limits
        },
    ] {
        assert!(write_orca_input_checked_with_limits(&data, &limits).is_err());
    }
    let raw = r#"{"comments":[],"keywords":[],"blocks":[],"coords":{"type":"xyz","charge":0,"multiplicity":1,"atoms":[{"element":"He","x":-0.0,"y":1,"z":2,"frozen":[false,false,false]}]}}"#;
    let exact = OrcaWriteLimits {
        max_input_bytes: raw.len(),
        ..Default::default()
    };
    assert!(write_orca_input_json_checked_with_limits(raw, &exact).is_ok());
    assert!(
        write_orca_input_json_checked_with_limits(
            raw,
            &OrcaWriteLimits {
                max_input_bytes: raw.len() - 1,
                ..exact
            }
        )
        .is_err()
    );
}
