use chematic_mol::{
    OrcaInputError, OrcaInputParseLimits, parse_orca_cartesian, parse_orca_cartesian_with_limits,
    parse_orca_input, parse_orca_input_with_limits,
};

#[test]
fn unit_aware_projection_matches_all_authored_fields() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/orca-cartesian-2026-10-22/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let source = case["source"].as_str().unwrap();
        let result = parse_orca_cartesian(source);
        if case["expected"] == "refuse" {
            assert!(result.is_err(), "{}", case["id"]);
            if case["core_refuse"] == true {
                assert!(parse_orca_input(source).is_err(), "{}", case["id"]);
            }
            continue;
        }
        let geometry = result.unwrap();
        let json = geometry.to_json();
        for key in [
            "elements",
            "frozen",
            "atom_extras",
            "charge",
            "multiplicity",
            "coordinate_units",
            "units_explicit",
        ] {
            assert_eq!(json[key], case["fields"][key], "{} {key}", case["id"]);
        }
        for (row, atom) in geometry.atoms.iter().enumerate() {
            for (axis, value) in [atom.x, atom.y, atom.z].iter().enumerate() {
                assert_eq!(
                    format!("{:016x}", value.to_bits()),
                    case["fields"]["coordinate_bits"][row][axis]
                        .as_str()
                        .unwrap()
                );
            }
        }
    }
}

#[test]
fn atom_limit_is_enforced_before_parsing_the_next_row() {
    let source = "* xyz 0 1\nHe 1 2 3\nnot a valid atom\n*\n";
    let limits = OrcaInputParseLimits {
        max_atoms: 1,
        ..Default::default()
    };
    assert!(matches!(
        parse_orca_input_with_limits(source, &limits),
        Err(OrcaInputError::ResourceLimit {
            resource: "atoms",
            actual: 2,
            limit: 1
        })
    ));
}

#[test]
fn parser_limits_cover_checked_projection_and_initial_bytes() {
    let source = "! Bohrs\n* xyz 0 1\nHe 2.0 0 0\n*\n";
    let exact = OrcaInputParseLimits {
        max_input_bytes: source.len(),
        max_atoms: 1,
        ..Default::default()
    };
    assert!(parse_orca_cartesian_with_limits(source, &exact).is_ok());
    for limits in [
        OrcaInputParseLimits {
            max_input_bytes: source.len() - 1,
            ..exact
        },
        OrcaInputParseLimits {
            max_atoms: 0,
            ..exact
        },
        OrcaInputParseLimits {
            max_lines: 2,
            ..exact
        },
        OrcaInputParseLimits {
            max_line_bytes: 4,
            ..exact
        },
        OrcaInputParseLimits {
            max_keywords: 0,
            ..exact
        },
    ] {
        assert!(parse_orca_cartesian_with_limits(source, &limits).is_err());
    }
}
