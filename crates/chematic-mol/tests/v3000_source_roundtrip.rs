use chematic_mol::{
    MolParseError, V3000SourceLimits, roundtrip_mol_v3000_source,
    roundtrip_mol_v3000_source_with_limits,
};
use serde_json::Value;
const PACKET: &str = include_str!("../../../validation/v3000_source_roundtrip/cases.json");
fn first() -> String {
    serde_json::from_str::<Value>(PACKET).unwrap()["cases"][0]["input"]
        .as_str()
        .unwrap()
        .to_owned()
}

#[test]
fn retains_all_source_fields_in_three_line_endings() {
    let packet: Value = serde_json::from_str(PACKET).unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let source = case["input"].as_str().unwrap();
        for input in [
            source.to_owned(),
            source.replace('\n', "\r\n"),
            source.trim_end_matches('\n').to_owned(),
        ] {
            assert_eq!(
                roundtrip_mol_v3000_source(&input).unwrap(),
                input,
                "{}",
                case["id"]
            );
        }
    }
}

#[test]
fn rejects_ambiguous_or_incomplete_source_instead_of_retaining_it() {
    let input = first();
    let extra = input.replace(
        "M  V30 END CTAB",
        "M  V30 END CTAB\nM  V30 BEGIN CTAB\nM  V30 END CTAB",
    );
    for bad in [
        input.replace("M  END", ""),
        input.clone() + "$$$$\n",
        input.clone() + &input,
        extra,
        input.replace("END CTAB", "END CTAB-"),
        input.replace("M  V30 2 O", "M  V30 1 O"),
        input.replace("M  V30 1 C", "M  V30 0 C"),
        input.replace("M  V30 1 1 1 2", "M  V30 1 1 1 99"),
        input.replace("0.123456789012345", "NaN"),
        input.replace("M  V30 END CTAB", "M  V30 END CTAB\nM  V30 unvalidated"),
    ] {
        assert!(roundtrip_mol_v3000_source(&bad).is_err(), "{bad}");
    }
}

#[test]
fn resource_limits_are_inclusive_and_fail_without_output() {
    let input = first();
    let defaults = V3000SourceLimits::default();
    let exact = V3000SourceLimits {
        max_input_bytes: input.len(),
        max_atoms: 2,
        max_lines: input.lines().count(),
        max_line_bytes: input.lines().map(str::len).max().unwrap(),
    };
    assert_eq!(
        roundtrip_mol_v3000_source_with_limits(&input, &exact).unwrap(),
        input
    );
    for limits in [
        V3000SourceLimits {
            max_input_bytes: input.len() - 1,
            ..defaults
        },
        V3000SourceLimits {
            max_atoms: 1,
            ..defaults
        },
        V3000SourceLimits {
            max_lines: input.lines().count() - 1,
            ..defaults
        },
        V3000SourceLimits {
            max_line_bytes: 10,
            ..defaults
        },
    ] {
        assert!(matches!(
            roundtrip_mol_v3000_source_with_limits(&input, &limits),
            Err(MolParseError::ResourceLimit { .. })
        ));
    }
}
