use chematic_mol::{
    ConversionLimits, ConversionPolicy, Mol2ParseLimits, PreservationStatus, convert_with_report,
    project_mol2_to_cjson, project_mol2_to_cjson_with_limits,
};
use serde_json::Value;
const PACKET: &str = include_str!("../../../validation/mol2_cjson_common/cases.json");
fn cases() -> Vec<Value> {
    serde_json::from_str::<Value>(PACKET).unwrap()["cases"]
        .as_array()
        .unwrap()
        .clone()
}
fn same_bits(actual: &[Value], expected: impl Iterator<Item = String>) {
    let wanted: Vec<_> = expected
        .map(|v| v.parse::<f64>().unwrap().to_bits())
        .collect();
    assert_eq!(
        actual
            .iter()
            .map(|v| v.as_f64().unwrap().to_bits())
            .collect::<Vec<_>>(),
        wanted
    );
}
#[test]
fn selected_fields_survive_and_ambiguous_or_invalid_data_fail() {
    for case in cases() {
        let source = case["input"].as_str().unwrap();
        for input in [
            source.to_owned(),
            source.replace('\n', "\r\n"),
            source.trim_end_matches('\n').to_owned(),
        ] {
            let result = project_mol2_to_cjson(&input);
            if case["outcome"] != "preserved" {
                assert!(result.is_err(), "{}", case["id"]);
                continue;
            }
            let value: Value = serde_json::from_str(&result.unwrap()).unwrap();
            let expected = &case["expected"];
            assert_eq!(value["properties"]["name"], expected["name"]);
            assert_eq!(value["atoms"]["elements"]["number"], expected["numbers"]);
            let charges = value["atoms"]["formalCharges"]
                .as_array()
                .cloned()
                .unwrap_or(vec![Value::from(0); 2]);
            assert_eq!(Value::Array(charges), expected["formal_charges"]);
            same_bits(
                value["atoms"]["coords"]["3d"].as_array().unwrap(),
                expected["coordinates"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .flat_map(|r| r.as_array().unwrap())
                    .map(|v| v.as_str().unwrap().to_owned()),
            );
            assert_eq!(
                value["bonds"]["connections"]["index"],
                serde_json::json!([0, 1])
            );
            assert_eq!(value["bonds"]["order"], serde_json::json!([1]));
            if let Some(method) = expected["method"].as_str() {
                let sets = value["partialCharges"].as_object().unwrap();
                assert_eq!(sets.len(), 1);
                same_bits(
                    sets[method].as_array().unwrap(),
                    expected["partial_charges"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|v| v.as_str().unwrap().to_owned()),
                );
            } else {
                assert!(value.get("partialCharges").is_none());
            }
        }
    }
}
#[test]
fn report_keeps_loss_policy_and_marks_missing_declared_values_unverified() {
    for case in cases().iter().filter(|c| c["outcome"] == "refused") {
        let report = convert_with_report(
            case["input"].as_str().unwrap(),
            "mol2",
            "cjson",
            ConversionPolicy::AllowLoss,
        )
        .unwrap();
        assert!(report.output.is_some());
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == "partial_charges"
                    && matches!(
                        o.status,
                        PreservationStatus::Dropped | PreservationStatus::Unverified
                    ))
        );
        assert!(
            convert_with_report(
                case["input"].as_str().unwrap(),
                "mol2",
                "cjson",
                ConversionPolicy::Strict
            )
            .unwrap()
            .output
            .is_none()
        );
    }
    let packet = cases();
    let zero = packet
        .iter()
        .find(|c| c["id"] == "unassigned_zero_placeholders")
        .unwrap();
    let report = convert_with_report(
        zero["input"].as_str().unwrap(),
        "mol2",
        "cjson",
        ConversionPolicy::AllowLoss,
    )
    .unwrap();
    assert!(
        report.observations.iter().any(
            |o| o.field == "unassigned_zero_columns" && o.status == PreservationStatus::Dropped
        )
    );
}
#[test]
fn custom_limits_are_inclusive_for_input_atoms_bonds_and_output() {
    let packet = cases();
    let input = packet[0]["input"].as_str().unwrap();
    let output = project_mol2_to_cjson(input).unwrap();
    let exact = ConversionLimits {
        mol2: Mol2ParseLimits {
            max_input_bytes: input.len(),
            max_atoms: 2,
            max_bonds: 1,
            ..Default::default()
        },
        max_output_bytes: output.len(),
    };
    assert_eq!(
        project_mol2_to_cjson_with_limits(input, &exact).unwrap(),
        output
    );
    for limits in [
        ConversionLimits {
            mol2: Mol2ParseLimits {
                max_input_bytes: input.len() - 1,
                ..exact.mol2
            },
            ..exact
        },
        ConversionLimits {
            mol2: Mol2ParseLimits {
                max_atoms: 1,
                ..exact.mol2
            },
            ..exact
        },
        ConversionLimits {
            mol2: Mol2ParseLimits {
                max_bonds: 0,
                ..exact.mol2
            },
            ..exact
        },
        ConversionLimits {
            max_output_bytes: output.len() - 1,
            ..exact
        },
    ] {
        assert!(project_mol2_to_cjson_with_limits(input, &limits).is_err());
    }
}
