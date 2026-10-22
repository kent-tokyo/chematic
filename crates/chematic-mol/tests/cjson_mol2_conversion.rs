use chematic_core::AtomIdx;
use chematic_mol::{
    ConversionLimits, ConversionPolicy as Policy, PreservationStatus as Status,
    convert_with_report, convert_with_report_with_limits, parse_cjson_record,
    parse_mol2_source_record,
};
use serde_json::{Value, json};

const SOURCE: &str = include_str!("../../../validation/cjson_preservation/partial-charges.cjson");
fn source() -> Value {
    serde_json::from_str(SOURCE).unwrap()
}
fn status(report: &chematic_mol::ConversionReport, field: &str) -> Status {
    report
        .observations
        .iter()
        .find(|o| o.field == field)
        .unwrap()
        .status
}

#[test]
fn typed_metadata_has_exact_values_and_nested_opaque_paths() {
    let mut value = source();
    value["partialCharges"]["other/method~"] = json!([-0.0, -1.4175225513503917, 5e-324, f64::MAX]);
    let text = value.to_string();
    let record = parse_cjson_record(&text).unwrap();
    assert_eq!(record.source(), text);
    assert_eq!(record.name(), Some("acetamide charges"));
    assert_eq!(record.opaque_fields(), ["/properties/vendor"]);
    assert_eq!(record.partial_charge_datasets().len(), 2);
    let values = &record.partial_charge_datasets()["other/method~"];
    for (a, b) in values
        .iter()
        .zip([-0.0_f64, -1.4175225513503917, 5e-324, f64::MAX])
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    assert_eq!(
        record.to_json()["partial_charge_datasets"]["other/method~"]
            .as_array()
            .unwrap()
            .len(),
        4
    );
}

#[test]
fn malformed_known_metadata_is_rejected_but_legacy_reader_remains_compatible() {
    for bad in [
        Value::Null,
        json!([]),
        json!({"method":null}),
        json!({"method":[]}),
        json!({"method":[0,0,0]}),
        json!({"method":[0,0,0,0,0]}),
        json!({"method":[0,0,0,"NaN"]}),
        json!({"method":[0,0,0,null]}),
    ] {
        let mut value = source();
        value["partialCharges"] = bad;
        assert!(parse_cjson_record(&value.to_string()).is_err(), "{value}");
        assert!(chematic_mol::parse_cjson(&value.to_string()).is_ok());
    }
    for bad in [
        Value::Null,
        json!(1),
        json!([]),
        json!({"name":null}),
        json!({"name":1}),
    ] {
        let mut value = source();
        value["properties"] = bad;
        assert!(parse_cjson_record(&value.to_string()).is_err(), "{value}");
    }
    assert!(parse_cjson_record(&SOURCE.replace("0.7, 0.25", "1e400, 0.25")).is_err());
    assert!(
        parse_cjson_record(&SOURCE.replace(
            "\"USER_CHARGES\":",
            "\"USER_CHARGES\": [0,0,0,0], \"USER_CHARGES\":"
        ))
        .is_err()
    );
}

#[test]
fn mol2_candidate_preserves_graph_name_charge_method_and_all_float_bits() {
    for text in [SOURCE.to_owned(), SOURCE.replace('\n', "\r\n")] {
        let before = parse_cjson_record(&text).unwrap();
        let report = convert_with_report(&text, "cjson", "mol2", Policy::AllowLoss).unwrap();
        assert_eq!(report.output_format, "mol2");
        assert_eq!(report.status, "completed_with_loss");
        assert!(report.candidate_graph_verified);
        let after = parse_mol2_source_record(report.output.as_ref().unwrap()).unwrap();
        assert_eq!(after.name(), before.name());
        assert_eq!(after.charge_type(), Some("USER_CHARGES"));
        assert_eq!(report.candidate_atom_mapping.len(), 4);
        assert_eq!(report.candidate_bond_mapping.len(), 3);
        for map in &report.candidate_atom_mapping {
            assert_eq!(map.source_id as usize, map.source_index);
            let a = before.coordinates()[map.source_index];
            let b = after.coordinates()[map.candidate_index];
            assert_eq!(
                [a.0.to_bits(), a.1.to_bits(), a.2.to_bits()],
                [b.0.to_bits(), b.1.to_bits(), b.2.to_bits()]
            );
            assert_eq!(
                before.partial_charge_datasets()["USER_CHARGES"][map.source_index].to_bits(),
                after.atoms()[map.candidate_index]
                    .partial_charge
                    .unwrap()
                    .to_bits()
            );
        }
        for field in [
            "projected_graph",
            "coordinates",
            "partial_charges",
            "charge_method",
            "molecule_name",
        ] {
            assert_eq!(status(&report, field), Status::Preserved, "{field}");
        }
        assert_eq!(status(&report, "opaque_fields"), Status::Dropped);
        assert_eq!(status(&report, "mol2_metadata"), Status::Inferred);
        let strict = convert_with_report(&text, "cjson", "mol2", Policy::Strict).unwrap();
        assert_eq!(strict.status, "refused");
        assert!(strict.output.is_none());
        assert!(strict.candidate_graph_verified);
        assert!(convert_with_report(&text, "cjson", "mol2", Policy::PreserveSource).is_err());
    }
}

#[test]
fn missing_geometry_or_ambiguous_dataset_refuses_even_allow_loss() {
    let mut missing = source();
    missing["atoms"].as_object_mut().unwrap().remove("coords");
    let mut multiple = source();
    multiple["partialCharges"]["SECOND"] = json!([0, 0, 0, 0]);
    let mut empty = source();
    empty["atoms"]["elements"]["number"] = json!([]);
    empty["atoms"]["formalCharges"] = json!([]);
    empty["atoms"]["isotopes"] = json!([]);
    empty["atoms"]["coords"]["3d"] = json!([]);
    empty["bonds"]["connections"]["index"] = json!([]);
    empty["bonds"]["order"] = json!([]);
    empty["partialCharges"]["USER_CHARGES"] = json!([]);
    for value in [missing, multiple, empty] {
        let report =
            convert_with_report(&value.to_string(), "cjson", "mol2", Policy::AllowLoss).unwrap();
        assert_eq!(report.status, "refused");
        assert!(report.output.is_none());
        assert!(!report.candidate_graph_verified);
        assert!(report.candidate_atom_mapping.is_empty());
    }
}

#[test]
fn unrepresentable_names_and_method_tokens_never_inject_mol2_syntax() {
    for name in [
        "",
        " leading",
        "trailing ",
        "bad\n@<TRIPOS>ATOM",
        "a\rname",
        "@<tripos>ATOM",
        "a\0b",
    ] {
        let mut value = source();
        value["properties"]["name"] = json!(name);
        assert!(parse_cjson_record(&value.to_string()).is_ok());
        let report =
            convert_with_report(&value.to_string(), "cjson", "mol2", Policy::AllowLoss).unwrap();
        assert_eq!(report.status, "refused", "{name:?}");
        assert!(report.output.is_none());
    }
    for method in [
        "",
        "NO_CHARGES",
        "no_charges",
        "two words",
        "x\n@<TRIPOS>ATOM",
        "@<tripos>ATOM",
        "a\0b",
    ] {
        let mut value = source();
        value["partialCharges"] = json!({method:[0,0,0,0]});
        let report =
            convert_with_report(&value.to_string(), "cjson", "mol2", Policy::AllowLoss).unwrap();
        assert_eq!(report.status, "refused", "{method:?}");
        assert!(report.output.is_none());
    }
    let mut value = source();
    value["properties"]["name"] = json!("アセトアミド 研究用");
    let report =
        convert_with_report(&value.to_string(), "cjson", "mol2", Policy::AllowLoss).unwrap();
    assert_eq!(status(&report, "molecule_name"), Status::Preserved);
}

#[test]
fn no_partial_dataset_does_not_create_charge_columns() {
    for datasets in [None, Some(json!({}))] {
        let mut value = source();
        value.as_object_mut().unwrap().remove("partialCharges");
        if let Some(datasets) = datasets {
            value["partialCharges"] = datasets;
        }
        value.as_object_mut().unwrap().remove("properties");
        let report =
            convert_with_report(&value.to_string(), "cjson", "mol2", Policy::AllowLoss).unwrap();
        assert!(report.candidate_graph_verified);
        let after = parse_mol2_source_record(report.output.as_ref().unwrap()).unwrap();
        assert_eq!(after.charge_type(), Some("NO_CHARGES"));
        assert!(after.atoms().iter().all(|a| a.partial_charge.is_none()));
        assert!(
            !report
                .observations
                .iter()
                .any(|o| o.field == "partial_charges")
        );
        assert_eq!(status(&report, "molecule_name"), Status::Inferred);
    }
}

#[test]
fn formal_charges_use_unity_independently_of_partial_charges() {
    // Original synthetic representation tests (project license), not claims
    // of chemically valid charge states. Exercise every supported i8 value.
    for charge in i8::MIN..=i8::MAX {
        for number in [7, 8, 11] {
            let value = json!({"atoms":{"elements":{"number":[number]},
                "formalCharges":[charge], "coords":{"3d":[1,-0.0,5e-324]}},
                "partialCharges":{"USER_CHARGES":[-0.7]}});
            let report =
                convert_with_report(&value.to_string(), "cjson", "mol2", Policy::AllowLoss)
                    .unwrap();
            assert!(
                report.candidate_graph_verified,
                "element {number}, charge {charge}"
            );
            assert_eq!(status(&report, "formal_charges"), Status::Preserved);
            assert_eq!(status(&report, "partial_charges"), Status::Preserved);
            assert_eq!(status(&report, "coordinates"), Status::Preserved);
            let output = report.output.as_ref().unwrap();
            let after = parse_mol2_source_record(output).unwrap();
            assert_eq!(after.molecule().atom(AtomIdx(0)).charge, charge);
            assert_eq!(after.atoms()[0].partial_charge, Some(-0.7));
            assert_eq!(after.charge_type(), Some("USER_CHARGES"));
            assert_eq!(report.candidate_atom_mapping.len(), 1);
            assert_eq!(report.candidate_atom_mapping[0].source_id, 0);
            if charge != 0 {
                assert!(output.contains(&format!("1 1\ncharge {charge}\n")));
                // Open Babel 3.2.1 reads charges before BOND, but skips them
                // when the same valid UNITY section follows BOND.
                assert!(
                    output.find("@<TRIPOS>UNITY_ATOM_ATTR").unwrap()
                        < output.find("@<TRIPOS>BOND").unwrap()
                );
            } else {
                assert!(!output.contains("@<TRIPOS>UNITY_ATOM_ATTR"));
            }
        }
    }
}

#[test]
fn unity_formal_charges_cover_sparse_ids_without_partial_charge_columns() {
    let value = json!({"atoms":{"elements":{"number":[11,6,17]},
        "formalCharges":[1,0,-1],"coords":{"3d":[0,0,0,1,0,0,2,0,0]}}});
    let source = value.to_string();
    let report = convert_with_report(&source, "cjson", "mol2", Policy::AllowLoss).unwrap();
    assert!(report.candidate_graph_verified);
    assert_eq!(status(&report, "formal_charges"), Status::Preserved);
    let output = report.output.as_ref().unwrap();
    assert!(output.contains("@<TRIPOS>UNITY_ATOM_ATTR\n1 1\ncharge 1\n3 1\ncharge -1\n"));
    let after = parse_mol2_source_record(output).unwrap();
    assert_eq!(after.charge_type(), Some("NO_CHARGES"));
    assert!(after.atoms().iter().all(|a| a.partial_charge.is_none()));
    assert_eq!(report.candidate_atom_mapping.len(), 3);
    let strict = convert_with_report(&source, "cjson", "mol2", Policy::Strict).unwrap();
    assert_eq!(strict.status, "refused");
    assert!(strict.output.is_none());
    assert!(strict.candidate_graph_verified);
    assert_eq!(status(&strict, "formal_charges"), Status::Preserved);

    let exact = ConversionLimits {
        max_output_bytes: output.len(),
        ..Default::default()
    };
    assert!(
        convert_with_report_with_limits(&source, "cjson", "mol2", Policy::AllowLoss, &exact)
            .is_ok()
    );
    let short = ConversionLimits {
        max_output_bytes: output.len() - 1,
        ..exact
    };
    assert!(
        convert_with_report_with_limits(&source, "cjson", "mol2", Policy::AllowLoss, &short)
            .is_err()
    );
}

#[test]
fn isotope_loss_disables_correspondence_and_charge_verification() {
    let value = json!({"atoms":{"elements":{"number":[7]}, "formalCharges":[1], "coords":{"3d":[1,2,3]}},
        "partialCharges":{"USER_CHARGES":[-0.7]}});
    let report =
        convert_with_report(&value.to_string(), "cjson", "mol2", Policy::AllowLoss).unwrap();
    assert!(report.candidate_graph_verified);
    let after = parse_mol2_source_record(report.output.as_ref().unwrap()).unwrap();
    assert_eq!(after.molecule().atom(AtomIdx(0)).charge, 1);
    assert_eq!(after.atoms()[0].partial_charge, Some(-0.7));
    for charge in [1, -1] {
        let mut changed = value.clone();
        changed["atoms"]["formalCharges"] = json!([charge]);
        changed["atoms"]["isotopes"] = json!([15]);
        let report =
            convert_with_report(&changed.to_string(), "cjson", "mol2", Policy::AllowLoss).unwrap();
        assert!(!report.candidate_graph_verified);
        assert!(report.candidate_atom_mapping.is_empty());
        assert_eq!(status(&report, "isotopes"), Status::Dropped);
        assert_eq!(status(&report, "formal_charges"), Status::Unverified);
        assert_eq!(status(&report, "partial_charges"), Status::Unverified);
    }
}

#[test]
fn partial_and_coordinate_precision_cover_finite_exponents() {
    let mut bits = 0x243f_6a88_85a3_08d3_u64;
    for _ in 0..256 {
        bits ^= bits << 13;
        bits ^= bits >> 7;
        bits ^= bits << 17;
        for n in [f64::from_bits(bits), -0.0, 5e-324, f64::MAX] {
            if !n.is_finite() {
                continue;
            }
            let value = json!({"atoms":{"elements":{"number":[6]},"coords":{"3d":[n,-0.0,-n]}},"partialCharges":{"USER_CHARGES":[n]}});
            let report =
                convert_with_report(&value.to_string(), "cjson", "mol2", Policy::AllowLoss)
                    .unwrap();
            assert_eq!(status(&report, "partial_charges"), Status::Preserved, "{n}");
            assert_eq!(status(&report, "coordinates"), Status::Preserved, "{n}");
        }
    }
}

#[test]
fn output_limits_and_smiles_metadata_losses_are_explicit() {
    let limits = ConversionLimits {
        max_output_bytes: 10,
        ..Default::default()
    };
    assert!(
        convert_with_report_with_limits(SOURCE, "cjson", "mol2", Policy::AllowLoss, &limits)
            .is_err()
    );
    let report = convert_with_report(SOURCE, "cjson", "smiles", Policy::AllowLoss).unwrap();
    for field in ["partial_charges", "charge_method", "molecule_name"] {
        assert_eq!(status(&report, field), Status::Dropped);
    }
}
