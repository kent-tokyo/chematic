//! Original synthetic title-only fixtures, under the project license.
use chematic_mol::{
    ConversionLimits, ConversionPolicy as Policy, PreservationStatus as Status,
    convert_with_report, convert_with_report_with_limits, parse_cjson, parse_cjson_record,
    parse_mol2_source_record, write_cjson_record,
};
use serde_json::{Value, json};

fn source(root: Option<&str>, legacy: Option<&str>) -> String {
    let mut value = json!({"atoms":{"elements":{"number":[2]},"coords":{"3d":[0,0,0]}}});
    if let Some(name) = root {
        value["name"] = json!(name);
    }
    if let Some(name) = legacy {
        value["properties"] = json!({"name":name});
    }
    value.to_string()
}

#[test]
fn root_legacy_and_agreeing_names_remain_verbatim() {
    for name in ["helium standard", "ヘリウム β", "a\\b\"c"] {
        for text in [
            source(Some(name), None),
            source(None, Some(name)),
            source(Some(name), Some(name)),
        ] {
            let record = parse_cjson_record(&text).unwrap();
            assert_eq!(record.name(), Some(name));
            assert!(!record.has_conflicting_names());
            assert!(record.opaque_fields().is_empty());
            assert_eq!(write_cjson_record(&record), text);
            let report = convert_with_report(&text, "cjson", "mol2", Policy::AllowLoss).unwrap();
            assert!(report.candidate_graph_verified);
            let after = parse_mol2_source_record(report.output.as_ref().unwrap()).unwrap();
            assert_eq!(after.name(), Some(name));
            assert!(
                report
                    .observations
                    .iter()
                    .any(|o| o.field == "molecule_name" && o.status == Status::Preserved)
            );
            assert!(
                convert_with_report(&text, "cjson", "mol2", Policy::Strict)
                    .unwrap()
                    .output
                    .is_none()
            );
        }
    }
}

#[test]
fn conflict_retains_source_and_refuses_mol2_without_selecting_a_title() {
    let text = source(Some("root title"), Some("legacy title"));
    let record = parse_cjson_record(&text).unwrap();
    assert_eq!(record.name(), None);
    assert!(record.has_conflicting_names());
    assert_eq!(record.to_json()["name_conflict"], true);
    for policy in [Policy::Strict, Policy::AllowLoss, Policy::PreserveSource] {
        let same = convert_with_report(&text, "cjson", "cjson", policy).unwrap();
        assert_eq!(same.output.as_deref(), Some(text.as_str()));
    }
    for policy in [Policy::Strict, Policy::AllowLoss] {
        let report = convert_with_report(&text, "cjson", "mol2", policy).unwrap();
        assert_eq!(report.status, "refused");
        assert!(report.output.is_none());
        assert!(!report.candidate_graph_verified);
        assert!(report.candidate_atom_mapping.is_empty());
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == "molecule_name" && o.status == Status::Unverified)
        );
    }
    let smiles = convert_with_report(&text, "cjson", "smiles", Policy::AllowLoss).unwrap();
    assert!(
        smiles
            .observations
            .iter()
            .any(|o| o.field == "molecule_name" && o.status == Status::Dropped)
    );
}

#[test]
fn malformed_root_names_are_rejected_without_changing_legacy_graph_reader() {
    let mut value: Value = serde_json::from_str(&source(None, None)).unwrap();
    for bad in [Value::Null, json!(3), json!(false), json!([]), json!({})] {
        value["name"] = bad;
        assert!(parse_cjson_record(&value.to_string()).is_err());
        assert!(parse_cjson(&value.to_string()).is_ok());
    }
    assert!(
        parse_cjson_record(r#"{"atoms":{"elements":{"number":[2]}},"name":"a","n\u0061me":"b"}"#)
            .is_err()
    );
}

#[test]
fn unrepresentable_root_titles_keep_exact_source_and_never_inject_mol2() {
    for name in [
        "",
        " title",
        "title ",
        "head\n@<TRIPOS>ATOM",
        "a\0b",
        "@<tripos>ATOM",
    ] {
        let text = source(Some(name), None);
        let record = parse_cjson_record(&text).unwrap();
        assert_eq!(record.name(), Some(name));
        assert_eq!(write_cjson_record(&record), text);
        let report = convert_with_report(&text, "cjson", "mol2", Policy::AllowLoss).unwrap();
        assert_eq!(report.status, "refused");
        assert!(report.output.is_none());
    }
}

#[test]
fn mol2_to_cjson_emits_standard_name_and_matching_legacy_alias_within_output_limit() {
    let mol2 = "@<TRIPOS>MOLECULE\nヘリウム β\n1 0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 He 0 0 0 He\n@<TRIPOS>BOND\n";
    let report = convert_with_report(mol2, "mol2", "cjson", Policy::AllowLoss).unwrap();
    let output = report.output.unwrap();
    let value: Value = serde_json::from_str(&output).unwrap();
    assert_eq!(value["name"], "ヘリウム β");
    assert_eq!(value["properties"]["name"], value["name"]);
    let after = parse_cjson_record(&output).unwrap();
    assert!(!after.has_conflicting_names());
    assert_eq!(after.name(), Some("ヘリウム β"));
    let exact = ConversionLimits {
        max_output_bytes: output.len(),
        ..Default::default()
    };
    assert!(
        convert_with_report_with_limits(mol2, "mol2", "cjson", Policy::AllowLoss, &exact).is_ok()
    );
    let short = ConversionLimits {
        max_output_bytes: output.len() - 1,
        ..exact
    };
    assert!(
        convert_with_report_with_limits(mol2, "mol2", "cjson", Policy::AllowLoss, &short).is_err()
    );
}
