use chematic_core::{AtomIdx, BondOrder};
use chematic_mol::{
    CjsonParseLimits, ConversionLimits, ConversionPolicy, PreservationStatus, convert_with_report,
    convert_with_report_json, parse_cjson, parse_cjson_record, parse_cjson_record_with_limits,
    write_cjson_record,
};
use serde_json::{Value, json};
const SOURCE: &str = include_str!("../../../validation/cjson_preservation/charged.cjson");

#[test]
fn exact_source_and_inspected_projection() {
    for source in [SOURCE.to_string(), SOURCE.replace('\n', "\r\n")] {
        let record = parse_cjson_record(&source).unwrap();
        assert_eq!(write_cjson_record(&record), source);
        assert_eq!(record.molecule().atom_count(), 4);
        assert_eq!(record.molecule().bond_count(), 2);
        assert_eq!(record.molecule().atom(AtomIdx(0)).isotope, Some(13));
        assert_eq!(record.molecule().atom(AtomIdx(1)).charge, 1);
        assert_eq!(record.molecule().atom(AtomIdx(2)).charge, -1);
        assert_eq!(record.molecule().atom(AtomIdx(3)).isotope, Some(2));
        assert_eq!(
            record.coordinates()[0].0.to_bits(),
            (-1.4175225513503917f64).to_bits()
        );
        assert_eq!(record.coordinates()[0].1.to_bits(), (-0.0f64).to_bits());
        assert_eq!(record.coordinates()[0].2.to_bits(), 1);
        assert_eq!(
            record.opaque_fields(),
            [
                "/atoms/coords/2d",
                "/atoms/elements/symbol",
                "/atoms/partialCharges",
                "/name",
                "/vendor~1a~0b"
            ]
        );
        assert!(
            record.formal_charges_explicit()
                && record.isotopes_explicit()
                && record.bond_orders_explicit()
        );
        let inspected = record.to_json();
        assert_eq!(inspected["schema"], "chematic.cjson-record.v1");
        assert_eq!(inspected["source"], source);
    }
}

#[test]
fn invalid_known_fields_are_never_silently_repaired() {
    let base: Value = serde_json::from_str(SOURCE).unwrap();
    for (pointer, values) in [
        (
            "/chemicalJson",
            vec![json!(0), json!(2), json!(1.0), Value::Null, json!("1")],
        ),
        (
            "/atoms/formalCharges",
            vec![
                json!([]),
                json!([0]),
                json!([0, 0, 0, 0, 0]),
                Value::Null,
                json!([0, 0, 0, 128]),
                json!([0, 0, 0, 1.5]),
            ],
        ),
        (
            "/atoms/isotopes",
            vec![
                json!([]),
                json!([0]),
                Value::Null,
                json!([0, 0, 0, -1]),
                json!([0, 0, 0, 65536]),
            ],
        ),
        ("/atoms/coords", vec![Value::Null, json!([])]),
        (
            "/atoms/coords/3d",
            vec![
                json!([]),
                json!([0, 0, 0]),
                Value::Null,
                json!([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "nan"]),
            ],
        ),
        (
            "/bonds/order",
            vec![
                json!([]),
                json!([1]),
                Value::Null,
                json!([2.2, 1]),
                json!([0, 1]),
                json!([5, 1]),
                json!([1.5000000001, 1]),
            ],
        ),
        (
            "/bonds/connections/index",
            vec![
                json!([0]),
                json!([0, 0]),
                json!([0, 1, 1, 0]),
                json!([0, 4]),
                json!([-1, 0]),
                json!([0, 1.0]),
            ],
        ),
        (
            "/atoms/elements/number",
            vec![json!([0]), json!([119]), json!([6.0]), Value::Null],
        ),
    ] {
        for value in values {
            let mut invalid = base.clone();
            *invalid.pointer_mut(pointer).unwrap() = value;
            assert!(
                parse_cjson_record(&invalid.to_string()).is_err(),
                "{pointer}: {invalid}"
            );
        }
    }
}

#[test]
fn duplicate_decoded_keys_are_rejected_recursively() {
    for source in [
        r#"{"atoms":{"elements":{"number":[6]}},"atoms":{"elements":{"number":[7]}}}"#,
        r#"{"atoms":{"elements":{"number":[6],"n\u0075mber":[7]}}}"#,
        r#"{"atoms":{"elements":{"number":[6]}},"vendor":{"x":1,"x":2}}"#,
        r#"{"atoms":{"elements":{"number":[6]}},"vendor":[{"x":1,"x":2}]}"#,
    ] {
        assert!(parse_cjson_record(source).is_err());
    }
}

#[test]
fn legacy_reader_remains_permissive() {
    let source = r#"{"atoms":{"elements":{"number":[6,8]},"formalCharges":[]},"bonds":{"connections":{"index":[0,1]},"order":[2.2]}}"#;
    let (legacy, _) = parse_cjson(source).unwrap();
    assert_eq!(legacy.bonds().next().unwrap().1.order, BondOrder::Double);
    assert!(parse_cjson_record(source).is_err());
}

#[test]
fn absent_optional_fields_and_local_aromatic_alias() {
    for order in ["1.5", "4"] {
        let source = format!(
            r#"{{"atoms":{{"elements":{{"number":[6,6]}}}},"bonds":{{"connections":{{"index":[0,1]}},"order":[{order}]}}}}"#
        );
        let record = parse_cjson_record(&source).unwrap();
        assert_eq!(record.declared_version(), None);
        assert!(!record.formal_charges_explicit());
        assert!(!record.isotopes_explicit());
        assert!(!record.coordinates_present());
        assert_eq!(
            record.molecule().bonds().next().unwrap().1.order,
            BondOrder::Aromatic
        );
        assert!(!record.molecule().atom(AtomIdx(0)).aromatic);
    }
}

#[test]
fn record_limits_cover_opaque_fields_and_graph() {
    for limits in [
        CjsonParseLimits {
            max_input_bytes: SOURCE.len() - 1,
            ..Default::default()
        },
        CjsonParseLimits {
            max_atoms: 3,
            ..Default::default()
        },
        CjsonParseLimits {
            max_bonds: 1,
            ..Default::default()
        },
        CjsonParseLimits {
            max_array_items: 11,
            ..Default::default()
        },
        CjsonParseLimits {
            max_string_bytes: 3,
            ..Default::default()
        },
        CjsonParseLimits {
            max_json_depth: 2,
            ..Default::default()
        },
    ] {
        assert!(parse_cjson_record_with_limits(SOURCE, &limits).is_err());
    }
}

#[test]
fn checked_self_copy_is_exact_for_every_policy() {
    for policy in [
        ConversionPolicy::Strict,
        ConversionPolicy::AllowLoss,
        ConversionPolicy::PreserveSource,
    ] {
        let report = convert_with_report(SOURCE, ".CJSON", "cjson", policy).unwrap();
        assert_eq!(report.status, "completed");
        assert_eq!(report.output.as_deref(), Some(SOURCE));
        assert_eq!(report.verification_scope, "exact_source_bytes");
        assert!(!report.candidate_graph_verified);
        assert!(report.candidate_atom_mapping.is_empty());
        assert!(report.observations.iter().all(|o| matches!(
            o.status,
            PreservationStatus::Preserved | PreservationStatus::OpaqueRetained
        )));
    }
}

#[test]
fn smiles_charge_isotope_mapping_and_loss_accounting() {
    let report = convert_with_report(SOURCE, "cjson", "smi", ConversionPolicy::AllowLoss).unwrap();
    assert_eq!(report.status, "completed_with_loss");
    assert!(report.candidate_graph_verified);
    let candidate = chematic_smiles::parse(report.output.as_ref().unwrap()).unwrap();
    let source = parse_cjson_record(SOURCE).unwrap();
    for map in &report.candidate_atom_mapping {
        assert_eq!(map.source_id as usize, map.source_index);
        let a = source.molecule().atom(AtomIdx(map.source_index as u32));
        let b = candidate.atom(AtomIdx(map.candidate_index as u32));
        assert_eq!(
            (a.element, a.charge, a.isotope),
            (b.element, b.charge, b.isotope)
        );
    }
    assert_eq!(report.candidate_atom_mapping.len(), 4);
    assert_eq!(report.candidate_bond_mapping.len(), 2);
    for field in ["source_text", "coordinates", "opaque_fields"] {
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == field && o.status == PreservationStatus::Dropped)
        );
    }
    let strict = convert_with_report(SOURCE, "cjson", "smiles", ConversionPolicy::Strict).unwrap();
    assert_eq!(strict.status, "refused");
    assert!(strict.output.is_none());
    assert!(strict.candidate_graph_verified);
}

#[test]
fn checked_errors_and_output_envelope_limits() {
    for (target, policy) in [
        ("cml", "allow_loss"),
        ("smiles", "preserve_source"),
        ("cjson", "typo"),
    ] {
        assert!(
            convert_with_report_json(
                SOURCE,
                "cjson",
                target,
                policy,
                &ConversionLimits::default()
            )
            .is_err()
        );
    }
    let limits = ConversionLimits {
        max_output_bytes: SOURCE.len(),
        ..Default::default()
    };
    assert!(convert_with_report_json(SOURCE, "cjson", "cjson", "strict", &limits).is_err());
}

#[test]
fn branch_mapping_follows_emitted_order_and_reports_missing_defaults() {
    let source = r#"{"atoms":{"elements":{"number":[6,8,7,6]},"isotopes":[13,18,15,12]},"bonds":{"connections":{"index":[0,1,0,2,1,3]}}}"#;
    let report =
        convert_with_report(source, "cjson", "smiles", ConversionPolicy::AllowLoss).unwrap();
    assert!(report.candidate_graph_verified);
    assert_eq!(
        report
            .candidate_atom_mapping
            .iter()
            .map(|m| m.source_index)
            .collect::<Vec<_>>(),
        [0, 1, 3, 2]
    );
    let candidate = chematic_smiles::parse(report.output.as_ref().unwrap()).unwrap();
    for map in &report.candidate_atom_mapping {
        assert_eq!(
            candidate.atom(AtomIdx(map.candidate_index as u32)).isotope,
            Some([13, 18, 15, 12][map.source_index])
        );
    }
    for field in ["formal_charge_defaults", "bond_order_defaults"] {
        assert!(
            report
                .observations
                .iter()
                .any(|o| o.field == field && o.status == PreservationStatus::Inferred)
        );
    }
    assert!(
        !report
            .observations
            .iter()
            .any(|o| o.field == "coordinates" || o.field == "opaque_fields")
    );
}
