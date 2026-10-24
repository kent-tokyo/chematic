use chematic_mol::{
    ConversionError, ConversionLimits, ConversionPolicy, convert_with_report, roundtrip_cml_source,
    roundtrip_cml_source_with_limits,
};

#[test]
fn standalone_source_retention_matches_report_without_claiming_a_graph() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/cml_source_preservation/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let source = case["source"].as_str().unwrap();
        assert_eq!(roundtrip_cml_source(source).unwrap(), source);
        let report =
            convert_with_report(source, "cml", "cml", ConversionPolicy::PreserveSource).unwrap();
        assert_eq!(report.output.unwrap(), source);
        assert!(!report.candidate_graph_verified);
        assert!(report.candidate_atom_mapping.is_empty());
    }
}

#[test]
fn exact_input_and_output_boundaries_are_independent_of_report_encoding() {
    let source = "<molecule xmlns=\"http://www.xml-cml.org/schema\"><name>&quot;α&quot;</name></molecule>\r\n";
    let mut limits = ConversionLimits::default();
    limits.mol2.max_input_bytes = source.len();
    limits.max_output_bytes = source.len();
    assert_eq!(
        roundtrip_cml_source_with_limits(source, &limits).unwrap(),
        source
    );
    limits.max_output_bytes -= 1;
    assert!(matches!(
        roundtrip_cml_source_with_limits(source, &limits),
        Err(ConversionError::ResourceLimit("output bytes"))
    ));
    limits.max_output_bytes += 1;
    limits.mol2.max_input_bytes -= 1;
    assert!(matches!(
        roundtrip_cml_source_with_limits(source, &limits),
        Err(ConversionError::ResourceLimit("input bytes"))
    ));
}

#[test]
fn element_depth_and_line_budgets_still_apply_without_projection() {
    let source = "<molecule xmlns=\"http://www.xml-cml.org/schema\"><atom/><atom/></molecule>";
    let mut limits = ConversionLimits::default();
    limits.mol2.max_atoms = 1;
    assert!(matches!(
        roundtrip_cml_source_with_limits(source, &limits),
        Err(ConversionError::ResourceLimit("CML atom/bond elements"))
    ));
    limits.mol2.max_atoms = 2;
    assert_eq!(
        roundtrip_cml_source_with_limits(source, &limits).unwrap(),
        source
    );
    limits.mol2.max_line_bytes = source.len() - 1;
    assert!(matches!(
        roundtrip_cml_source_with_limits(source, &limits),
        Err(ConversionError::ResourceLimit("line bytes"))
    ));
    let deep = format!(
        "<molecule xmlns=\"http://www.xml-cml.org/schema\">{}{}</molecule>",
        "<v>".repeat(64),
        "</v>".repeat(64)
    );
    assert!(matches!(
        roundtrip_cml_source(&deep),
        Err(ConversionError::ResourceLimit("CML XML depth"))
    ));
}
