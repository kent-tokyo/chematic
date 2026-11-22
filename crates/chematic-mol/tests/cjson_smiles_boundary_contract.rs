use chematic_mol::{
    ConversionError, ConversionLimits, ConversionPolicy, PreservationStatus, convert_with_report,
    convert_with_report_with_limits,
};

#[test]
fn aromatic_cjson_projection_discloses_unverified_chemical_aromaticity() {
    let source = r#"{"atoms":{"elements":{"number":[6,6]}},"bonds":{"connections":{"index":[0,1]},"order":[1.5]}}"#;
    let report =
        convert_with_report(source, "cjson", "smiles", ConversionPolicy::AllowLoss).unwrap();
    assert!(
        report
            .observations
            .iter()
            .any(|o| o.field == "aromaticity_interpretation"
                && o.status == PreservationStatus::Unverified
                && o.stage == "projection")
    );
}

#[test]
fn cjson_smiles_projection_respects_atom_and_output_limits() {
    let source = serde_json::json!({"atoms":{"elements":{"number":vec![6; 513]}}}).to_string();
    assert!(matches!(
        convert_with_report(&source, "cjson", "smiles", ConversionPolicy::AllowLoss),
        Err(ConversionError::ResourceLimit(
            "SMILES projection atoms (512)"
        ))
    ));
    let source =
        r#"{"atoms":{"elements":{"number":[6,6,8]}},"bonds":{"connections":{"index":[0,1,1,2]}}}"#;
    let limits = ConversionLimits {
        max_output_bytes: 1,
        ..Default::default()
    };
    assert!(matches!(
        convert_with_report_with_limits(
            source,
            "cjson",
            "smiles",
            ConversionPolicy::AllowLoss,
            &limits
        ),
        Err(ConversionError::ResourceLimit("output bytes"))
    ));
    assert!(matches!(
        convert_with_report_with_limits(
            source,
            "cjson",
            "cjson",
            ConversionPolicy::Strict,
            &limits
        ),
        Err(ConversionError::ResourceLimit("output bytes"))
    ));
}
