use super::*;

#[test]
fn candidate_verification_rejects_graph_coordinates_and_name_corruption() {
    let source = r#"{"name":"alpha","atoms":{"elements":{"number":[6]},"coords":{"3d":[1,2,3]}}}"#;
    let without_coords = r#"{"name":"alpha","atoms":{"elements":{"number":[6]}}}"#;
    let record = crate::parse_cjson_record(source).unwrap();
    let output = project_cjson_to_cml(source).unwrap();
    let absent = project_cjson_to_cml(without_coords).unwrap();
    let limits = ConversionLimits::default();
    assert!(verify(&record, &output, &limits).is_ok());

    for (candidate, reason) in [
        (
            output.replace("elementType=\"C\"", "elementType=\"N\""),
            "indexed graph",
        ),
        (absent, "coordinate presence"),
        (
            output.replace(
                &format!("x3=\"{:.17e}\"", 1.0),
                &format!("x3=\"{:.17e}\"", 4.0),
            ),
            "coordinate bits",
        ),
        (
            output.replace("<name>alpha</name>", "<name>beta</name>"),
            "verify name",
        ),
        (
            "<molecule><atomArray><atom elementType=\"Xx\"/></atomArray></molecule>".into(),
            "reparse candidate",
        ),
        (output.replace("</molecule>", ""), "reparse candidate"),
    ] {
        assert_ne!(candidate, output, "{reason}");
        let error = verify(&record, &candidate, &limits).unwrap_err();
        assert!(error.to_string().contains(reason), "{reason}: {error}");
    }
}
