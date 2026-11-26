use chematic_mol::{ExtxyzWriter, XyzProperty, XyzPropertyKind, XyzValue, parse_extxyz};

#[test]
fn invalid_header_metadata_does_not_write_a_partial_frame() {
    let base = parse_extxyz("1\n\nC 0 0 0\n").unwrap();
    let mut cases = Vec::new();
    let mut bad_key = base.clone();
    bad_key.info.push(("bad key".into(), "value".into()));
    cases.push(bad_key);
    let mut bad_value = base.clone();
    bad_value.info.push(("note".into(), "control\u{1}".into()));
    cases.push(bad_value);
    let mut bad_property = base;
    bad_property.properties.push(XyzProperty {
        name: "bad name".into(),
        kind: XyzPropertyKind::Real,
        count: 1,
        values: vec![vec![XyzValue::Real(1.0)]],
    });
    cases.push(bad_property);

    for frame in cases {
        let mut output = Vec::new();
        let error = ExtxyzWriter::new(&mut output)
            .write_frame(&frame)
            .unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(output.is_empty());
    }
}
