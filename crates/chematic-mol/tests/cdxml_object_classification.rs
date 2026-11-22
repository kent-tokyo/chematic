use chematic_mol::{CdxmlDocument, CdxmlObjectKind};

#[test]
fn document_exposes_presentation_object_classes_and_retains_unknown_xml() {
    let source = "<CDXML><page id=\"p\"><group/><text/><caption/><graphic/><curve/><table/><scheme/><bracket_attachment/><vendor-object id=\"x\"/></page></CDXML>";
    let document = CdxmlDocument::parse(source).unwrap();
    let objects = document.pages[0].child_objects(&[]);
    let expected = [
        CdxmlObjectKind::Group,
        CdxmlObjectKind::Text,
        CdxmlObjectKind::Caption,
        CdxmlObjectKind::Graphic,
        CdxmlObjectKind::Curve,
        CdxmlObjectKind::Table,
        CdxmlObjectKind::Scheme,
        CdxmlObjectKind::BracketAttachment,
        CdxmlObjectKind::Unsupported("vendor-object".into()),
    ];
    assert_eq!(
        objects.iter().map(|o| o.kind()).collect::<Vec<_>>(),
        expected
    );
    assert!(objects.last().unwrap().raw_xml.contains("vendor-object"));
}
