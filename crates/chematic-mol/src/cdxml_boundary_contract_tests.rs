use crate::cdxml::CdxmlParseLimits;
use crate::cdxml_document::*;
use serde_json::json;
const DOCUMENT: &str = "<CDXML>\n<page id=\"p\">\n<group id=\"g\">\n<text id=\"t\"/>\n</group>\n</page>\n<page id=\"q\">\n<caption id=\"c\"/>\n</page>\n</CDXML>\n";
#[test]
fn cdxml_document_structure_and_edit_errors_are_atomic() {
    let doc = CdxmlDocument::parse(DOCUMENT).unwrap();
    for input in [
        "<CDXML>\n<CDXML>\n</CDXML>",
        "</CDXML>",
        "<page/>",
        "<CDXML>\n<page>\n<group>\n</page>",
        "<CDXML>\n<page>\n</group>",
        "<CDXML>\n<page>",
    ] {
        assert!(
            !CdxmlDocument::parse(input)
                .unwrap_err()
                .to_string()
                .is_empty()
        );
    }
    for limits in [
        CdxmlParseLimits {
            max_input_bytes: 1,
            ..Default::default()
        },
        CdxmlParseLimits {
            max_lines: 1,
            ..Default::default()
        },
        CdxmlParseLimits {
            max_line_bytes: 1,
            ..Default::default()
        },
        CdxmlParseLimits {
            max_atoms: 0,
            max_bonds: 0,
            ..Default::default()
        },
    ] {
        assert!(
            CdxmlDocument::parse_with_limits(DOCUMENT, &limits)
                .unwrap_err()
                .to_string()
                .contains("limit")
        );
    }
    for edit in [
        CdxmlEdit::ReplaceObjectPath {
            page_id: "p".into(),
            path: vec![],
            raw_xml: "<text/>".into(),
        },
        CdxmlEdit::ReplaceObjectPath {
            page_id: "p".into(),
            path: vec![99],
            raw_xml: "<text/>".into(),
        },
        CdxmlEdit::RemoveObject {
            page_id: "q".into(),
            object_index: 99,
        },
        CdxmlEdit::ReplaceObject {
            page_id: "p".into(),
            object_index: 0,
            raw_xml: "".into(),
        },
        CdxmlEdit::SetObjectAttribute {
            page_id: "p".into(),
            object_index: 0,
            key: "".into(),
            value: "x".into(),
        },
    ] {
        assert!(!doc.apply(&edit).unwrap_err().to_string().is_empty());
        assert_eq!(doc.write(), DOCUMENT);
    }
    assert!(doc.apply_json_edit("{").is_err());
    let changed = doc
        .apply(&CdxmlEdit::ReplaceObjectPath {
            page_id: "q".into(),
            path: vec![0],
            raw_xml: "<caption id=\"replacement\"/>".into(),
        })
        .unwrap();
    assert!(changed.write().contains("replacement"));
    assert!(changed.write().contains("id=\"g\""));
    assert_eq!(doc.write(), DOCUMENT);
}

#[test]
fn cdxml_typed_presentation_requires_finite_correctly_shaped_values() {
    let doc = CdxmlDocument::parse(DOCUMENT).unwrap();
    let object = &doc.pages[0].children[0];
    assert!(object.transform().unwrap().is_none());
    assert!(object.z_order().unwrap().is_none());
    assert!(object.text_style().unwrap().is_none());
    for value in [
        json!(4),
        json!("bad 0 0 1 0 0"),
        json!("1 0 0 1 0"),
        json!("1 0 0 1 NaN 0"),
    ] {
        let mut object = object.clone();
        object.attributes.insert("Matrix".into(), value);
        assert!(object.transform().is_err());
    }
    let mut object = object.clone();
    object
        .attributes
        .insert("Matrix".into(), json!("1,0,0,1,2,3"));
    assert_eq!(
        object.transform().unwrap().unwrap().matrix,
        [1., 0., 0., 1., 2., 3.]
    );
    for value in [json!(4), json!("bad"), json!("NaN")] {
        object.attributes.insert("Size".into(), value);
        assert!(object.text_style().is_err());
    }
    for value in [json!(4), json!("bad")] {
        object.attributes.insert("ZOrder".into(), value);
        assert!(object.z_order().is_err());
    }
    let page = &doc.pages[0];
    assert!(page.bounding_box().unwrap().is_none());
    for value in [
        json!(4),
        json!("bad 0 1 2"),
        json!("0 1 2"),
        json!("0 1 2 NaN"),
    ] {
        let mut page = page.clone();
        page.attributes.insert("BoundingBox".into(), value);
        assert!(page.bounding_box().is_err());
    }
}
