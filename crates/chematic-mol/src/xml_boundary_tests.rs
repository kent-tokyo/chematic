use crate::{cdxml::*, cml::*, mrv::*};

fn cml(body: &str) -> String {
    format!("<cml>\n<molecule>\n{body}\n</molecule>\n</cml>\n")
}
fn cdxml(body: &str) -> String {
    format!(
        "<CDXML>\n<page id=\"p1\">\n<fragment id=\"f1\">\n{body}\n</fragment>\n</page>\n</CDXML>\n"
    )
}

#[test]
fn strict_cml_rejects_invalid_references_coordinates_and_xml() {
    for (input, diagnostic) in [
        ("<cml/>".into(), "missing molecule"),
        (cml(""), "no atom"),
        ("<cml><molecule></cml>".into(), "malformed"),
        (cml("<atom id=\"a1\" elementType=\"Qq\"/>"), "element"),
        (
            cml("<atom id=\"a1\" elementType=\"C\" x2=\"NaN\" y2=\"0\"/>"),
            "finite",
        ),
        (
            cml("<atom id=\"a1\" elementType=\"C\"/>\n<bond atomRefs2=\"a1 absent\" order=\"1\"/>"),
            "ref",
        ),
        (
            cml("<atom id=\"a1\" elementType=\"C\"/>\n<bond atomRefs2=\"a1\" order=\"1\"/>"),
            "atomRefs2",
        ),
        (
            cml(
                "<atom id=\"a1\" elementType=\"C\"/>\n<atom id=\"a2\" elementType=\"C\"/>\n<bond atomRefs2=\"a1 a2\" order=\"bad\"/>",
            ),
            "order",
        ),
    ] {
        let error = parse_cml_strict(&input).err().unwrap();
        assert!(error.to_string().contains(diagnostic), "{input}: {error}");
    }
}

#[test]
fn cml_limits_apply_before_graph_construction() {
    let input = cml(
        "<atom id=\"a1\" elementType=\"C\"/>\n<atom id=\"a2\" elementType=\"C\"/>\n<bond atomRefs2=\"a1 a2\" order=\"1\"/>",
    );
    for limits in [
        CmlParseLimits {
            max_input_bytes: 0,
            ..Default::default()
        },
        CmlParseLimits {
            max_line_bytes: 1,
            ..Default::default()
        },
        CmlParseLimits {
            max_lines: 1,
            ..Default::default()
        },
        CmlParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        CmlParseLimits {
            max_bonds: 0,
            ..Default::default()
        },
        CmlParseLimits {
            max_xml_elements: 1,
            ..Default::default()
        },
    ] {
        let error = parse_cml_with_limits(&input, &limits).err().unwrap();
        assert!(matches!(error, CmlError::ResourceLimit { .. }));
        assert!(error.to_string().contains("limit"));
    }
}

#[test]
fn mrv_rejects_queries_and_ambiguous_atom_or_bond_records() {
    for (body, diagnostic) in [
        (
            "<atomArray><atom id=\"a1\" elementType=\"Qq\"/></atomArray>",
            "element",
        ),
        (
            "<atomArray><atom id=\"a1\"/><atom id=\"a1\"/></atomArray>",
            "duplicate",
        ),
        (
            "<atomArray><atom id=\"a1\" elementType=\"R\"/></atomArray>",
            "R-group",
        ),
        (
            "<atomArray><atom id=\"a1\" elementType=\"*\"/></atomArray>",
            "wildcard",
        ),
        (
            "<atomArray><atom id=\"a1\" mrvStereoGroup=\"1\"/></atomArray>",
            "stereo group",
        ),
        (
            "<atomArray><atom id=\"a1\" radical=\"bad\"/></atomArray>",
            "radical",
        ),
        (
            "<atomArray><atom id=\"a1\"/></atomArray><bondArray><bond atomRefs2=\"a1 missing\"/></bondArray>",
            "ref",
        ),
        (
            "<atomArray><atom id=\"a1\"/></atomArray><bondArray><bond atomRefs2=\"a1\"/></bondArray>",
            "atomRefs2",
        ),
        (
            "<atomArray><atom id=\"a1\"/><atom id=\"a2\"/></atomArray><bondArray><bond atomRefs2=\"a1 a2\" order=\"bad\"/></bondArray>",
            "order",
        ),
        (
            "<atomArray><atom id=\"a1\"/><atom id=\"a2\"/></atomArray><bondArray><bond atomRefs2=\"a1 a2\" queryType=\"any\"/></bondArray>",
            "query",
        ),
    ] {
        let input = cml(body);
        let error = parse_mrv(&input).err().unwrap();
        assert!(error.to_string().contains(diagnostic), "{input}: {error}");
    }
    let input = cml(
        "<atomArray><atom id=\"a1\"/><atom id=\"a2\"/></atomArray><bondArray><bond atomRefs2=\"a1 a2\"><bondStereo>BAD</bondStereo></bond></bondArray>",
    );
    assert!(matches!(
        parse_mrv(&input),
        Err(MrvError::InvalidBondStereo { .. })
    ));
    for limits in [
        MrvParseLimits {
            max_input_bytes: 1,
            ..Default::default()
        },
        MrvParseLimits {
            max_depth: 1,
            ..Default::default()
        },
        MrvParseLimits {
            max_attr_len: 0,
            ..Default::default()
        },
    ] {
        assert!(
            parse_mrv_with_limits(&input, &limits)
                .err()
                .unwrap()
                .to_string()
                .contains("limit")
        );
    }
    assert!(
        parse_mrv("<!DOCTYPE cml><cml/>")
            .err()
            .unwrap()
            .to_string()
            .contains("DOCTYPE")
    );
}

#[test]
fn cdxml_rejects_invalid_atom_coordinates_and_bond_endpoints() {
    for (body, diagnostic) in [
        ("<n id=\"1\" Element=\"999\" p=\"0 0\"/>", "atomic number"),
        ("<n id=\"1\" p=\"0\"/>", "coords"),
        ("<n id=\"1\" p=\"0 0\"/>\n<b id=\"2\" B=\"1\"/>", "missing"),
        (
            "<n id=\"1\" p=\"0 0\"/>\n<b id=\"2\" B=\"1\" E=\"9\"/>",
            "ref",
        ),
    ] {
        let input = cdxml(body);
        let error = parse_cdxml(&input).err().unwrap();
        assert!(error.to_string().contains(diagnostic), "{input}: {error}");
    }
    let input =
        cdxml("<n id=\"1\" p=\"0 0\"/>\n<n id=\"2\" p=\"1 0\"/>\n<b id=\"3\" B=\"1\" E=\"2\"/>");
    for limits in [
        CdxmlParseLimits {
            max_input_bytes: 0,
            ..Default::default()
        },
        CdxmlParseLimits {
            max_line_bytes: 1,
            ..Default::default()
        },
        CdxmlParseLimits {
            max_lines: 1,
            ..Default::default()
        },
        CdxmlParseLimits {
            max_attribute_bytes: 0,
            ..Default::default()
        },
        CdxmlParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        CdxmlParseLimits {
            max_bonds: 0,
            ..Default::default()
        },
        CdxmlParseLimits {
            max_fragments: 0,
            ..Default::default()
        },
    ] {
        let error = parse_cdxml_with_limits(&input, &limits).err().unwrap();
        assert!(error.to_string().contains("limit"), "{error}");
    }
}
