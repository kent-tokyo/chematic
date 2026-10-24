use chematic_mol::{
    ConversionError, ConversionLimits, project_cjson_to_cml, project_cjson_to_cml_with_limits,
};

fn packet() -> serde_json::Value {
    serde_json::from_str(include_str!(
        "../../../validation/cjson_cml_projection/cases.json"
    ))
    .unwrap()
}

#[test]
fn source_graph_xyz_and_name_are_retained_in_cml() {
    for case in packet()["valid"].as_array().unwrap() {
        let output = project_cjson_to_cml(case["input"].as_str().unwrap()).unwrap();
        let (mol, xyz) = chematic_mol::parse_cml_with_3d(&output).unwrap();
        let expected = &case["expected"];
        assert_eq!(
            mol.atom_count(),
            expected["numbers"].as_array().unwrap().len()
        );
        for (i, atom) in mol.atoms() {
            let n = i.0 as usize;
            assert_eq!(
                u64::from(atom.element.atomic_number()),
                expected["numbers"][n].as_u64().unwrap()
            );
            assert_eq!(
                i64::from(atom.charge),
                expected["charges"][n].as_i64().unwrap()
            );
            assert_eq!(
                u64::from(atom.isotope.unwrap_or(0)),
                expected["isotopes"][n].as_u64().unwrap()
            );
        }
        assert_eq!(
            mol.bond_count(),
            expected["bonds"].as_array().unwrap().len()
        );
        if expected["coordinates"].is_null() {
            assert!(xyz.is_none());
        } else {
            let expected = expected["coordinates"]
                .as_array()
                .unwrap()
                .iter()
                .flat_map(|r| r.as_array().unwrap())
                .map(|v| v.as_f64().unwrap().to_bits())
                .collect::<Vec<_>>();
            assert_eq!(
                xyz.unwrap()
                    .into_iter()
                    .flatten()
                    .map(f64::to_bits)
                    .collect::<Vec<_>>(),
                expected
            );
        }
        let xml = roxmltree_for_test(&output);
        assert_eq!(xml, expected["name"].as_str().map(str::to_owned));
    }
}

// Read the name without reusing the writer's escaping implementation.
fn roxmltree_for_test(output: &str) -> Option<String> {
    roxmltree::Document::parse(output)
        .unwrap()
        .root_element()
        .children()
        .find(|n| n.has_tag_name(("http://www.xml-cml.org/schema", "name")))
        .map(|n| n.text().unwrap_or("").to_owned())
}

#[test]
fn refuses_unrepresented_source_information_and_malformed_graphs() {
    for case in packet()["invalid"].as_array().unwrap() {
        assert!(
            project_cjson_to_cml(case["input"].as_str().unwrap()).is_err(),
            "{}",
            case["id"]
        );
    }
}

#[test]
fn xml_name_cannot_introduce_graph_nodes_or_invalid_xml_characters() {
    let name = "</name><atom id='forged' elementType='H'/><name>&amp;\r\n\t";
    let source = serde_json::json!({"name":name,"atoms":{"elements":{"number":[6]}}});
    let output = project_cjson_to_cml(&source.to_string()).unwrap();
    assert_eq!(
        chematic_mol::parse_cml_with_3d(&output)
            .unwrap()
            .0
            .atom_count(),
        1
    );
    assert_eq!(roxmltree_for_test(&output).as_deref(), Some(name));
    for character in ['\0', '\u{1f}', '\u{fffe}', '\u{ffff}'] {
        let mut invalid = source.clone();
        invalid["name"] = serde_json::json!(format!("name{character}"));
        assert!(project_cjson_to_cml(&invalid.to_string()).is_err());
    }
}

#[test]
fn preserves_binary64_extremes_and_enforces_output_limit() {
    let points = [-0.0, f64::from_bits(1), f64::MAX];
    let source =
        serde_json::json!({"atoms":{"elements":{"number":[6]},"coords":{"3d":points}}}).to_string();
    let output = project_cjson_to_cml(&source).unwrap();
    let (_, coords) = chematic_mol::parse_cml_with_3d(&output).unwrap();
    assert_eq!(
        coords.unwrap()[0].map(f64::to_bits),
        points.map(f64::to_bits)
    );
    let limits = ConversionLimits {
        max_output_bytes: output.len() - 1,
        ..Default::default()
    };
    assert!(matches!(
        project_cjson_to_cml_with_limits(&source, &limits),
        Err(ConversionError::ResourceLimit("output bytes"))
    ));
    let limits = ConversionLimits {
        max_output_bytes: output.len(),
        ..Default::default()
    };
    assert_eq!(
        project_cjson_to_cml_with_limits(&source, &limits).unwrap(),
        output
    );
}
