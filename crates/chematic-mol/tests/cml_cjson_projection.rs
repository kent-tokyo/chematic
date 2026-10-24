use chematic_mol::{ConversionLimits, project_cml_to_cjson, project_cml_to_cjson_with_limits};
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/cml_cjson_projection/cases.json"
    ))
    .unwrap()
}
fn variants(text: &str) -> Vec<String> {
    vec![
        text.into(),
        text.replace("><", ">\n<"),
        text.replace("><", ">\r\n<"),
    ]
}
#[test]
fn profile_preserves_and_refuses_across_xml_layouts() {
    let p = packet();
    for c in p["valid"].as_array().unwrap() {
        for text in variants(c["input"].as_str().unwrap()) {
            let out = project_cml_to_cjson(&text).unwrap_or_else(|e| panic!("{}: {e}", c["id"]));
            let r = chematic_mol::parse_cjson_record(&out).unwrap();
            let (source, xyz) = chematic_mol::parse_cml_with_3d(&text).unwrap();
            assert_eq!(source.atom_count(), r.molecule().atom_count());
            assert_eq!(source.bond_count(), r.molecule().bond_count());
            assert_eq!(xyz.is_some(), r.coordinates_present());
            if let Some(xyz) = xyz {
                for (a, &(x, y, z)) in xyz.iter().zip(r.coordinates()) {
                    assert_eq!(a.map(f64::to_bits), [x, y, z].map(f64::to_bits));
                }
            }
        }
    }
    for c in p["invalid"].as_array().unwrap() {
        for text in variants(c["input"].as_str().unwrap()) {
            assert!(project_cml_to_cjson(&text).is_err(), "{}", c["id"]);
        }
    }
}
#[test]
fn limits_and_numeric_extremes_are_checked() {
    let p = packet();
    let text = p["valid"][0]["input"].as_str().unwrap();
    for kind in 0..6 {
        let mut l = ConversionLimits::default();
        match kind {
            0 => l.mol2.max_input_bytes = text.len() - 1,
            1 => l.mol2.max_atoms = 1,
            2 => l.mol2.max_bonds = 0,
            3 => l.max_output_bytes = 1,
            4 => l.mol2.max_lines = 0,
            _ => l.mol2.max_line_bytes = 1,
        }
        assert!(project_cml_to_cjson_with_limits(text, &l).is_err());
    }
    let text = "<molecule><atom id='a' elementType='O' x3='-0' y3='5e-324' z3='1.7976931348623157e308'/></molecule>";
    let out = project_cml_to_cjson(text).unwrap();
    let record = chematic_mol::parse_cjson_record(&out).unwrap();
    assert_eq!(record.coordinates()[0].0.to_bits(), (-0.0f64).to_bits());
    assert_eq!(record.coordinates()[0].1.to_bits(), 1);
    assert_eq!(record.coordinates()[0].2, f64::MAX);
}
#[test]
fn name_segments_and_legacy_xml_are_supported_without_inventing_a_name() {
    let out = project_cml_to_cjson("<molecule title='ABC'><name>A<!--ignored--><![CDATA[B]]><?pi ignored?>C</name><atom id='a' elementType='O'/></molecule>").unwrap();
    assert_eq!(serde_json::from_str::<Value>(&out).unwrap()["name"], "ABC");
    let out = project_cml_to_cjson("<molecule><atom id='a' elementType='O'/></molecule>").unwrap();
    assert!(
        serde_json::from_str::<Value>(&out)
            .unwrap()
            .get("name")
            .is_none()
    );
}

#[test]
fn explicit_h_nodes_and_isotope_zero_aliases_do_not_silently_change() {
    let text = "<molecule><atom id='c' elementType='C' hydrogenCount='1'/><atom id='d' elementType='H' isotopeNumber='2'/><bond atomRefs2='c d' order='1'/></molecule>";
    // CML total=1 includes the explicit D node, leaving zero non-node H;
    // CJSON would infer three, so returning ordinary CJSON would change it.
    let error = project_cml_to_cjson(text).err().unwrap().to_string();
    assert!(error.contains("hydrogen count"));
    assert!(
        project_cml_to_cjson("<molecule>&#160;<atom id='a' elementType='O'/></molecule>").is_err()
    );
    for mass in ["0", "00", "+0"] {
        let text =
            format!("<molecule><atom id='a' elementType='O' isotopeNumber='{mass}'/></molecule>");
        assert!(project_cml_to_cjson(&text).is_err());
    }
}

#[test]
fn unsupported_bond_orders_report_the_parsed_order() {
    for (token, expected) in [
        ("0", "Zero"),
        ("any", "QueryAny"),
        ("S/D", "QuerySingleOrDouble"),
        ("S/A", "QuerySingleOrAromatic"),
        ("D/A", "QueryDoubleOrAromatic"),
        ("A", "Aromatic"),
    ] {
        let text = format!(
            "<molecule><atom id='a' elementType='C'/><atom id='b' elementType='N'/><bond atomRefs2='a b' order='{token}'/></molecule>"
        );
        let error = project_cml_to_cjson(&text).unwrap_err().to_string();
        assert!(
            error.contains(&format!("cannot preserve bond order {expected}")),
            "{token}: {error}"
        );
    }
}
