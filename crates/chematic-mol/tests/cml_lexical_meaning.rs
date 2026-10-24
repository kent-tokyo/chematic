use chematic_core::{Molecule, implicit_hcount};
use chematic_mol::{
    CmlParseLimits, ConversionPolicy, convert_with_report, parse_cml, parse_cml_strict,
    parse_cml_strict_with_limits, parse_cml_with_3d, parse_cml_with_limits,
};
use serde_json::{Value, json};

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/cml_lexical_meaning/cases.json"
    ))
    .unwrap()
}

fn fields(mol: &Molecule) -> Value {
    let mut bonds: Vec<_> = mol
        .bonds()
        .map(|(_, b)| {
            [
                b.atom1.0.min(b.atom2.0),
                b.atom1.0.max(b.atom2.0),
                b.order.order_value().unwrap() as u32,
            ]
        })
        .collect();
    bonds.sort();
    json!({
        "numbers":mol.atoms().map(|(_, a)| a.element.atomic_number()).collect::<Vec<_>>(),
        "charges":mol.atoms().map(|(_, a)| a.charge).collect::<Vec<_>>(),
        "isotopes":mol.atoms().map(|(_, a)| a.isotope.unwrap_or(0)).collect::<Vec<_>>(),
        "hydrogens":mol.atoms().map(|(i, _)| implicit_hcount(mol, i)).collect::<Vec<_>>(),
        "atom_maps":mol.atoms().map(|(_, a)| a.atom_map.unwrap_or(0)).collect::<Vec<_>>(),
        "bonds":bonds,
    })
}

#[test]
fn formatting_and_opaque_xml_cannot_change_the_molecular_projection() {
    for c in packet()["cases"].as_array().unwrap() {
        let text = c["input"].as_str().unwrap();
        let mut expected = c["expected"].clone();
        let coords = expected
            .as_object_mut()
            .unwrap()
            .remove("coords2d")
            .unwrap();
        for reader in [parse_cml, parse_cml_strict] {
            let (mol, actual_coords) = reader(text).unwrap();
            assert_eq!(fields(&mol), expected, "{}", c["id"]);
            assert_eq!(json!(actual_coords), coords, "{}", c["id"]);
        }
        let (mol, points) = parse_cml_with_3d(text).unwrap();
        assert_eq!(fields(&mol), expected, "{}", c["id"]);
        assert!(points.is_none());
        let report =
            convert_with_report(text, "cml", "cml", ConversionPolicy::PreserveSource).unwrap();
        assert_eq!(report.output.as_deref(), Some(text));
        assert!(!report.candidate_graph_verified);
    }
}

#[test]
fn raw_input_limits_apply_before_opaque_content_is_discarded() {
    let source = packet()["cases"][11]["input"].as_str().unwrap().to_string();
    for reader in [parse_cml_with_limits, parse_cml_strict_with_limits] {
        let exact = CmlParseLimits {
            max_input_bytes: source.len(),
            ..Default::default()
        };
        assert!(reader(&source, &exact).is_ok());
        assert!(
            reader(
                &source,
                &CmlParseLimits {
                    max_input_bytes: source.len() - 1,
                    ..exact
                }
            )
            .is_err()
        );
        assert!(
            reader(
                &source,
                &CmlParseLimits {
                    max_lines: 2,
                    ..Default::default()
                }
            )
            .is_err()
        );
        let long = source.replace("opaque example", &"x".repeat(4096));
        assert!(
            reader(
                &long,
                &CmlParseLimits {
                    max_line_bytes: 1024,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
}

#[test]
fn lenient_boundary_and_strict_rejections_remain_distinct() {
    assert_eq!(parse_cml("<cml/>").unwrap().0.atom_count(), 0);
    for text in [
        "<molecule><!--\n<atom elementType='C'/>\n--></molecule>",
        "<molecule><![CDATA[\n<atom elementType='C'/>\n]]></molecule>",
    ] {
        assert!(parse_cml_strict(text).is_err());
        assert_eq!(parse_cml(text).unwrap().0.atom_count(), 0);
    }
    for suffix in ["<!-- unfinished", "<![CDATA[unfinished", "<?unfinished"] {
        let text = format!("<molecule><atom id='a1' elementType='O'/>{suffix}");
        assert!(parse_cml_strict(&text).is_err());
    }
}
