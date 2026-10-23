use chematic_mol::{CmlError, CmlParseLimits, parse_cml_with_3d, parse_cml_with_3d_with_limits};
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/cml_namespace_projection/cases.json"
    ))
    .unwrap()
}
fn layouts(text: &str) -> [String; 3] {
    let lf = text.replace("><", ">\n<");
    [text.into(), lf.clone(), lf.replace('\n', "\r\n")]
}

#[test]
fn namespace_spelling_and_foreign_annotations_cannot_change_graph() {
    for c in packet()["valid"].as_array().unwrap() {
        for text in layouts(c["input"].as_str().unwrap()) {
            let (mol, coords) = parse_cml_with_3d(&text).unwrap();
            let expected = &c["expected"];
            let xyz: Vec<Vec<f64>> =
                serde_json::from_value(expected["coordinates"].clone()).unwrap();
            let coords = coords.unwrap();
            assert_eq!(coords.len(), xyz.len());
            for (actual, wanted) in coords.iter().flatten().zip(xyz.iter().flatten()) {
                assert_eq!(actual.to_bits(), wanted.to_bits(), "{}", c["id"]);
            }
            assert_eq!(
                serde_json::to_value(
                    mol.atoms()
                        .map(|(_, a)| a.element.atomic_number())
                        .collect::<Vec<_>>()
                )
                .unwrap(),
                expected["numbers"]
            );
            assert_eq!(
                serde_json::to_value(mol.atoms().map(|(_, a)| a.charge).collect::<Vec<_>>())
                    .unwrap(),
                expected["formal_charges"]
            );
            assert_eq!(
                serde_json::to_value(
                    mol.atoms()
                        .map(|(_, a)| a.isotope.unwrap_or(0))
                        .collect::<Vec<_>>()
                )
                .unwrap(),
                expected["isotopes"]
            );
            assert_eq!(
                mol.bond_count(),
                expected["bonds"].as_array().unwrap().len()
            );
            for ((_, bond), wanted) in mol.bonds().zip(expected["bonds"].as_array().unwrap()) {
                assert_eq!(u64::from(bond.atom1.0), wanted[0].as_u64().unwrap());
                assert_eq!(u64::from(bond.atom2.0), wanted[1].as_u64().unwrap());
                assert_eq!(
                    bond.order.order_value().unwrap() as f64,
                    wanted[2].as_f64().unwrap()
                );
            }
        }
    }
}

#[test]
fn invalid_xml_scope_ids_and_bonds_produce_no_graph() {
    for c in packet()["invalid"].as_array().unwrap() {
        for text in layouts(c["input"].as_str().unwrap()) {
            assert!(parse_cml_with_3d(&text).is_err(), "{}", c["id"]);
        }
    }
}

#[test]
fn foreign_subtrees_are_still_counted_in_xml_resource_budgets() {
    let source = packet()["valid"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "foreign-ancestor")
        .unwrap()["input"]
        .as_str()
        .unwrap()
        .to_owned();
    let exact = CmlParseLimits {
        max_atoms: 2,
        max_bonds: 1,
        ..Default::default()
    };
    assert_eq!(
        parse_cml_with_3d_with_limits(&source, &exact)
            .unwrap()
            .0
            .atom_count(),
        2
    );
    for limits in [
        CmlParseLimits {
            max_atoms: 1,
            ..exact
        },
        CmlParseLimits {
            max_bonds: 0,
            ..exact
        },
        CmlParseLimits {
            max_xml_elements: 5,
            ..exact
        },
    ] {
        assert!(matches!(
            parse_cml_with_3d_with_limits(&source, &limits),
            Err(CmlError::ResourceLimit { .. })
        ));
    }
}
