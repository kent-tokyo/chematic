//! Shared adversarial parser contract.
//!
//! This is intentionally driven by the same fixture consumed by the Python
//! and Node-hosted WASM tests. Every common topology parser must reject the
//! malformed input instead of panicking, accepting a partial molecule, or
//! silently falling back to an empty result.

use serde_json::Value;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cross_binding_contract.json"
));

#[test]
fn every_common_topology_parser_rejects_shared_adversarial_inputs() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let cases = document["adversarial"]
        .as_array()
        .expect("adversarial array");
    assert_eq!(cases.len(), 8);

    for case in cases {
        let id = case["id"].as_str().unwrap();
        let format = case["format"].as_str().unwrap();
        let input = case["input"].as_str().unwrap();
        let accepted = match format {
            "smiles" => chematic_smiles::parse(input).is_ok(),
            "mol" => chematic_mol::parse_mol(input).is_ok(),
            "mol_v3000" => chematic_mol::parse_mol_v3000(input).is_ok(),
            "mol2" => chematic_mol::parse_mol2(input).is_ok(),
            "cml" => chematic_mol::parse_cml(input).is_ok(),
            "cjson" => chematic_mol::parse_cjson(input).is_ok(),
            "moljson" => chematic_mol::parse_moljson(input).is_ok(),
            "cdxml" => chematic_mol::parse_cdxml(input).is_ok(),
            other => panic!("unregistered adversarial parser format: {other}"),
        };
        assert!(!accepted, "adversarial fixture was accepted: {id}");
    }
}

#[test]
fn extxyz_binding_contract_matches_shared_fixture() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["extxyz_contract"];
    assert_eq!(contract["schema_version"], 1);
    let frame = chematic_mol::parse_extxyz(contract["input"].as_str().unwrap()).unwrap();
    let expected = &contract["expected"];

    let expected_coords: Vec<Vec<f64>> =
        serde_json::from_value(expected["coords"].clone()).unwrap();
    let expected_coords = expected_coords
        .into_iter()
        .map(|row| (row[0], row[1], row[2]))
        .collect::<Vec<_>>();
    assert_eq!(frame.coords(), expected_coords);
    assert_eq!(
        frame.lattice.map(|lattice| lattice.to_vec()),
        Some(serde_json::from_value(expected["lattice"].clone()).unwrap())
    );
    let expected_info: Vec<(String, String)> = expected["info"]
        .as_object()
        .unwrap()
        .iter()
        .map(|(key, value)| (key.clone(), value.as_str().unwrap().to_string()))
        .collect();
    assert_eq!(frame.info, expected_info);

    let forces = frame
        .properties
        .iter()
        .find(|property| property.name == "forces")
        .unwrap();
    assert_eq!(forces.values.len(), 3);
    assert_eq!(
        forces.values[0],
        vec![
            chematic_mol::XyzValue::Real(0.1),
            chematic_mol::XyzValue::Real(0.0),
            chematic_mol::XyzValue::Real(0.0)
        ]
    );
    let tags = frame
        .properties
        .iter()
        .find(|property| property.name == "tag")
        .unwrap();
    assert_eq!(
        tags.values,
        vec![
            vec![chematic_mol::XyzValue::Integer(1)],
            vec![chematic_mol::XyzValue::Integer(2)],
            vec![chematic_mol::XyzValue::Integer(2)],
        ]
    );
}

#[test]
fn extxyz_writer_binding_contract_roundtrips_shared_fixture() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["extxyz_contract"];
    let frame = chematic_mol::parse_extxyz(contract["input"].as_str().unwrap()).unwrap();
    let written = chematic_mol::write_extxyz(&frame).expect("extxyz writer");
    let reparsed = chematic_mol::parse_extxyz(&written).expect("reparse writer output");
    assert_eq!(reparsed.coords(), frame.coords());
    assert_eq!(reparsed.lattice, frame.lattice);
    assert_eq!(reparsed.properties, frame.properties);
    assert_eq!(reparsed.info, frame.info);
}

#[test]
fn semantic_expansion_binding_contract_matches_shared_fixture() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["semantic_expansion_contract"];
    assert_eq!(contract["schema_version"], 1);

    for case in contract["cases"].as_array().unwrap() {
        let model = chematic_mol::SemanticModel::from_json(&case["model"]).unwrap();
        let selected = model.apply_json_command(&case["command"]).unwrap();
        if let Some(expected) = case.get("expected_selected_alternative") {
            assert_eq!(
                selected.to_json()["r_groups"][0]["selected_alternative"],
                *expected
            );
        }
        if let Some(expected) = case.get("expected_repeat_count") {
            assert_eq!(
                selected.to_json()["polymer_units"][0]["repeat_count"],
                *expected
            );
        }
        let base = chematic_smiles::parse(case["base_smiles"].as_str().unwrap()).unwrap();
        let expanded = selected.expand(&base).unwrap();
        assert_eq!(
            expanded.to_json()["source_to_expanded"],
            case["expected_source_to_expanded"],
            "{}",
            case["id"].as_str().unwrap()
        );
    }
}
