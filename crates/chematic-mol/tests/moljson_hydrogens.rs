use chematic_core::{AtomIdx, implicit_hcount};
use chematic_mol::{MolJsonError, parse_moljson, write_moljson};
use serde_json::{Value, json};

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/moljson_hydrogens/cases.json"
    ))
    .unwrap()
}

#[test]
fn authored_non_node_counts_and_graph_identity_survive_three_round_trips() {
    for case in packet()["cases"].as_array().unwrap() {
        let expected = chematic_smiles::parse(case["expected_smiles"].as_str().unwrap()).unwrap();
        assert_eq!(
            chematic_smiles::canonical_smiles(&parse_moljson(&write_moljson(&expected)).unwrap()),
            chematic_smiles::canonical_smiles(&expected),
        );
        let mut mol = parse_moljson(case["input"].as_str().unwrap()).unwrap();
        for _ in 0..3 {
            assert_eq!(
                chematic_smiles::canonical_smiles(&mol),
                chematic_smiles::canonical_smiles(&expected),
                "{}",
                case["id"]
            );
            for (i, count) in case["expected_hydrogens"]
                .as_array()
                .unwrap()
                .iter()
                .enumerate()
            {
                assert_eq!(
                    u64::from(implicit_hcount(&mol, AtomIdx(i as u32))),
                    count.as_u64().unwrap(),
                    "{} atom {i}",
                    case["id"]
                );
            }
            mol = parse_moljson(&write_moljson(&mol)).unwrap();
        }
    }
}

#[test]
fn invalid_present_counts_fail_without_defaulting_or_wrapping() {
    for value in packet()["invalid_hydrogens"].as_array().unwrap() {
        let input = json!({"atoms":[{"id":"a1","element":"C","hydrogens":value}],"bonds":[]});
        assert!(
            matches!(parse_moljson(&input.to_string()), Err(MolJsonError::InvalidJson(detail)) if detail.contains("hydrogens")),
            "{value}"
        );
    }
}

#[test]
fn missing_count_keeps_inference_and_storage_boundary_is_exact() {
    for (field, expected) in [(None, None), (Some(0), Some(0)), (Some(255), Some(255))] {
        let mut atom = json!({"id":"a1","element":"C"});
        if let Some(value) = field {
            atom["hydrogens"] = json!(value);
        }
        let mol = parse_moljson(&json!({"atoms":[atom],"bonds":[]}).to_string()).unwrap();
        assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, expected);
    }
}
