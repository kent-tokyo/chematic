use serde_json::Value;
use std::io::Write;
use std::process::{Command, Stdio};

fn packet() -> Value {
    serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/moljson_graph_input/cases.json"
    )))
    .unwrap()
}

fn invoke(text: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .args([
            "convert",
            "--input-format",
            "moljson",
            "--output-format",
            "moljson",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(text.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

#[test]
fn cli_retains_supported_atom_and_bond_fields() {
    for case in packet()["valid"].as_array().unwrap() {
        let output = invoke(case["input"].as_str().unwrap());
        assert!(
            output.status.success(),
            "{}: {}",
            case["id"],
            String::from_utf8_lossy(&output.stderr)
        );
        let actual: Value = serde_json::from_slice(&output.stdout).unwrap();
        let expected = &case["expected"];
        assert_eq!(
            actual["atoms"].as_array().unwrap().len(),
            expected["atoms"].as_array().unwrap().len()
        );
        for (a, e) in actual["atoms"]
            .as_array()
            .unwrap()
            .iter()
            .zip(expected["atoms"].as_array().unwrap())
        {
            for field in ["element", "charge", "isotope", "hydrogens", "aromatic"] {
                assert_eq!(a[field], e[field], "{} {field}", case["id"]);
            }
        }
        assert_eq!(
            actual["bonds"].as_array().unwrap().len(),
            expected["bonds"].as_array().unwrap().len()
        );
        for (b, e) in actual["bonds"]
            .as_array()
            .unwrap()
            .iter()
            .zip(expected["bonds"].as_array().unwrap())
        {
            assert_eq!(b["order"].as_f64(), e["order"].as_f64(), "{}", case["id"]);
            assert_eq!(b["aromatic"], e["aromatic"]);
            for (field, index) in [("source_id", "source"), ("target_id", "target")] {
                let i = e[index].as_u64().unwrap() + 1;
                assert_eq!(b[field].as_str().unwrap(), format!("a{i}"));
            }
        }
    }
}

#[test]
fn cli_refuses_invalid_graphs_without_partial_molecule_output() {
    for case in packet()["invalid"].as_array().unwrap() {
        let output = invoke(case["input"].as_str().unwrap());
        assert_eq!(output.status.code(), Some(2), "{}", case["id"]);
        assert!(output.stdout.is_empty(), "{}", case["id"]);
        assert!(
            String::from_utf8_lossy(&output.stderr).starts_with("chematic: MolJSON:"),
            "{}",
            case["id"]
        );
    }
}
