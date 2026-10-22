use serde_json::Value;
use std::io::Write;
use std::process::{Command, Stdio};
fn run(input: &str, target: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .args([
            "convert",
            "--input-format",
            "cjson",
            "--output-format",
            target,
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
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}
#[test]
fn cjson_common_conversion_refuses_invalid_declared_data_without_stdout() {
    let packet: Value =
        serde_json::from_str(include_str!("../../../validation/cjson_checked/cases.json")).unwrap();
    for case in packet["invalid"].as_array().unwrap() {
        for target in [
            "smiles",
            "mol",
            "mol_v3000",
            "mol2",
            "cml",
            "cjson",
            "moljson",
        ] {
            let result = run(case["input"].as_str().unwrap(), target);
            assert!(!result.status.success(), "{} {target}", case["id"]);
            assert!(result.stdout.is_empty());
            assert!(!result.stderr.is_empty());
        }
    }
}
#[test]
fn cjson_common_conversion_still_accepts_complete_and_compatibility_inputs() {
    let packet: Value =
        serde_json::from_str(include_str!("../../../validation/cjson_checked/cases.json")).unwrap();
    for group in ["valid", "compatibility"] {
        for case in packet[group].as_array().unwrap() {
            let result = run(case["input"].as_str().unwrap(), "cjson");
            assert!(result.status.success(), "{}", case["id"]);
            let output: Value = serde_json::from_slice(&result.stdout).unwrap();
            assert_eq!(
                output["atoms"]["elements"]["number"],
                case["expected"]["numbers"]
            );
        }
    }
}
