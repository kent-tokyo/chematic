use serde_json::Value;
use std::io::Write;
use std::process::{Command, Stdio};

#[test]
fn cdxml_cli_preserves_supported_bonds_and_refuses_unsupported_declarations() {
    let packet: Value = serde_json::from_str(include_str!(
        "../../../validation/cdxml_bond_meaning/cases.json"
    ))
    .unwrap();
    for c in packet["cases"].as_array().unwrap() {
        let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
            .args([
                "convert",
                "--input-format",
                "cdxml",
                "--output-format",
                "smiles",
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
            .write_all(c["input"].as_str().unwrap().as_bytes())
            .unwrap();
        let output = child.wait_with_output().unwrap();
        if c["group"] == "refusal" {
            assert!(!output.status.success(), "{}", c["id"]);
            assert!(output.stdout.is_empty());
            assert!(
                String::from_utf8(output.stderr)
                    .unwrap()
                    .contains("unsupported CDXML bond order")
            );
        } else {
            assert!(output.status.success(), "{}: {:?}", c["id"], output.stderr);
            let expected = chematic_smiles::canonical_smiles(
                &chematic_smiles::parse(c["expected"]["smiles"].as_str().unwrap()).unwrap(),
            );
            assert_eq!(
                String::from_utf8(output.stdout).unwrap().trim(),
                expected,
                "{}",
                c["id"]
            );
        }
    }
}
