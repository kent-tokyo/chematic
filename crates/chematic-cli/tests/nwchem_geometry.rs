use std::io::Write;
use std::process::{Command, Stdio};
fn invoke(args: &[&str], text: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .args(args)
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
fn nwchem_geometry_shared_packet_and_json_editing() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/nwchem_geometry/cases.json"
    ))
    .unwrap();
    for case in packet["valid"].as_array().unwrap() {
        let source = case["input"].as_str().unwrap();
        let output = invoke(&["nwchem-geometry", "--inspect"], source);
        assert!(output.status.success(), "{}: {:?}", case["id"], output);
        let record: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(record, case["expected"]);
        let output = invoke(&["nwchem-geometry", "--from-json"], &record.to_string());
        assert!(output.status.success());
        assert_eq!(
            chematic_mol::parse_nwchem_geometry(std::str::from_utf8(&output.stdout).unwrap())
                .unwrap(),
            chematic_mol::parse_nwchem_geometry(source).unwrap()
        );
    }
    for case in packet["invalid"].as_array().unwrap() {
        let output = invoke(&["nwchem-geometry"], case["input"].as_str().unwrap());
        assert!(!output.status.success(), "{}", case["id"]);
        assert!(output.stdout.is_empty());
    }
}
#[test]
fn nwchem_geometry_refuses_unknown_json_and_input_limit() {
    assert!(
        !invoke(&["nwchem-geometry", "--from-json"], "{\"unknown\":1}")
            .status
            .success()
    );
    assert!(
        !invoke(
            &["nwchem-geometry", "--max-input-bytes", "1"],
            "geometry\nend\n"
        )
        .status
        .success()
    );
}
