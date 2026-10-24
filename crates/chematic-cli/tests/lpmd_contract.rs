use std::io::Write;
use std::process::{Command, Stdio};
fn run(args: &[&str], input: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .arg("lpmd")
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
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}
#[test]
fn packet_and_json_edit_roundtrip() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/lpmd_preservation/cases.json"
    ))
    .unwrap();
    for c in packet["valid"].as_array().unwrap() {
        let out = run(&["--inspect"], c["text"].as_str().unwrap());
        assert!(out.status.success(), "{:?}", out);
        let record: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
        assert_eq!(record, c["record"]);
        let out = run(&["--from-json"], &record.to_string());
        assert!(out.status.success());
        let out = run(&["--inspect"], std::str::from_utf8(&out.stdout).unwrap());
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&out.stdout).unwrap(),
            record
        );
    }
    for c in packet["refused"].as_array().unwrap() {
        let out = run(&[], c["text"].as_str().unwrap());
        assert_eq!(out.status.code(), Some(2));
        assert!(out.stdout.is_empty());
    }
    for c in packet["projection_refused"].as_array().unwrap() {
        let out = run(&["--angstrom"], c["text"].as_str().unwrap());
        assert_eq!(out.status.code(), Some(2));
        assert!(out.stdout.is_empty());
    }
}
#[test]
fn projection_and_limit_and_conflicting_modes() {
    let input = "LPMD 2.0 L\nHDR SYM X Y Z\n1\n4 1 0 -1 8 2 1 0 16\nC .125 .25 .5\n";
    let out = run(&["--angstrom"], input);
    assert!(out.status.success());
    assert_eq!(
        serde_json::from_slice::<serde_json::Value>(&out.stdout).unwrap(),
        serde_json::json!([[[0.75, 2.125, 8.5]]])
    );
    assert_eq!(
        run(&["--inspect", "--from-json"], input).status.code(),
        Some(2)
    );
    assert_eq!(
        run(&["--max-input-bytes", "10"], input).status.code(),
        Some(2)
    );
}
