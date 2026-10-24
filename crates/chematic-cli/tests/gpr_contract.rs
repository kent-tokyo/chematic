use std::io::Write;
use std::process::{Command, Output, Stdio};
fn run(text: &str, args: &[&str]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .arg("gpr")
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
fn all_fields_roundtrip_through_inspect_and_edited_json() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/gpr_preservation/cases.json"
    ))
    .unwrap();
    for c in packet["valid"].as_array().unwrap() {
        let output = run(c["text"].as_str().unwrap(), &["--inspect"]);
        assert!(output.status.success(), "{} {:?}", c["id"], output);
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&output.stdout).unwrap(),
            c["expected"]
        );
        let text = run(
            std::str::from_utf8(&output.stdout).unwrap(),
            &["--from-json"],
        );
        assert!(text.status.success());
        let output = run(std::str::from_utf8(&text.stdout).unwrap(), &["--inspect"]);
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&output.stdout).unwrap(),
            c["expected"]
        );
    }
}
#[test]
fn malformed_records_exit_two_without_partial_output() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/gpr_preservation/cases.json"
    ))
    .unwrap();
    for c in packet["refused"].as_array().unwrap() {
        let out = run(c["text"].as_str().unwrap(), &[]);
        assert_eq!(out.status.code(), Some(2), "{}", c["id"]);
        assert!(out.stdout.is_empty());
        assert!(
            String::from_utf8(out.stderr)
                .unwrap()
                .starts_with("chematic: GPR line ")
        );
    }
    assert!(!run("", &["--inspect", "--from-json"]).status.success());
}
