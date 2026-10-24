use serde_json::Value;
use std::io::Write;
use std::process::{Command, Output, Stdio};
fn run(text: &str, args: &[&str]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .arg("molden")
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
fn shared_geometry_vibration_records_are_inspected_edited_and_written() {
    let p: Value = serde_json::from_str(include_str!(
        "../../../validation/molden_vibrations/cases.json"
    ))
    .unwrap();
    for c in p["cases"].as_array().unwrap() {
        let out = run(c["input"].as_str().unwrap(), &["--inspect"]);
        assert!(out.status.success(), "{}: {:?}", c["id"], out.stderr);
        let mut record: Value = serde_json::from_slice(&out.stdout).unwrap();
        assert_eq!(record, c["expected"]);
        if record["frequencies"].is_array() {
            record["frequencies"][0] = (-123.45678901234567).into();
        }
        let text = run(&record.to_string(), &["--from-json"]);
        assert!(text.status.success(), "{:?}", text.stderr);
        let back = run(std::str::from_utf8(&text.stdout).unwrap(), &["--inspect"]);
        assert_eq!(
            serde_json::from_slice::<Value>(&back.stdout).unwrap(),
            record
        );
    }
}
#[test]
fn invalid_input_never_produces_partial_stdout() {
    let p: Value = serde_json::from_str(include_str!(
        "../../../validation/molden_vibrations/cases.json"
    ))
    .unwrap();
    for c in p["refusals"].as_array().unwrap() {
        let out = run(c["input"].as_str().unwrap(), &[]);
        assert!(!out.status.success(), "{}", c["id"]);
        assert!(out.stdout.is_empty());
    }
    assert!(!run("{}", &["--inspect", "--from-json"]).status.success());
    assert!(
        !run(
            "[Molden Format]\n[Atoms] AU\nH 1 1 0 0 0\n",
            &["--max-input-bytes", "1"]
        )
        .status
        .success()
    );
}
