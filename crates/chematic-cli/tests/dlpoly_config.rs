use serde_json::Value;
use std::io::Write;
use std::process::{Command, Output, Stdio};
fn run(text: &str, args: &[&str]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .arg("dlpoly-config")
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
fn inspect_and_edit_every_shared_record() {
    let p: Value =
        serde_json::from_str(include_str!("../../../validation/dlpoly_config/cases.json")).unwrap();
    for c in p["cases"].as_array().unwrap() {
        let result = run(c["input"].as_str().unwrap(), &["--inspect"]);
        assert!(result.status.success(), "{}: {:?}", c["id"], result.stderr);
        let mut record: Value = serde_json::from_slice(&result.stdout).unwrap();
        assert_eq!(record, c["expected"]);
        record["atoms"][0]["position"][0] = 1.2345678901234567.into();
        let written = run(&record.to_string(), &["--from-json"]);
        assert!(written.status.success(), "{:?}", written.stderr);
        let back = run(
            std::str::from_utf8(&written.stdout).unwrap(),
            &["--inspect"],
        );
        assert_eq!(
            serde_json::from_slice::<Value>(&back.stdout).unwrap(),
            record
        );
    }
}
#[test]
fn refusal_produces_no_partial_stdout() {
    let p: Value =
        serde_json::from_str(include_str!("../../../validation/dlpoly_config/cases.json")).unwrap();
    for c in p["refusals"].as_array().unwrap() {
        let result = run(c["input"].as_str().unwrap(), &[]);
        assert!(!result.status.success(), "{}", c["id"]);
        assert!(result.stdout.is_empty());
    }
    assert!(!run("{}", &["--inspect", "--from-json"]).status.success());
    assert!(
        !run("title\n0 0\nX\n0 0 0\n", &["--max-input-bytes", "1"])
            .status
            .success()
    );
}
