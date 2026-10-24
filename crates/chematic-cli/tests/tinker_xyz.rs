use serde_json::Value;
use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run(text: &str, args: &[&str]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .arg("tinker-xyz")
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
fn cli_inspects_edits_and_writes_all_shared_records() {
    let packet: Value =
        serde_json::from_str(include_str!("../../../validation/tinker_xyz/cases.json")).unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let result = run(case["input"].as_str().unwrap(), &["--inspect"]);
        assert!(
            result.status.success(),
            "{}: {:?}",
            case["id"],
            result.stderr
        );
        let mut record: Value = serde_json::from_slice(&result.stdout).unwrap();
        assert_eq!(record, case["expected"]);
        record["atoms"][0]["atom_type"] = 90001.into();
        let written = run(&record.to_string(), &["--from-json"]);
        assert!(written.status.success(), "{:?}", written.stderr);
        let inspected = run(
            std::str::from_utf8(&written.stdout).unwrap(),
            &["--inspect"],
        );
        assert_eq!(
            serde_json::from_slice::<Value>(&inspected.stdout).unwrap(),
            record
        );
    }
}

#[test]
fn cli_refuses_invalid_input_without_partial_output() {
    let packet: Value =
        serde_json::from_str(include_str!("../../../validation/tinker_xyz/cases.json")).unwrap();
    for case in packet["refusals"].as_array().unwrap() {
        let result = run(case["input"].as_str().unwrap(), &[]);
        assert!(!result.status.success(), "{}", case["id"]);
        assert!(result.stdout.is_empty(), "{}", case["id"]);
    }
    assert!(!run("{}", &["--inspect", "--from-json"]).status.success());
    assert!(
        !run("1\n1 C 0 0 0 1\n", &["--max-input-bytes", "1"])
            .status
            .success()
    );
    assert!(
        !run(
            "{\"title\":\"\",\"cell\":null,\"atoms\":[],\"charge\":1}",
            &["--from-json"]
        )
        .status
        .success()
    );
}
