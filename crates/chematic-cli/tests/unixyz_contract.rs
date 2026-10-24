use std::io::Write;
use std::process::{Command, Output, Stdio};
fn run(text: &str, args: &[&str]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .arg("unixyz")
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
fn packet_roundtrips_and_inspects_through_actual_cli() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/unixyz_preservation/cases.json"
    ))
    .unwrap();
    for case in packet["valid"].as_array().unwrap() {
        let inspect = run(case["text"].as_str().unwrap(), &["--inspect"]);
        assert!(inspect.status.success(), "{} {:?}", case["id"], inspect);
        let record: serde_json::Value = serde_json::from_slice(&inspect.stdout).unwrap();
        assert_eq!(record, case["expected"]);
        let rewritten = run(
            std::str::from_utf8(&inspect.stdout).unwrap(),
            &["--from-json"],
        );
        assert!(rewritten.status.success());
        let result = run(
            std::str::from_utf8(&rewritten.stdout).unwrap(),
            &["--inspect"],
        );
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&result.stdout).unwrap(),
            record
        );
    }
}
#[test]
fn refusals_are_exit_two_without_partial_stdout() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/unixyz_preservation/cases.json"
    ))
    .unwrap();
    for case in packet["refused"].as_array().unwrap() {
        let output = run(case["text"].as_str().unwrap(), &[]);
        assert_eq!(output.status.code(), Some(2), "{}", case["id"]);
        assert!(output.stdout.is_empty());
        assert!(
            std::str::from_utf8(&output.stderr)
                .unwrap()
                .starts_with("chematic: UniChem XYZ line ")
        );
    }
    assert!(!run("", &["--inspect", "--from-json"]).status.success());
}
