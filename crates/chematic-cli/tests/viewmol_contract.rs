use std::io::Write;
use std::process::{Command, Output, Stdio};
fn run(text: &str, args: &[&str]) -> Output {
    let mut c = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .arg("viewmol")
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    c.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();
    c.wait_with_output().unwrap()
}
fn packet() -> serde_json::Value {
    serde_json::from_str(include_str!(
        "../../../validation/viewmol_preservation/cases.json"
    ))
    .unwrap()
}
#[test]
fn inspect_regenerate_and_project_through_actual_cli() {
    for c in packet()["valid"].as_array().unwrap() {
        let out = run(c["text"].as_str().unwrap(), &["--inspect"]);
        assert!(out.status.success(), "{} {:?}", c["id"], out);
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&out.stdout).unwrap(),
            c["expected"]
        );
        let text = run(std::str::from_utf8(&out.stdout).unwrap(), &["--from-json"]);
        assert!(text.status.success());
        let after = run(std::str::from_utf8(&text.stdout).unwrap(), &["--inspect"]);
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&after.stdout).unwrap(),
            c["expected"]
        );
        let p = run(
            std::str::from_utf8(&out.stdout).unwrap(),
            &["--from-json", "--angstrom"],
        );
        assert!(p.status.success());
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&p.stdout).unwrap(),
            c["expected_angstrom"]
        );
    }
}
#[test]
fn invalid_input_and_unsafe_projection_exit_two_without_partial_stdout() {
    for c in packet()["refused"].as_array().unwrap() {
        let out = run(c["text"].as_str().unwrap(), &[]);
        assert_eq!(out.status.code(), Some(2), "{}", c["id"]);
        assert!(out.stdout.is_empty());
        assert!(
            String::from_utf8(out.stderr)
                .unwrap()
                .starts_with("chematic: ViewMol line ")
        );
    }
    for c in packet()["projection_refused"].as_array().unwrap() {
        assert!(run(c["text"].as_str().unwrap(), &[]).status.success());
        let out = run(c["text"].as_str().unwrap(), &["--angstrom"]);
        assert_eq!(out.status.code(), Some(2));
        assert!(out.stdout.is_empty());
    }
    assert!(!run("", &["--inspect", "--angstrom"]).status.success());
    assert!(!run("", &["--inspect", "--from-json"]).status.success());
}
