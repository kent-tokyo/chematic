use std::io::Write;
use std::process::{Command, Stdio};
fn run(text: &str, extra: &[&str]) -> std::process::Output {
    let mut c = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .args([
            "convert",
            "--input-format",
            "cml",
            "--output-format",
            "cjson",
        ])
        .args(extra)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    c.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();
    c.wait_with_output().unwrap()
}
#[test]
fn common_cli_uses_the_checked_projector() {
    let p: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/cml_cjson_projection/cases.json"
    ))
    .unwrap();
    for c in p["valid"].as_array().unwrap() {
        let s = c["input"].as_str().unwrap();
        let out = run(s, &[]);
        assert!(
            out.status.success(),
            "{}: {}",
            c["id"],
            String::from_utf8_lossy(&out.stderr)
        );
        assert_eq!(
            String::from_utf8(out.stdout).unwrap().trim_end(),
            chematic_mol::project_cml_to_cjson(s).unwrap()
        );
    }
    for c in p["invalid"].as_array().unwrap() {
        let out = run(c["input"].as_str().unwrap(), &[]);
        assert!(!out.status.success(), "{}", c["id"]);
        assert!(out.stdout.is_empty());
    }
}
#[test]
fn refusal_does_not_overwrite_existing_output() {
    let file = std::env::temp_dir().join(format!("chematic-cml-cjson-{}.json", std::process::id()));
    std::fs::write(&file, "keep").unwrap();
    let out = run(
        "<molecule><atom id='a' elementType='C' hydrogenCount='0'/></molecule>",
        &["--output", file.to_str().unwrap()],
    );
    assert!(!out.status.success());
    assert!(out.stdout.is_empty());
    assert_eq!(std::fs::read_to_string(&file).unwrap(), "keep");
    std::fs::remove_file(file).unwrap();
}
