use std::io::Write;
use std::process::{Command, Stdio};
fn run(args: &[&str], text: &str) -> std::process::Output {
    let mut p = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .arg("gromos96")
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    p.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();
    p.wait_with_output().unwrap()
}
#[test]
fn all_frames_raw_fields_and_explicit_units_through_cli() {
    let p: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/gromos96_preservation/cases.json"
    ))
    .unwrap();
    for c in p["valid"].as_array().unwrap() {
        let out = run(
            &[
                "--inspect",
                "--length-unit",
                c["record"]["length_unit"].as_str().unwrap(),
            ],
            c["text"].as_str().unwrap(),
        );
        assert!(out.status.success(), "{:?}", out);
        let r: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
        assert_eq!(r, c["record"]);
        let out = run(&["--from-json"], &r.to_string());
        assert!(out.status.success());
        let back = run(
            &[
                "--inspect",
                "--length-unit",
                r["length_unit"].as_str().unwrap(),
            ],
            std::str::from_utf8(&out.stdout).unwrap(),
        );
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&back.stdout).unwrap(),
            r
        );
    }
    for c in p["refused"].as_array().unwrap() {
        let out = run(&[], c["text"].as_str().unwrap());
        assert_eq!(out.status.code(), Some(2));
        assert!(out.stdout.is_empty());
    }
    for c in p["projection_refused"].as_array().unwrap() {
        let out = run(
            &[
                "--angstrom",
                "--length-unit",
                c["record"]["length_unit"].as_str().unwrap(),
            ],
            c["text"].as_str().unwrap(),
        );
        assert_eq!(out.status.code(), Some(2));
        assert!(out.stdout.is_empty());
    }
}
#[test]
fn unknown_units_and_invalid_modes_never_overwrite_output() {
    let p: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/gromos96_preservation/cases.json"
    ))
    .unwrap();
    let text = p["valid"][0]["text"].as_str().unwrap();
    assert_eq!(run(&["--angstrom"], text).status.code(), Some(2));
    assert_eq!(
        run(&["--inspect", "--angstrom"], text).status.code(),
        Some(2)
    );
    assert_eq!(
        run(&["--from-json", "--length-unit", "nm"], text)
            .status
            .code(),
        Some(2)
    );
    let file = std::env::temp_dir().join(format!("chematic-g96-refusal-{}", std::process::id()));
    std::fs::write(&file, "keep").unwrap();
    let out = run(&["--angstrom", "-o", file.to_str().unwrap()], text);
    assert_eq!(out.status.code(), Some(2));
    assert_eq!(std::fs::read_to_string(&file).unwrap(), "keep");
    std::fs::remove_file(file).unwrap();
}
