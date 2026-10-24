use std::io::Write;
use std::process::{Command, Stdio};

fn convert(source: &str, path: Option<&std::path::Path>) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_chematic"));
    command.args([
        "convert",
        "--input-format",
        "cjson",
        "--output-format",
        "cml",
    ]);
    if let Some(path) = path {
        command.arg("--output").arg(path);
    }
    let mut child = command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(source.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

#[test]
fn common_command_uses_the_same_checked_projection() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/cjson_cml_projection/cases.json"
    ))
    .unwrap();
    for case in packet["valid"].as_array().unwrap() {
        let source = case["input"].as_str().unwrap();
        let output = convert(source, None);
        assert!(
            output.status.success(),
            "{}: {}",
            case["id"],
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            output.stdout,
            chematic_mol::project_cjson_to_cml(source)
                .unwrap()
                .as_bytes()
        );
    }
}

#[test]
fn refused_source_never_emits_output_or_overwrites_a_file() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/cjson_cml_projection/cases.json"
    ))
    .unwrap();
    let path = std::env::temp_dir().join(format!(
        "chematic-cjson-cml-{}.sentinel",
        std::process::id()
    ));
    std::fs::write(&path, "keep").unwrap();
    for case in packet["invalid"].as_array().unwrap() {
        for target in [None, Some(path.as_path())] {
            let output = convert(case["input"].as_str().unwrap(), target);
            assert!(!output.status.success(), "{}", case["id"]);
            assert!(output.stdout.is_empty());
            assert!(!output.stderr.is_empty());
        }
    }
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "keep");
    std::fs::remove_file(path).unwrap();
}
