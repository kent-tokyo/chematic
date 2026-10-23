use std::io::Write;
use std::process::{Command, Stdio};

fn convert(source: &str, output: Option<&std::path::Path>) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_chematic"));
    command.args([
        "convert",
        "--input-format",
        ".CML",
        "--output-format",
        "cml",
    ]);
    if let Some(path) = output {
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
fn common_command_returns_exact_source_without_adding_a_newline() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/cml_source_preservation/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let source = case["source"].as_str().unwrap();
        for text in [source, source.trim_end_matches(['\r', '\n'])] {
            let output = convert(text, None);
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert_eq!(output.stdout, text.as_bytes());
        }
    }
}

#[test]
fn malformed_later_xml_does_not_emit_partial_text_or_overwrite_a_file() {
    let source = "<molecule xmlns=\"http://www.xml-cml.org/schema\"><atom/></molecule><bad/>";
    let path = std::env::temp_dir().join(format!(
        "chematic-cml-common-{}.sentinel",
        std::process::id()
    ));
    std::fs::write(&path, "keep existing file").unwrap();
    for target in [None, Some(path.as_path())] {
        let output = convert(source, target);
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(!output.stderr.is_empty());
    }
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        "keep existing file"
    );
    std::fs::remove_file(path).unwrap();
}
