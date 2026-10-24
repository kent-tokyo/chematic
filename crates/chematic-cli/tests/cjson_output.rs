use std::io::Write;
use std::process::{Command, Stdio};

fn convert(source: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .args([
            "convert",
            "--input-format",
            "smiles",
            "--output-format",
            "cjson",
        ])
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
fn unsupported_graph_meaning_never_produces_success_or_partial_stdout() {
    for source in [
        "F/C=C/F",
        "N[C@@](F)(Cl)Br",
        "[C]",
        "[CH4]",
        "c1ccccc1",
        "[CH4:7]",
        "*C",
        "[Pt@SP1](Cl)(Br)(I)F",
    ] {
        let output = convert(source);
        assert!(!output.status.success(), "{source}");
        assert!(output.stdout.is_empty(), "{source}");
        assert!(
            String::from_utf8(output.stderr)
                .unwrap()
                .contains("cannot preserve"),
            "{source}"
        );
    }
}

#[test]
fn supported_topology_still_converts() {
    let output = convert("CCO");
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(
        value["atoms"]["elements"]["number"],
        serde_json::json!([6, 6, 8])
    );
    assert_eq!(value["bonds"]["order"], serde_json::json!([1, 1]));
    assert!(value["atoms"].get("coords").is_none());
}
