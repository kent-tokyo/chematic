//! Exercise the executable boundary, including stdin, stdout, files and errors.
use serde_json::{Value, json};
use std::io::Write;
use std::process::{Command, Output, Stdio};

fn invoke(args: &[&str], stdin: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
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
        .write_all(stdin.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}
fn json_output(args: &[&str], stdin: &str) -> Value {
    let out = invoke(args, stdin);
    assert!(
        out.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_slice(&out.stdout).unwrap()
}

#[test]
fn conversion_round_trips_all_advertised_molecular_formats() {
    for format in ["v3000", "cml", "cjson", "moljson", "cdxml"] {
        let out = invoke(
            &[
                "convert",
                "--input-format",
                "smiles",
                "--output-format",
                format,
            ],
            "CCO",
        );
        assert!(
            out.status.success(),
            "{format}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let encoded = String::from_utf8(out.stdout).unwrap();
        let out = invoke(
            &[
                "convert",
                "--input-format",
                format,
                "--output-format",
                "smiles",
            ],
            &encoded,
        );
        assert!(
            out.status.success(),
            "{format}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let decoded = String::from_utf8(out.stdout).unwrap();
        let report = json_output(&["parse", decoded.trim()], "");
        assert_eq!(report["formula"], "C2H6O", "{format}");
        assert_eq!(report["bonds"], 2, "{format}");
    }
}

#[test]
fn molecule_and_reaction_commands_emit_inspectable_json() {
    for (args, path, expected) in [
        (vec!["parse", "CCO"], "/atoms", json!(3)),
        (vec!["descriptors", "CCO"], "/formula", json!("C2H6O")),
        (vec!["cxsmiles", "CCO"], "/atom_count", json!(3)),
        (
            vec!["fingerprint", "CCO", "--algorithm", "ecfp6"],
            "/algorithm",
            json!("ecfp6"),
        ),
        (vec!["similarity", "CCO", "CCO"], "/similarity", json!(1.)),
        (vec!["substructure", "CCO", "[O]"], "/matches", json!([[2]])),
        (
            vec!["standardize", "CCO.[Na+]", "--largest-fragment-only"],
            "/output/atoms",
            json!(3),
        ),
        (vec!["reaction", "CCO>>CC=O"], "/products", json!(1)),
        (
            vec!["reaction-match", "CCO>>CC=O", "CO>>C=O"],
            "/matched",
            json!(true),
        ),
        (
            vec!["reaction-balance", "CCO>>CCO"],
            "/balanced",
            json!(true),
        ),
        (
            vec!["reaction-fingerprint", "CCO>>CC=O", "--mode", "or"],
            "/mode",
            json!("or"),
        ),
        (
            vec!["reaction-similarity", "CCO>>CC=O", "CCO>>CC=O"],
            "/similarity",
            json!(1.),
        ),
    ] {
        let output = json_output(&args, "");
        assert_eq!(
            output.pointer(path).unwrap(),
            &expected,
            "{args:?}: {output}"
        );
    }
    let report = json_output(&["report", "CCO"], "");
    assert!(report.is_object());
    assert!(report.to_string().contains("C2H6O"));
}

#[test]
fn batch_commands_keep_valid_records_and_input_order() {
    for (command, input) in [
        ("batch-descriptors", "CCO\ninvalid?\n"),
        ("batch-canonicalize", "CCO\ninvalid?\n"),
        ("batch-fingerprints", "CCO\ninvalid?\n"),
        ("batch-standardize", "CCO\ninvalid?\n"),
        ("batch-similarity", "CCO\tCCO\ninvalid?\tCCO\n"),
        ("batch-substructure", "CCO\t[O]\ninvalid?\t[O]\n"),
        ("batch-reactions", "CCO>>CC=O\ninvalid?\n"),
    ] {
        let output = json_output(&[command], input);
        assert_eq!(output["valid_count"], 1, "{command}: {output}");
        assert_eq!(output["error_count"], 1, "{command}: {output}");
        assert_eq!(output["records"][0]["input_index"], 0);
        assert_eq!(output["records"][1]["input_index"], 1);
        assert!(output["records"][1]["error"].is_string());
    }
    let output = json_output(&["batch-report"], "CCO\n");
    assert!(output.is_object());
    assert!(output.to_string().contains("C2H6O"));
}

#[test]
fn conversion_supports_files_and_errors_use_stderr() {
    let dir = std::env::temp_dir().join(format!("chematic-cli-contracts-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.smi");
    let output = dir.join("output.mol");
    std::fs::write(&input, "CCO").unwrap();
    let result = invoke(
        &[
            "convert",
            "--input-format",
            "smiles",
            "--output-format",
            "mol",
            "--input",
            input.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ],
        "",
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(result.stdout.is_empty());
    let mol = std::fs::read_to_string(&output).unwrap();
    let back = invoke(
        &[
            "convert",
            "--input-format",
            "mol",
            "--output-format",
            "smiles",
        ],
        &mol,
    );
    assert!(back.status.success());
    let smiles = String::from_utf8(back.stdout).unwrap();
    let parsed = json_output(&["parse", smiles.trim()], "");
    assert_eq!(parsed["formula"], "C2H6O");
    assert_eq!(parsed["atoms"], 3);
    assert_eq!(parsed["bonds"], 2);
    for args in [
        vec!["parse", "invalid?"],
        vec!["fingerprint", "CCO", "--algorithm", "unknown"],
        vec!["reaction-fingerprint", "CCO>>CC=O", "--mode", "unknown"],
        vec!["batch-descriptors", "--max-records", "0"],
        vec![
            "convert",
            "--input-format",
            "unknown",
            "--output-format",
            "smiles",
        ],
        vec![
            "convert",
            "--input-format",
            "smiles",
            "--output-format",
            "mol",
            "--input",
            dir.to_str().unwrap(),
        ],
    ] {
        let result = invoke(&args, "CCO\n");
        assert!(!result.status.success(), "{args:?}");
        assert!(result.stdout.is_empty());
        assert!(!result.stderr.is_empty());
    }
    std::fs::remove_dir_all(dir).unwrap();
}
