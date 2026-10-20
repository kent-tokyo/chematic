use std::io::{ErrorKind, Write};
use std::process::{Command, Output, Stdio};

const SOURCE: &str = include_str!("../../../validation/mol2_preservation/amide.mol2");

fn run_formats(input: &str, source: &str, target: &str, args: &[&str]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .args([
            "convert-with-report",
            "--input-format",
            source,
            "--output-format",
            target,
        ])
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let write_result = child.stdin.take().unwrap().write_all(input.as_bytes());
    let output = child.wait_with_output().unwrap();
    match write_result {
        Ok(()) => {}
        // Invalid arguments can make the CLI exit before it reads stdin.
        Err(error) if error.kind() == ErrorKind::BrokenPipe && !output.status.success() => {}
        Err(error) => panic!("failed to send CLI input: {error}"),
    }
    output
}

fn run_to(input: &str, target: &str, args: &[&str]) -> Output {
    run_formats(input, "mol2", target, args)
}

fn run(input: &str, args: &[&str]) -> Output {
    run_to(input, "smiles", args)
}

#[test]
fn strict_report_is_machine_readable_and_signals_refusal() {
    let output = run(SOURCE, &[]);
    assert!(!output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["status"], "refused");
    assert!(report["output"].is_null());
    assert!(
        String::from_utf8(output.stderr)
            .unwrap()
            .contains("conversion refused")
    );
}

#[test]
fn explicit_loss_policy_returns_the_shared_rust_report() {
    let output = run(SOURCE, &["--policy", "allow_loss"]);
    assert!(output.status.success(), "{:?}", output.stderr);
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let expected = chematic_mol::convert_with_report(
        SOURCE,
        "mol2",
        "smiles",
        chematic_mol::ConversionPolicy::AllowLoss,
    )
    .unwrap();
    assert_eq!(report, serde_json::to_value(expected).unwrap());
}

#[test]
fn malformed_input_does_not_overwrite_existing_report_file() {
    let path = std::env::temp_dir().join(format!(
        "chematic-invalid-report-{}.json",
        std::process::id()
    ));
    std::fs::write(&path, "sentinel").unwrap();
    let output = run("invalid", &["--output", path.to_str().unwrap()]);
    let remaining = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(path).unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(remaining, "sentinel");
}

#[test]
fn source_retention_outputs_a_report_envelope() {
    let output = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .args([
            "convert-with-report",
            "--input-format",
            "mol2",
            "--output-format",
            "mol2",
            "--policy",
            "preserve_source",
            "--input",
        ])
        .arg(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/mol2_preservation/amide.mol2"
        ))
        .output()
        .unwrap();
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["output"], SOURCE);
    assert_eq!(report["verification_scope"], "exact_source_bytes");
}

#[test]
fn cjson_report_preserves_geometry_and_partial_charges() {
    let output = run_to(SOURCE, "cjson", &["--policy", "allow_loss"]);
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let candidate: serde_json::Value =
        serde_json::from_str(report["output"].as_str().unwrap()).unwrap();
    assert_eq!(candidate["atoms"]["coords"]["3d"][0], 0.123456789);
    assert!(
        report["observations"]
            .as_array()
            .unwrap()
            .iter()
            .any(|o| o["field"] == "coordinates" && o["status"] == "preserved")
    );
    assert!(
        report["observations"]
            .as_array()
            .unwrap()
            .iter()
            .any(|o| o["field"] == "partial_charges" && o["status"] == "preserved")
    );
    assert_eq!(run_to(SOURCE, "cjson", &[]).status.code(), Some(2));
}

#[test]
fn cjson_input_retention_refusal_and_malformed_file_safety() {
    let source = include_str!("../../../validation/cjson_preservation/charged.cjson");
    let output = run_formats(source, "cjson", "cjson", &[]);
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["output"], source);
    let strict = run_formats(source, "cjson", "smiles", &[]);
    assert_eq!(strict.status.code(), Some(2));
    let allowed = run_formats(source, "cjson", "smiles", &["--policy", "allow_loss"]);
    assert!(allowed.status.success());
    let report: serde_json::Value = serde_json::from_slice(&allowed.stdout).unwrap();
    assert_eq!(report["candidate_graph_verified"], true);
    let path = std::env::temp_dir().join(format!(
        "chematic-invalid-cjson-report-{}.json",
        std::process::id()
    ));
    std::fs::write(&path, "sentinel").unwrap();
    let bad = source.replace("[0, 1, -1, 0]", "[]");
    let output = run_formats(
        &bad,
        "cjson",
        "cjson",
        &["--output", path.to_str().unwrap()],
    );
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "sentinel");
    std::fs::remove_file(path).unwrap();
}

#[test]
fn cjson_to_mol2_preserves_declared_partial_charges_and_refuses_ambiguous_sets() {
    let source = include_str!("../../../validation/cjson_preservation/partial-charges.cjson");
    let allowed = run_formats(source, "cjson", "mol2", &["--policy", "allow_loss"]);
    assert!(allowed.status.success());
    let report: serde_json::Value = serde_json::from_slice(&allowed.stdout).unwrap();
    assert_eq!(report["candidate_graph_verified"], true);
    let molecule =
        chematic_mol::parse_mol2_source_record(report["output"].as_str().unwrap()).unwrap();
    assert_eq!(molecule.name(), Some("acetamide charges"));
    assert_eq!(molecule.atoms()[2].partial_charge, Some(-0.7));
    assert_eq!(
        run_formats(source, "cjson", "mol2", &[]).status.code(),
        Some(2)
    );
    let mut value: serde_json::Value = serde_json::from_str(source).unwrap();
    value["partialCharges"]["SECOND"] = serde_json::json!([0, 0, 0, 0]);
    let refused = run_formats(
        &value.to_string(),
        "cjson",
        "mol2",
        &["--policy", "allow_loss"],
    );
    assert_eq!(refused.status.code(), Some(2));
    let report: serde_json::Value = serde_json::from_slice(&refused.stdout).unwrap();
    assert!(report["output"].is_null());
}
