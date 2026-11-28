use std::io::Write;
use std::process::{Command, Stdio};

const WATER: &str =
    "%chk=a.chk\n# HF/6-31G\n\nwater\n\n0 1\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\nH 0.93 0.0 -0.24\n\n";

fn run(args: &[&str], text: &str) -> std::process::Output {
    let mut p = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .arg("gjf-document")
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
fn gjf_document_round_trips_source_text_exactly() {
    let out = run(&[], WATER);
    assert!(out.status.success(), "{out:?}");
    assert_eq!(std::str::from_utf8(&out.stdout).unwrap(), WATER);
}

#[test]
fn gjf_document_inspect_reports_each_job_geometry_and_source_range() {
    let out = run(&["--inspect"], WATER);
    assert!(out.status.success(), "{out:?}");
    let v: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(v["schema"], "chematic.gjf-document.v1");
    let job = &v["jobs"][0];
    assert_eq!(job["elements"], serde_json::json!(["O", "H", "H"]));
    assert_eq!(job["charge"], 0);
    assert_eq!(job["multiplicity"], 1);
    assert_eq!(job["coords"][1], serde_json::json!([0.0, 0.0, 0.96]));
    assert_eq!(job["source_start"], 0);
    assert_eq!(job["source_end"], WATER.len() as u64);
}

#[test]
fn gjf_document_enforces_the_job_limit() {
    let out = run(&["--max-jobs", "0"], WATER);
    assert_eq!(out.status.code(), Some(2));
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("exceeding limit 0"), "{err}");
}
