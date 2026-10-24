use serde_json::Value;
use std::io::Write;
use std::process::{Command, Stdio};

#[test]
fn cml_layout_and_opaque_xml_do_not_change_cli_graph_projection() {
    let packet: Value = serde_json::from_str(include_str!(
        "../../../validation/cml_lexical_meaning/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
            .args([
                "convert",
                "--input-format",
                "cml",
                "--output-format",
                "moljson",
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
            .write_all(case["input"].as_str().unwrap().as_bytes())
            .unwrap();
        let result = child.wait_with_output().unwrap();
        assert!(
            result.status.success(),
            "{}: {:?}",
            case["id"],
            result.stderr
        );
        let value: Value = serde_json::from_slice(&result.stdout).unwrap();
        let atoms = value["atoms"].as_array().unwrap();
        assert_eq!(
            atoms.len(),
            case["expected"]["numbers"].as_array().unwrap().len()
        );
        for (i, atom) in atoms.iter().enumerate() {
            let number = match atom["element"].as_str().unwrap() {
                "C" => 6,
                "N" => 7,
                "O" => 8,
                other => panic!("{other}"),
            };
            assert_eq!(number, case["expected"]["numbers"][i]);
            assert_eq!(atom["charge"], case["expected"]["charges"][i]);
            assert_eq!(atom["isotope"], case["expected"]["isotopes"][i]);
            assert_eq!(atom["hydrogens"], case["expected"]["hydrogens"][i]);
        }
        let bonds = value["bonds"].as_array().unwrap();
        assert_eq!(
            bonds.len(),
            case["expected"]["bonds"].as_array().unwrap().len()
        );
        for (i, bond) in bonds.iter().enumerate() {
            let expected = &case["expected"]["bonds"][i];
            assert_eq!(
                bond["source_id"],
                atoms[expected[0].as_u64().unwrap() as usize]["id"]
            );
            assert_eq!(
                bond["target_id"],
                atoms[expected[1].as_u64().unwrap() as usize]["id"]
            );
            assert_eq!(
                bond["order"].as_f64().unwrap(),
                expected[2].as_f64().unwrap()
            );
        }
    }
}
