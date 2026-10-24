use serde_json::Value;
use std::io::Write;
use std::process::{Command, Stdio};

#[test]
fn cli_mol_and_sdf_alias_keep_isotopes_and_neighboring_atom_fields() {
    let packet: Value = serde_json::from_str(include_str!(
        "../../../validation/v2000_isotope_output/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let before = chematic_smiles::parse(smiles).unwrap();
        for target in ["mol", "sdf"] {
            let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
                .args([
                    "convert",
                    "--input-format",
                    "smiles",
                    "--output-format",
                    target,
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
                .write_all(smiles.as_bytes())
                .unwrap();
            let output = child.wait_with_output().unwrap();
            assert!(
                output.status.success(),
                "{}: {:?}",
                case["id"],
                output.stderr
            );
            let text = String::from_utf8(output.stdout).unwrap();
            let (after, _) = chematic_mol::parse_mol(&text).unwrap();
            assert_eq!(after.atom_count(), before.atom_count(), "{}", case["id"]);
            assert_eq!(after.bond_count(), before.bond_count(), "{}", case["id"]);
            for ((_, a), (_, b)) in after.atoms().zip(before.atoms()) {
                assert_eq!(
                    (a.element, a.charge, a.isotope, a.atom_map),
                    (b.element, b.charge, b.isotope, b.atom_map),
                    "{}",
                    case["id"]
                );
            }
        }
    }
}
