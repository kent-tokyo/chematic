use serde_json::Value;
use std::io::Write;
use std::process::{Command, Output, Stdio};

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/v3000_atom_scalars/cases.json"
    ))
    .unwrap()
}
fn convert(input: &str, target: &str, args: &[&str]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
        .args([
            "convert",
            "--input-format",
            "v3000",
            "--output-format",
            target,
        ])
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
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

#[test]
fn fields_survive_cli_smiles_and_v3000_regeneration() {
    for c in packet()["cases"].as_array().unwrap() {
        let input = c["input"].as_str().unwrap();
        let output = convert(input, "smiles", &[]);
        assert!(output.status.success(), "{:?}", output.stderr);
        let text = String::from_utf8(output.stdout).unwrap();
        let mol = chematic_smiles::parse(text.trim()).unwrap();
        let atom = mol.atoms().next().unwrap().1;
        let e = &c["expected"];
        assert_eq!(i64::from(atom.charge), e["charge"].as_i64().unwrap());
        assert_eq!(
            u64::from(atom.isotope.unwrap_or(0)),
            e["isotope"].as_u64().unwrap()
        );
        assert_eq!(
            u64::from(atom.atom_map.unwrap_or(0)),
            e["atom_map"].as_u64().unwrap()
        );
        let output = convert(input, "v3000", &[]);
        assert!(output.status.success());
        let text = String::from_utf8(output.stdout).unwrap();
        let (after, _, _) = chematic_mol::parse_mol_v3000_with_coords(&text).unwrap();
        let after = after.atoms().next().unwrap().1;
        assert_eq!(
            (after.element, after.charge, after.isotope, after.atom_map),
            (atom.element, atom.charge, atom.isotope, atom.atom_map)
        );
        if c["scope"] == "quoted_compatibility" {
            assert!(text.contains("LABEL=\"memo CHG=2 MASS=13 end\""));
        }
    }
}

#[test]
fn bad_scalars_never_emit_a_partial_result_or_overwrite_a_file() {
    let dir = std::env::temp_dir().join(format!("chematic-v3000-scalars-{}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("protected.smi");
    std::fs::write(&path, "keep original\n").unwrap();
    for c in packet()["invalid"].as_array().unwrap() {
        for target in ["smiles", "v3000"] {
            let output = convert(c["input"].as_str().unwrap(), target, &[]);
            assert!(!output.status.success());
            assert!(output.stdout.is_empty());
            assert!(String::from_utf8_lossy(&output.stderr).contains(c["field"].as_str().unwrap()));
        }
        let output = convert(
            c["input"].as_str().unwrap(),
            "smiles",
            &["--output", path.to_str().unwrap()],
        );
        assert!(!output.status.success());
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "keep original\n");
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
