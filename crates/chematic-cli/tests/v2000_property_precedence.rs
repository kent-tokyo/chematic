use chematic_core::Molecule;
use serde_json::{Value, json};
use std::io::Write;
use std::process::{Command, Stdio};

fn fields(mol: &Molecule) -> Value {
    let mut atoms: Vec<_> = mol
        .atoms()
        .map(|(_, a)| {
            [
                i64::from(a.element.atomic_number()),
                i64::from(a.charge),
                i64::from(a.isotope.unwrap_or(0)),
                i64::from(a.atom_map.unwrap_or(0)),
            ]
        })
        .collect();
    atoms.sort();
    let mut bonds: Vec<_> = mol
        .bonds()
        .map(|(_, b)| {
            let a = mol.atom(b.atom1).element.atomic_number();
            let c = mol.atom(b.atom2).element.atomic_number();
            [a.min(c), a.max(c), b.order.order_value().unwrap() as u8]
        })
        .collect();
    bonds.sort();
    json!({"atoms": atoms, "bonds": bonds})
}

#[test]
fn cli_smiles_obeys_global_charge_and_isotope_precedence() {
    let packet: Value = serde_json::from_str(include_str!(
        "../../../validation/v2000_property_precedence/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let mut child = Command::new(env!("CARGO_BIN_EXE_chematic"))
            .args([
                "convert",
                "--input-format",
                "mol",
                "--output-format",
                "smiles",
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
        let output = child.wait_with_output().unwrap();
        assert!(
            output.status.success(),
            "{}: {:?}",
            case["id"],
            output.stderr
        );
        let text = String::from_utf8(output.stdout).unwrap();
        let after = chematic_smiles::parse(text.trim()).unwrap();
        let expected = &case["expected"];
        let mut atoms: Vec<_> = (0..expected["numbers"].as_array().unwrap().len())
            .map(|i| {
                ["numbers", "charges", "isotopes", "atom_maps"]
                    .map(|k| expected[k][i].as_i64().unwrap())
            })
            .collect();
        atoms.sort();
        let mut bonds: Vec<_> = expected["bonds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|b| {
                let a = expected["numbers"][b[0].as_u64().unwrap() as usize]
                    .as_u64()
                    .unwrap();
                let c = expected["numbers"][b[1].as_u64().unwrap() as usize]
                    .as_u64()
                    .unwrap();
                [a.min(c), a.max(c), b[2].as_u64().unwrap()]
            })
            .collect();
        bonds.sort();
        assert_eq!(
            fields(&after),
            json!({"atoms": atoms, "bonds": bonds}),
            "{}",
            case["id"]
        );
    }
}
