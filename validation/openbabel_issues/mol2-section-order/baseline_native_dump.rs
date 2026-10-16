//! Baseline-only dump harness: copy into chematic-mol/examples at the pinned base.
use chematic_mol::{parse_mol2_record, write_mol2, write_mol2_record};
use serde_json::json;

fn main() {
    let path = std::env::args().nth(1).expect("fixture JSON path");
    let fixture: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut rows = Vec::new();
    for case in fixture["cases"].as_array().unwrap() {
        let output = if let Some(text) = case["smiles"].as_str() {
            write_mol2(&chematic_smiles::parse(text).unwrap(), &[])
        } else {
            write_mol2_record(&parse_mol2_record(case["record"].as_str().unwrap()).unwrap())
        };
        let parsed = parse_mol2_record(&output).unwrap();
        rows.push(json!({"id":case["id"],"output":output,
            "formal_charges":parsed.molecule.atoms().map(|(_,a)|a.charge).collect::<Vec<_>>(),
            "atoms":parsed.atoms.len(),"bonds":parsed.bonds.len()}));
    }
    println!(
        "{}",
        json!({"schema":"chematic.mol2-section-order-native.v1","rows":rows})
    );
}
