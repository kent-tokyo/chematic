//! Actual editable/native writer output for the original section-order controls.
use chematic_mol::{
    Mol2SectionOrder, parse_mol2_record, write_mol2, write_mol2_record,
    write_mol2_record_with_section_order, write_mol2_with_section_order,
};
use serde_json::json;

fn main() {
    let path = std::env::args().nth(1).expect("fixture JSON path");
    let fixture: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut rows = Vec::new();
    for case in fixture["cases"].as_array().unwrap() {
        let default_output = if let Some(text) = case["smiles"].as_str() {
            write_mol2(&chematic_smiles::parse(text).unwrap(), &[])
        } else {
            write_mol2_record(&parse_mol2_record(case["record"].as_str().unwrap()).unwrap())
        };
        let output = if let Some(text) = case["smiles"].as_str() {
            write_mol2_with_section_order(
                &chematic_smiles::parse(text).unwrap(),
                &[],
                Mol2SectionOrder::AtomBondAdjacent,
            )
        } else {
            write_mol2_record_with_section_order(
                &parse_mol2_record(case["record"].as_str().unwrap()).unwrap(),
                Mol2SectionOrder::AtomBondAdjacent,
            )
        };
        let parsed = parse_mol2_record(&output).unwrap();
        rows.push(
            json!({"id":case["id"],"output":output, "default_output":default_output,
            "formal_charges":parsed.molecule.atoms().map(|(_,a)|a.charge).collect::<Vec<_>>(),
            "atoms":parsed.atoms.len(),"bonds":parsed.bonds.len()}),
        );
    }
    println!(
        "{}",
        json!({"schema":"chematic.mol2-section-order-native.v1","rows":rows})
    );
}
