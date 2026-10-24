//! Reproduce numeric output from all editable/generated native MOL2 writers.
use chematic_mol::mol2_source::write_mol2_source;
use chematic_mol::{
    Mol2SectionOrder, parse_mol2_record, write_mol2, write_mol2_record,
    write_mol2_record_with_section_order,
};
use serde_json::{Value, json};

fn main() {
    let path = std::env::args().nth(1).expect("fixture JSON path");
    let packet: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let rows: Vec<_> = packet["cases"].as_array().unwrap().iter().map(|case| {
        let record = parse_mol2_record(case["input"].as_str().unwrap()).unwrap();
        json!({
            "id": case["id"],
            "record": write_mol2_record(&record),
            "record_adjacent": write_mol2_record_with_section_order(&record, Mol2SectionOrder::AtomBondAdjacent),
            "graph": write_mol2(&record.molecule, &record.coords),
            "source_graph": write_mol2_source(&record.molecule, &record.coords),
        })
    }).collect();
    println!(
        "{}",
        json!({"schema":"chematic.mol2-numeric-output.v1", "rows":rows})
    );
}
