use std::io::{self, BufRead};

use chematic_chem::rgroup_decompose_labeled;
use chematic_smiles::parse;
use serde_json::{Value, json};

fn main() {
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let row: Value = serde_json::from_str(&line).unwrap();
        let core = row["core"].as_str().unwrap();
        let mols: Vec<_> = row["molecules"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| parse(value.as_str().unwrap()).unwrap())
            .collect();
        let refs: Vec<_> = mols.iter().collect();
        let result = rgroup_decompose_labeled(core, &refs).unwrap();
        let rows: Vec<_> = result
            .iter()
            .map(|row| {
                row.as_ref().map(|row| {
                    let mut value = serde_json::Map::new();
                    value.insert("Core".to_string(), json!(row.core_smiles));
                    for (label, group) in &row.r_groups {
                        value.insert(format!("R{label}"), json!(group));
                    }
                    Value::Object(value)
                })
            })
            .collect();
        println!("{}", json!(rows));
    }
}
