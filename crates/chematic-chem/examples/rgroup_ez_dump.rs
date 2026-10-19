//! Emit labelled R-group decompositions for the pinned E/Z holdout pairs.

use std::io::{self, Read};

use chematic_chem::rgroup::rgroup_decompose_labeled;
use chematic_smiles::parse;
use serde_json::{Value, json};

fn main() {
    let mut input = String::new();
    io::stdin()
        .read_to_string(&mut input)
        .expect("read fixture");
    let fixture: Value = serde_json::from_str(&input).expect("parse fixture");
    let core = fixture["core"].as_str().expect("core string");
    let substituents = fixture["substituents"]
        .as_array()
        .expect("substituent array");

    for (left_index, left) in substituents.iter().enumerate() {
        let left = left.as_str().expect("left substituent");
        for right in &substituents[left_index + 1..] {
            let right = right.as_str().expect("right substituent");
            let inputs = [
                format!("{left}/C=C/{right}"),
                format!("{left}/C=C\\{right}"),
            ];
            let molecules = inputs
                .iter()
                .map(|smiles| parse(smiles).expect("parse input molecule"))
                .collect::<Vec<_>>();
            let refs = molecules.iter().collect::<Vec<_>>();
            let rows = rgroup_decompose_labeled(core, &refs).expect("decompose pair");
            let rows = rows
                .into_iter()
                .map(|row| {
                    row.map(|row| json!({"Core": row.core_smiles, "R": row.r_groups}))
                        .expect("matched input")
                })
                .collect::<Vec<_>>();
            println!(
                "{}",
                json!({"left": left, "right": right, "inputs": inputs, "rows": rows})
            );
        }
    }
}
