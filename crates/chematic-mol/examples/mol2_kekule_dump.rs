use std::io::{self, Read};

fn main() {
    let mut text = String::new();
    io::stdin().read_to_string(&mut text).unwrap();
    let inputs: Vec<String> = serde_json::from_str(&text).unwrap();
    let rows: Vec<_> = inputs
        .iter()
        .map(|input| match chematic_mol::kekulize_mol2_record(input) {
            Ok(output) => serde_json::json!({"output":output}),
            Err(error) => serde_json::json!({"error":error.to_string()}),
        })
        .collect();
    println!("{}", serde_json::to_string(&rows).unwrap());
}
