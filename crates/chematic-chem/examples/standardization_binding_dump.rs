//! Emit the explicit shared standardization profile as JSONL for parity.

use chematic_chem::{StandardizeOptions, standardize};
use serde_json::json;
use std::io::{self, BufRead, Write};

fn main() {
    let stdout = io::stdout();
    let mut out = stdout.lock();
    let options = StandardizeOptions {
        largest_fragment_only: true,
        ..Default::default()
    };
    for (index, line) in io::stdin().lock().lines().enumerate() {
        let smiles = line.expect("read corpus line");
        let smiles = smiles.trim();
        if smiles.is_empty() {
            continue;
        }
        let record = match chematic_smiles::parse(smiles) {
            Ok(mol) => {
                let standardized = standardize(&mol, &options);
                json!({"index": index, "smiles": smiles, "status": "ok", "standardized_smiles": chematic_smiles::canonical_smiles(&standardized)})
            }
            Err(error) => {
                json!({"index": index, "smiles": smiles, "status": "error", "error": error.to_string()})
            }
        };
        serde_json::to_writer(&mut out, &record).expect("serialize standardization record");
        writeln!(out).expect("write standardization record");
    }
}
