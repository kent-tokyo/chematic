//! Emit RDKit-compatible ECFP4 raw bitInfo provenance as JSONL.

use serde_json::json;
use std::io::{self, BufRead, Write};

fn raw_bit_info(mol: &chematic_core::Molecule) -> Vec<(u32, Vec<(u32, u32)>)> {
    let result = chematic_fp::rdkit_morgan_ecfp4_experimental(mol).expect("handled by caller");
    let mut values: Vec<_> = result.raw_bit_info.into_iter().collect();
    for (_, entries) in &mut values { entries.sort_unstable(); }
    values.sort_unstable_by_key(|(identifier, _)| *identifier);
    values
}

fn main() {
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for (index, line) in io::stdin().lock().lines().enumerate() {
        let smiles = line.expect("read corpus line");
        let smiles = smiles.trim();
        if smiles.is_empty() { continue; }
        let record = match chematic_smiles::parse(smiles) {
            Ok(mol) => match chematic_fp::rdkit_morgan_ecfp4_experimental(&mol) {
                Ok(_) => json!({"index": index, "smiles": smiles, "status": "ok", "raw_bit_info": raw_bit_info(&mol)}),
                Err(error) => json!({"index": index, "smiles": smiles, "status": "error", "error": error.to_string()}),
            },
            Err(error) => json!({"index": index, "smiles": smiles, "status": "error", "error": error.to_string()}),
        };
        serde_json::to_writer(&mut out, &record).expect("serialize raw bitInfo record");
        writeln!(out).expect("write raw bitInfo record");
    }
}
