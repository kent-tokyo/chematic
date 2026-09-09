//! Emit RDKit-compatible ECFP4 bytes as JSONL for cross-binding parity.

use serde_json::json;
use std::io::{self, BufRead, Write};

fn bytes_hex(mol: &chematic_core::Molecule) -> String {
    let result = chematic_fp::rdkit_morgan_ecfp4_experimental(mol).expect("handled by caller");
    (0..256usize)
        .map(|byte_idx| {
            let mut byte = 0u8;
            for bit in 0..8usize {
                if result.fingerprint.get(byte_idx * 8 + bit) {
                    byte |= 1 << bit;
                }
            }
            format!("{byte:02x}")
        })
        .collect()
}

fn main() {
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for (index, line) in io::stdin().lock().lines().enumerate() {
        let smiles = line.expect("read corpus line");
        let smiles = smiles.trim();
        if smiles.is_empty() {
            continue;
        }
        let record = match chematic_smiles::parse(smiles) {
            Ok(mol) => match chematic_fp::rdkit_morgan_ecfp4_experimental(&mol) {
                Ok(_) => {
                    json!({"index": index, "smiles": smiles, "status": "ok", "ecfp4_hex": bytes_hex(&mol)})
                }
                Err(error) => {
                    json!({"index": index, "smiles": smiles, "status": "error", "error": error.to_string()})
                }
            },
            Err(error) => {
                json!({"index": index, "smiles": smiles, "status": "error", "error": error.to_string()})
            }
        };
        serde_json::to_writer(&mut out, &record).expect("serialize RDKit ECFP4 record");
        writeln!(out).expect("write RDKit ECFP4 record");
    }
}
