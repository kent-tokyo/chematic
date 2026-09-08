//! Emit native MACCS 166-bit keys as JSONL for the cross-binding parity lane.

use chematic_fp::maccs;
use serde_json::json;
use std::io::{self, BufRead, Write};

fn bytes_hex(mol: &chematic_core::Molecule) -> String {
    let fp = maccs(mol);
    (0..21usize)
        .map(|byte_idx| {
            let mut byte = 0u8;
            for bit in 0..8usize {
                let index = byte_idx * 8 + bit;
                if index < 166 && fp.get(index) {
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
            Ok(mol) => json!({"index": index, "smiles": smiles, "status": "ok", "maccs_hex": bytes_hex(&mol)}),
            Err(error) => json!({"index": index, "smiles": smiles, "status": "error", "error": error.to_string()}),
        };
        serde_json::to_writer(&mut out, &record).expect("serialize MACCS record");
        writeln!(out).expect("write MACCS record");
    }
}
