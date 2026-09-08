//! Emit native topological-torsion fingerprint bytes as JSONL for parity.

use chematic_fp::torsion_fp;
use serde_json::json;
use std::io::{self, BufRead, Write};

fn bytes_hex(mol: &chematic_core::Molecule) -> String {
    let fp = torsion_fp(mol);
    (0..256usize)
        .map(|byte_idx| {
            let mut byte = 0u8;
            for bit in 0..8usize {
                if fp.get(byte_idx * 8 + bit) {
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
            Ok(mol) => json!({"index": index, "smiles": smiles, "status": "ok", "torsion_hex": bytes_hex(&mol)}),
            Err(error) => json!({"index": index, "smiles": smiles, "status": "error", "error": error.to_string()}),
        };
        serde_json::to_writer(&mut out, &record).expect("serialize torsion record");
        writeln!(out).expect("write torsion record");
    }
}
