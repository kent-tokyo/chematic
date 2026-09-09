//! Emit native topological-path fingerprint bytes as JSONL for binding parity.

use chematic_fp::{TopoPathConfig, topo_path};
use serde_json::json;
use std::io::{self, BufRead, Write};

fn bytes_hex(mol: &chematic_core::Molecule) -> String {
    let fp = topo_path(mol, &TopoPathConfig::default());
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
            Ok(mol) => {
                json!({"index": index, "smiles": smiles, "status": "ok", "topo_path_hex": bytes_hex(&mol)})
            }
            Err(error) => {
                json!({"index": index, "smiles": smiles, "status": "error", "error": error.to_string()})
            }
        };
        serde_json::to_writer(&mut out, &record).expect("serialize topo-path record");
        writeln!(out).expect("write topo-path record");
    }
}
