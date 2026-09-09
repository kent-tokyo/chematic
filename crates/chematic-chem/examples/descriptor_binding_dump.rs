//! Emit the core descriptor contract as JSONL for the cross-binding parity lane.
//!
//! The input is one SMILES per line.  This is deliberately a small binding
//! contract (not an external accuracy oracle): Rust, Python, and Node/WASM
//! must agree on parse status and the values produced by the same chematic
//! implementation.

use chematic_chem::{hba_count, hbd_count, heavy_atom_count, molecular_weight, tpsa};
use serde_json::json;
use std::io::{self, BufRead, Write};

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
            Ok(mol) => json!({
                "index": index,
                "smiles": smiles,
                "status": "ok",
                "descriptors": {
                    "mw": molecular_weight(&mol),
                    "tpsa": tpsa(&mol),
                    "hbd": hbd_count(&mol),
                    "hba": hba_count(&mol),
                    "heavy_atoms": heavy_atom_count(&mol),
                }
            }),
            Err(error) => json!({
                "index": index,
                "smiles": smiles,
                "status": "error",
                "error": error.to_string(),
            }),
        };
        serde_json::to_writer(&mut out, &record).expect("serialize descriptor record");
        writeln!(out).expect("write descriptor record");
    }
}
