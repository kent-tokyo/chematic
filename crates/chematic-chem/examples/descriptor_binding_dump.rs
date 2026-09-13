//! Emit the core descriptor contract as JSONL for the cross-binding parity lane.
//!
//! The input is one SMILES per line.  This is deliberately a small binding
//! contract (not an external accuracy oracle): Rust, Python, and Node/WASM
//! must agree on parse status and the values produced by the same chematic
//! implementation.

use chematic_chem::{
    fsp3, hbd_count, logp_crippen, molar_refractivity, rdkit_aromatic_ring_count, rdkit_hba_count,
    rdkit_molecular_weight, tpsa,
};
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
                    "molecular_weight": rdkit_molecular_weight(&mol),
                    "tpsa": tpsa(&mol),
                    "hbd": hbd_count(&mol),
                    "hba": rdkit_hba_count(&mol),
                    "logp": logp_crippen(&mol),
                    "molar_refractivity": molar_refractivity(&mol),
                    "fsp3": fsp3(&mol),
                    "aromatic_ring_count": rdkit_aromatic_ring_count(&mol),
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
