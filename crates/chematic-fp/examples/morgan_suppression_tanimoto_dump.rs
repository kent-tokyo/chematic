//! Phase B Tanimoto-correlation-vs-RDKit before/after record (non-gating
//! reference measurement, same sampling frame as PR #120's
//! `ecfp_rdkit_environment_parity.py::tanimoto_correlation`: population from
//! this dump's rows, sample_size=300, seed=42, computed in
//! `scripts/ecfp_rdkit_suppression_tanimoto.py`).
//!
//! Dumps, per molecule, the folded 2048-bit on-bit list from BOTH:
//!   - "baseline": `ecfp4_rdkit_invariants()` (existing, pre-Phase-B path)
//!   - "suppression": `ecfp4_rdkit_environment_experimental()` (this PR's new path)
//!
//! for correlation against RDKit's own real `default.folded_on_bits` from
//! the same oracle rows used throughout this PR's validation.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example morgan_suppression_tanimoto_dump \
//!     -- <SMILES.csv> <out.jsonl>
//! ```

use chematic_fp::{BitVec2048, ecfp4_rdkit_environment_experimental, ecfp4_rdkit_invariants};
use chematic_smiles::parse;
use serde_json::json;
use std::fs;
use std::io::Write;

fn on_bits(fp: &BitVec2048) -> Vec<usize> {
    (0..2048).filter(|&i| fp.get(i)).collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args.get(1).cloned().unwrap_or_else(|| {
        panic!("usage: morgan_suppression_tanimoto_dump <SMILES.csv> <out.jsonl>")
    });
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "tanimoto.jsonl".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));
    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    let mut rows = 0usize;

    for line in content.lines() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        let row_id = rows;
        match parse(smi) {
            Ok(mol) => {
                let baseline_bits = on_bits(&ecfp4_rdkit_invariants(&mol));
                let suppression_bits = on_bits(&ecfp4_rdkit_environment_experimental(&mol));
                writeln!(
                    f,
                    "{}",
                    json!({
                        "row_id": row_id,
                        "smiles": smi,
                        "parse_ok": true,
                        "baseline_folded_bits": baseline_bits,
                        "suppression_folded_bits": suppression_bits,
                    })
                )
                .unwrap();
            }
            Err(e) => {
                writeln!(
                    f,
                    "{}",
                    json!({"row_id": row_id, "smiles": smi, "parse_ok": false, "error": e.to_string()})
                )
                .unwrap();
            }
        }
        rows += 1;
    }

    eprintln!("rows={rows} out={out_path}");
}
