//! Morgan/ECFP RDKit environment-parity diagnostic PR: dumps chematic's
//! real per-`(atom, radius)` Morgan trace (see
//! `crates/chematic-fp/src/ecfp_diagnostics.rs`) for a SMILES corpus, one
//! JSON object per molecule (not per atom, so the Python comparator can
//! join by SMILES key).
//!
//! Uses `EcfpInvariantMode::RdkitMorgan`, radius 2, 2048 bits, chirality
//! off -- matching the RDKit oracle's pinned Morgan options (see
//! `scripts/gen_ecfp_rdkit_environment_oracle.py`).
//!
//! A SMILES chematic fails to parse gets a `parse_ok: false` row instead of
//! being silently dropped.
//!
//! Diagnostic only. Requires the `diagnostics` feature:
//! ```text
//! cargo run -p chematic-fp --release --features diagnostics \
//!     --example morgan_rdkit_environment_trace -- <SMILES.csv> <out.jsonl>
//! ```

use chematic_fp::diagnostics::morgan_trace;
use chematic_fp::{EcfpConfig, EcfpInvariantMode};
use chematic_smiles::parse;
use serde_json::json;
use std::fs;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args.get(1).cloned().unwrap_or_else(|| {
        panic!("usage: morgan_rdkit_environment_trace <SMILES.csv> <out.jsonl>")
    });
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "trace.jsonl".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));
    let config = EcfpConfig::default(); // radius 2, nbits 2048, chirality off

    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    let mut input_lines = 0usize;
    let mut parse_fail = 0usize;
    let mut rows = 0usize;

    for line in content.lines() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        input_lines += 1;
        // 0-based position in this file's row stream -- must match the RDKit
        // oracle's row_id exactly, since both are built from the same
        // ordered SMILES list (see the comparator's module docstring). This
        // is the explicit sync check the position-only design relied on
        // implicitly before.
        let row_id = rows;

        let mol = match parse(smi) {
            Ok(m) => m,
            Err(e) => {
                parse_fail += 1;
                writeln!(
                    f,
                    "{}",
                    json!({"row_id": row_id, "smiles": smi, "parse_ok": false, "error": e.to_string()})
                )
                .unwrap();
                rows += 1;
                continue;
            }
        };

        let atomic_numbers: Vec<u8> = (0..mol.atom_count())
            .map(|i| {
                mol.atom(chematic_core::AtomIdx(i as u32))
                    .element
                    .atomic_number()
            })
            .collect();

        let trace = morgan_trace(&mol, &config, EcfpInvariantMode::RdkitMorgan);
        let trace_json: Vec<_> = trace
            .iter()
            .map(|e| {
                json!({
                    "atom_idx": e.atom_idx,
                    "radius": e.radius,
                    "raw_environment_id": e.raw_environment_id,
                    "folded_bit": e.folded_bit,
                    "emitted": e.emitted,
                    "atom_ball": e.atom_ball,
                })
            })
            .collect();

        writeln!(
            f,
            "{}",
            json!({
                "row_id": row_id,
                "smiles": smi,
                "parse_ok": true,
                "atom_count": mol.atom_count(),
                "atomic_numbers": atomic_numbers,
                "trace": trace_json,
            })
        )
        .unwrap();
        rows += 1;
    }

    eprintln!("input_lines={input_lines} parse_fail={parse_fail} rows={rows} out={out_path}");
}
