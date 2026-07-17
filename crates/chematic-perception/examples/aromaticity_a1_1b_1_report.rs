//! Aromaticity-A1-1b-1: run the actual production entry point
//! (`assign_aromaticity_rdkit_parity_experimental`) against the diagnosis
//! corpus and emit the same per-(molecule, atom) JSONL schema as
//! `aromaticity_a1_1b_0_report.rs`, for a direct `diff` against the
//! already-RDKit-verified `aromaticity_a1_1b_0_trace.jsonl`.
//!
//! Unlike the A1-1b-0 report, this passes the *raw* (possibly
//! aromatic-form) parsed molecule straight to the production function --
//! no manual pre-kekulization -- so it actually exercises the wiring layer
//! (`clear_aromatic_flags` + internal `kekulize`/`apply_kekule` +
//! `validate_aromaticity_invariants`) that A1-1b-1 added on top of the
//! already-verified `rdkit_parity_aromaticity` engine. A byte-identical
//! diff against the A1-1b-0 trace is the actual proof that this wiring is
//! mechanical (in particular, that clearing stale aromatic flags before
//! kekulizing does not perturb aromatic-form heteroaromatic input).
//!
//! Run:
//! ```text
//! cargo run -p chematic-perception --release \
//!     --example aromaticity_a1_1b_1_report \
//!     -- validation/aromaticity_a1_0_corpus.jsonl \
//!     > /tmp/aromaticity_a1_1b_1_trace.jsonl
//! diff /tmp/aromaticity_a1_1b_1_trace.jsonl validation/results/aromaticity_a1_1b_0_trace.jsonl
//! ```

use std::fs;

use chematic_perception::assign_aromaticity_rdkit_parity_experimental;
use serde_json::{Value, json};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let corpus_path = args
        .get(1)
        .map(String::as_str)
        .unwrap_or("validation/aromaticity_a1_0_corpus.jsonl");

    let corpus = fs::read_to_string(corpus_path)
        .unwrap_or_else(|e| panic!("failed to read {corpus_path}: {e}"));

    let mut rows_written = 0usize;
    let mut molecules_seen = 0usize;
    let mut failures = 0usize;

    for line in corpus.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let case: Value = serde_json::from_str(line).expect("valid corpus JSON line");
        let bucket = case["bucket"].as_str().unwrap_or("unknown").to_string();
        let case_id = case["case_id"].as_str().unwrap_or("").to_string();
        let smiles = case["smiles"].as_str().unwrap_or("").to_string();

        let raw = match chematic_smiles::parse(&smiles) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("SKIP {case_id} ({smiles}): parse error: {e}");
                failures += 1;
                continue;
            }
        };

        // The production entry point does its own kekulization internally
        // -- pass the raw parsed (aromatic-form) molecule directly, not a
        // manually pre-kekulized one.
        let model = match assign_aromaticity_rdkit_parity_experimental(&raw) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("SKIP {case_id} ({smiles}): {e}");
                failures += 1;
                continue;
            }
        };
        molecules_seen += 1;

        for (idx, _atom) in raw.atoms() {
            let row = json!({
                "bucket": bucket,
                "case_id": case_id,
                "smiles": smiles,
                "atom_idx": idx.0,
                "rdkit_parity_atom_aromatic": model.is_atom_aromatic(idx),
            });
            println!("{row}");
            rows_written += 1;
        }
    }

    eprintln!("molecules: {molecules_seen} (failures: {failures}), rows: {rows_written}");
}
