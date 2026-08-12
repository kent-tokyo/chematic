//! Aromaticity-A1-1b-0: run the RDKit-parity reference engine
//! (`rdkit_parity_aromaticity`) against the diagnosis corpus and emit one
//! JSONL row per (molecule, atom) for RDKit-joining.
//!
//! Diagnostic only -- `rdkit_parity_aromaticity` is not wired into
//! production. See `docs/rfcs/aromaticity_a1_rfc.md`'s "A1-1b-0" section.
//!
//! Run (requires the `diagnostics` feature):
//! ```text
//! cargo run -p chematic-perception --release --features diagnostics \
//!     --example aromaticity_a1_1b_0_report \
//!     -- validation/aromaticity_a1_0_corpus.jsonl \
//!     > validation/results/aromaticity_a1_1b_0_trace.jsonl
//! ```

use std::fs;

use chematic_perception::diagnostics::rdkit_parity_aromaticity;
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
        // Kekulize-then-perceive as step 0, matching RDKit's own pipeline
        // (Kekulize always runs before setAromaticity) -- required
        // precondition of rdkit_parity_aromaticity.
        let mol = match chematic_core::kekulize(&raw) {
            Ok(k) => chematic_core::apply_kekule(&raw, &k),
            Err(e) => {
                eprintln!("SKIP {case_id} ({smiles}): kekulize error: {e}");
                failures += 1;
                continue;
            }
        };
        molecules_seen += 1;

        let (aromatic_atoms, _aromatic_bonds) = rdkit_parity_aromaticity(&mol);

        for (idx, _atom) in mol.atoms() {
            let row = json!({
                "bucket": bucket,
                "case_id": case_id,
                "smiles": smiles,
                "atom_idx": idx.0,
                "rdkit_parity_atom_aromatic": aromatic_atoms.contains(&idx),
            });
            println!("{row}");
            rows_written += 1;
        }
    }

    eprintln!("molecules: {molecules_seen} (failures: {failures}), rows: {rows_written}");
}
