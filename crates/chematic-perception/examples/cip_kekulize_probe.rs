//! Diag/accurate-cip-audit: probes `chematic_core::kekulize()` across a full
//! SMILES corpus, to test (not assume) whether a given CIP residual
//! stereocenter's aromatic environment falls into the kekulize-hard-fail bug
//! already root-caused by the parallel `diag/aromaticity-rdkit-parity`
//! diagnosis (tropylium/imidazolium/pyridinium/pyrylium/tellurophene/
//! phosphole-shaped rings -- `atom_must_be_matched`'s charge-blind
//! lone-pair-donor rules and missing acceptor/Te rules,
//! `crates/chematic-core/src/kekulization.rs`).
//!
//! This is strictly more general than string-matching those 6 fixture SMILES:
//! it re-runs the real algorithm on every corpus molecule, so it also catches
//! *substituted* variants of the same ring classes that the 40-fixture corpus
//! never sampled.
//!
//! Diagnostic only -- does not change `kekulize()`'s behavior. No file under
//! `crates/*/src/**` is touched by this example.
//!
//! Run:
//! ```text
//! cargo run -p chematic-perception --release --example cip_kekulize_probe \
//!     -- ~/Downloads/SMILES.csv \
//!     > validation/results/cip_kekulize_probe.jsonl
//! ```

use std::fs;

use serde_json::json;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: cip_kekulize_probe <smiles.csv>");
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    for line in content.lines() {
        let smi = line.split(',').next().unwrap_or("").trim();
        if smi.is_empty() || smi.eq_ignore_ascii_case("smiles") {
            continue;
        }
        let mol = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(e) => {
                println!(
                    "{}",
                    json!({"smiles": smi, "parse_ok": false, "parse_error": e.to_string()})
                );
                continue;
            }
        };
        match chematic_core::kekulize(&mol) {
            Ok(_) => {
                println!(
                    "{}",
                    json!({"smiles": smi, "parse_ok": true, "kekulize_ok": true})
                );
            }
            Err(e) => {
                println!(
                    "{}",
                    json!({
                        "smiles": smi,
                        "parse_ok": true,
                        "kekulize_ok": false,
                        "kekulize_error": e.to_string(),
                    })
                );
            }
        }
    }
}
