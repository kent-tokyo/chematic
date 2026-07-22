//! K1 validation: scan a SMILES corpus and report `chematic_core::kekulize`
//! success/failure per line, one JSON row per molecule, to stdout.
//!
//! Used to produce a before/after diff for the `fix/kekulize-charge-aware-k1`
//! change to `atom_must_be_matched` (crates/chematic-core/src/kekulization.rs)
//! against the project's standard 5,000-molecule SMILES corpus
//! (`~/Downloads/SMILES.csv`, referenced in CLAUDE.md / scripts/bench5k.py).
//! Diagnostic only -- does not modify any production module.
//!
//! Run:
//! ```text
//! cargo run -p chematic-smiles --release --example kekulize_corpus_scan -- \
//!     ~/Downloads/SMILES.csv > /tmp/kekulize_corpus_AFTER.jsonl
//! ```

use std::env;
use std::fs;

/// Minimal JSON string escaping -- avoids pulling in serde_json just for this
/// one-off diagnostic dump (SMILES/error text only ever contains a small,
/// known set of special characters).
fn esc(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn main() {
    let path = env::args()
        .nth(1)
        .expect("usage: kekulize_corpus_scan <smiles_file>");
    let text = fs::read_to_string(&path).expect("failed to read corpus file");

    for (i, line) in text.lines().enumerate() {
        let smiles = line.trim();
        if smiles.is_empty() {
            continue;
        }
        match chematic_smiles::parse(smiles) {
            Ok(mol) => match chematic_core::kekulize(&mol) {
                Ok(_) => println!(
                    r#"{{"idx": {i}, "smiles": "{}", "parse_ok": true, "kekulize_ok": true, "error": null}}"#,
                    esc(smiles)
                ),
                Err(e) => println!(
                    r#"{{"idx": {i}, "smiles": "{}", "parse_ok": true, "kekulize_ok": false, "error": "{}"}}"#,
                    esc(smiles),
                    esc(&e.to_string())
                ),
            },
            Err(e) => println!(
                r#"{{"idx": {i}, "smiles": "{}", "parse_ok": false, "kekulize_ok": null, "error": "{}"}}"#,
                esc(smiles),
                esc(&e.to_string())
            ),
        };
    }
}
