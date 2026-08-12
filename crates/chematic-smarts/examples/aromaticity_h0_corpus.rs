//! Aromaticity-H0: build and persist the 94-molecule corpus of canonical
//! round-trip idempotency instability found while validating A1-1b-1
//! (`docs/rfcs/aromaticity_a1_rfc.md`'s "A1-1b-1" section), bucketed by which
//! aromaticity path(s) trigger it:
//!
//! - `experimental_only` (25): stable under `RdkitLike` (current default),
//!   unstable under `RdkitParityExperimental`.
//! - `baseline_only` (7): unstable under `RdkitLike`, stable under
//!   `RdkitParityExperimental`.
//! - `common` (62): unstable under both.
//!
//! Diagnostic only -- reads the full SMILES corpus, writes one frozen JSONL
//! file. No production code touched.
//!
//! Run:
//! ```text
//! cargo run -p chematic-smarts --release --example aromaticity_h0_corpus \
//!     -- ~/Downloads/SMILES.csv \
//!     > validation/aromaticity_h0_corpus.jsonl
//! ```

use std::fs;

use chematic_perception::apply_aromaticity_rdkit_parity_experimental;
use chematic_perception::{AromaticityAlgorithm, apply_aromaticity_ex};
use chematic_smiles::canonical_smiles;
use serde_json::json;

/// `true` if canonicalizing `mol`'s output and re-parsing it does not
/// reproduce the same canonical string (idempotency instability). Returns
/// `None` if the applied molecule can't be canonicalized/re-parsed at all
/// (treated as "not applicable", not folded into stable/unstable).
fn round_trip_unstable(applied: &chematic_core::Molecule) -> Option<bool> {
    let c1 = canonical_smiles(applied);
    let reparsed = chematic_smiles::parse(&c1).ok()?;
    let c2 = canonical_smiles(&reparsed);
    Some(c1 != c2)
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: aromaticity_h0_corpus <smiles.csv>");
    let content =
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    let mut experimental_only = Vec::new();
    let mut baseline_only = Vec::new();
    let mut common = Vec::new();
    let mut n_ok = 0usize;

    for line in content.lines() {
        let smi = line.split(',').next().unwrap_or("").trim();
        if smi.is_empty() || smi.eq_ignore_ascii_case("smiles") {
            continue;
        }
        let raw = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(_) => continue,
        };

        let default_applied = apply_aromaticity_ex(&raw, AromaticityAlgorithm::RdkitLike);
        let Some(default_unstable) = round_trip_unstable(&default_applied) else {
            continue;
        };

        let experimental_applied = match apply_aromaticity_rdkit_parity_experimental(&raw) {
            Ok(m) => m,
            Err(_) => continue, // the one known kekulize-gap molecule; not in this corpus
        };
        let Some(experimental_unstable) = round_trip_unstable(&experimental_applied) else {
            continue;
        };

        n_ok += 1;
        match (default_unstable, experimental_unstable) {
            (false, true) => experimental_only.push(smi.to_string()),
            (true, false) => baseline_only.push(smi.to_string()),
            (true, true) => common.push(smi.to_string()),
            (false, false) => {}
        }
    }

    eprintln!(
        "processed {n_ok} molecules: experimental_only={}, baseline_only={}, common={}, total_unstable_somewhere={}",
        experimental_only.len(),
        baseline_only.len(),
        common.len(),
        experimental_only.len() + baseline_only.len() + common.len()
    );

    for (bucket, smis) in [
        ("experimental_only", &experimental_only),
        ("baseline_only", &baseline_only),
        ("common", &common),
    ] {
        for (i, smi) in smis.iter().enumerate() {
            let row = json!({
                "case_id": format!("{bucket}-{i:03}"),
                "bucket": bucket,
                "smiles": smi,
            });
            println!("{row}");
        }
    }
}
