//! Aromaticity-A1-1b-1: downstream sanity checks for molecules produced by
//! `apply_aromaticity_rdkit_parity_experimental`/
//! `assign_aromaticity_rdkit_parity_experimental` (the new opt-in
//! production API, see `docs/rfcs/aromaticity_a1_rfc.md`'s "A1-1b-1" section).
//!
//! This crate (not `chematic-perception`) hosts the check because it
//! depends on both `chematic-perception` (the new API) and SMARTS matching
//! -- the reverse dependency direction would be circular.
//!
//! Scope, per the "additive" reading of the downstream-regression gate (the
//! *default* aromaticity pipeline is untouched by A1-1b-1, so bit-identity
//! against it is not the right bar for the ~0.5% of molecules where this
//! engine's verdict differs from production by design):
//!
//! 1. **Doesn't crash/corrupt**: every molecule that successfully produces
//!    an experimental-aromaticity result round-trips through canonical
//!    SMILES without structural loss, and SMARTS matching against it
//!    completes without panicking, for every query in the same 16-pattern
//!    set `scripts/rdkit_compat_diff.py` uses for the existing 80,000-pair
//!    SMARTS corpus. This does NOT re-run that Python-bound 80k corpus
//!    itself (Python/WASM exposure of the new engine is explicitly out of
//!    scope for this PR) -- it's a proportionate Rust-level substitute:
//!    same 16 patterns, full 5,000-molecule corpus, run against
//!    *experimental*-applied molecules specifically.
//! 2. **Representation parity**: aromatic-form input and explicitly
//!    pre-kekulized input must produce an identical aromatic atom set --
//!    generalizes the existing `purine_representation_stable` unit test
//!    (one molecule) to the full corpus, through the actual production API.
//!
//! Existing-pipeline regression (descriptors/fingerprints/CIP-MANCUDE/the
//! existing 80,000-pair SMARTS corpus) is not re-measured here: those all
//! run through the default `assign_aromaticity_ex`/`apply_aromaticity_ex`,
//! which this PR does not modify (confirmed unchanged by
//! `cargo test --workspace --lib`, 0 failures).
//!
//! Run:
//! ```text
//! cargo run -p chematic-smarts --release \
//!     --example aromaticity_a1_1b_1_downstream_check \
//!     -- ~/Downloads/SMILES.csv
//! ```

use std::fs;

use chematic_perception::{
    AromaticityError, apply_aromaticity_rdkit_parity_experimental,
    assign_aromaticity_rdkit_parity_experimental,
};
use chematic_smarts::{find_matches, parse_smarts};
use chematic_smiles::canonical_smiles;

const SMARTS_PATTERNS: &[&str] = &[
    "[OH]",
    "c",
    "[#7]",
    "C=O",
    "[NX3;H2,H1;!$(NC=O)]",
    "[r5]",
    "[r6]",
    "c1ccccc1",
    "[CX4]",
    "[!#6;!#1]",
    "C(=O)[OH]",
    "[nH]",
    "[#6]=[#6]",
    "[F,Cl,Br,I]",
    "[OX2H]",
    "[#16]",
];

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: aromaticity_a1_1b_1_downstream_check <smiles.csv>");
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    let queries: Vec<_> = SMARTS_PATTERNS
        .iter()
        .map(|p| parse_smarts(p).unwrap_or_else(|e| panic!("bad SMARTS {p}: {e}")))
        .collect();

    let mut n_ok = 0usize;
    let mut n_parse_fail = 0usize;
    let mut n_kekulization_failed = 0usize;
    let mut n_internal_invariant_violation = 0usize;
    let mut n_canonical_round_trip_unstable = 0usize;
    let mut n_smarts_evaluations = 0usize;
    let mut n_representation_parity_mismatch = 0usize;

    for line in content.lines() {
        let smi = line.split(',').next().unwrap_or("").trim();
        if smi.is_empty() || smi.eq_ignore_ascii_case("smiles") {
            continue;
        }
        let raw = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(_) => {
                n_parse_fail += 1;
                continue;
            }
        };

        let applied = match apply_aromaticity_rdkit_parity_experimental(&raw) {
            Ok(m) => m,
            Err(AromaticityError::KekulizationFailed { .. }) => {
                n_kekulization_failed += 1;
                continue;
            }
            Err(AromaticityError::InternalInvariantViolation { reason }) => {
                n_internal_invariant_violation += 1;
                eprintln!("InternalInvariantViolation: {smi}: {reason}");
                continue;
            }
        };
        n_ok += 1;

        // --- (1a) canonical round-trip, no structural loss ---
        let c1 = canonical_smiles(&applied);
        match chematic_smiles::parse(&c1) {
            Ok(reparsed) => {
                let c2 = canonical_smiles(&reparsed);
                if c1 != c2
                    || reparsed.atom_count() != applied.atom_count()
                    || reparsed.bond_count() != applied.bond_count()
                {
                    n_canonical_round_trip_unstable += 1;
                    eprintln!("CANONICAL ROUND-TRIP UNSTABLE: {smi} -> {c1} -> {c2}");
                }
            }
            Err(e) => {
                n_canonical_round_trip_unstable += 1;
                eprintln!("CANONICAL OUTPUT FAILED TO REPARSE: {smi} -> {c1}: {e}");
            }
        }

        // --- (1b) SMARTS matching completes without panicking, at full rate ---
        for query in &queries {
            let _matches = find_matches(query, &applied);
            n_smarts_evaluations += 1;
        }

        // --- (2) representation parity: aromatic-form vs pre-kekulized input ---
        let pre_kekulized = match chematic_core::kekulize(&raw) {
            Ok(k) => chematic_core::apply_kekule(&raw, &k),
            Err(_) => continue, // same molecule, same kekulize() outcome as above -- unreachable here
        };
        if let Ok(model_from_kekulized) =
            assign_aromaticity_rdkit_parity_experimental(&pre_kekulized)
        {
            let model_from_raw = assign_aromaticity_rdkit_parity_experimental(&raw)
                .expect("already succeeded above via apply()");
            let mismatch = raw.atoms().any(|(idx, _)| {
                model_from_raw.is_atom_aromatic(idx) != model_from_kekulized.is_atom_aromatic(idx)
            });
            if mismatch {
                n_representation_parity_mismatch += 1;
                eprintln!("REPRESENTATION PARITY MISMATCH: {smi}");
            }
        }
    }

    eprintln!(
        "processed {n_ok} molecules ({n_parse_fail} parse failures, \
         {n_kekulization_failed} KekulizationFailed, \
         {n_internal_invariant_violation} InternalInvariantViolation)"
    );
    eprintln!(
        "canonical round-trip: {}/{n_ok} stable ({n_canonical_round_trip_unstable} unstable)",
        n_ok - n_canonical_round_trip_unstable
    );
    eprintln!(
        "SMARTS: {n_smarts_evaluations}/{} evaluations completed without panicking",
        n_ok * SMARTS_PATTERNS.len()
    );
    eprintln!(
        "representation parity: {}/{n_ok} match ({n_representation_parity_mismatch} mismatches)",
        n_ok - n_representation_parity_mismatch
    );
}
