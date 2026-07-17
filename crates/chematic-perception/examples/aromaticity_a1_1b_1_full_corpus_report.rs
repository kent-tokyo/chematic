//! Aromaticity-A1-1b-1: run the actual production entry point
//! (`assign_aromaticity_rdkit_parity_experimental`) against a full SMILES
//! corpus and emit the same per-atom/per-bond JSONL schema as
//! `rdkit_parity_full_corpus.rs`, for a direct `diff` against that
//! already-RDKit-verified (100.0000% set agreement) trace, plus explicit
//! error-kind counts for the two `AromaticityError` variants.
//!
//! Byte-identical atom/bond rows here transitively inherit the A1-1b-0
//! full-corpus RDKit-parity result -- this only needs to prove the wiring
//! layer added on top (flag-clearing + internal kekulize + invariant check)
//! doesn't perturb it, not re-join against RDKit again.
//!
//! Run:
//! ```text
//! cargo run -p chematic-perception --release \
//!     --example aromaticity_a1_1b_1_full_corpus_report \
//!     -- ~/Downloads/SMILES.csv \
//!     > /tmp/aromaticity_a1_1b_1_full_corpus_trace.jsonl
//! ```

use std::fs;

use chematic_perception::{
    AromaticityError, apply_aromaticity_rdkit_parity_experimental,
    assign_aromaticity_rdkit_parity_experimental,
};
use serde_json::json;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .expect("usage: aromaticity_a1_1b_1_full_corpus_report <smiles.csv>");
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    let mut n_ok = 0usize;
    let mut n_parse_fail = 0usize;
    let mut n_kekulization_failed = 0usize;
    let mut n_internal_invariant_violation = 0usize;
    let mut n_apply_disagrees_with_assign = 0usize;

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

        let model = match assign_aromaticity_rdkit_parity_experimental(&raw) {
            Ok(m) => m,
            Err(AromaticityError::KekulizationFailed { reason }) => {
                n_kekulization_failed += 1;
                eprintln!("KekulizationFailed: {smi}: {reason}");
                continue;
            }
            Err(AromaticityError::InternalInvariantViolation { reason }) => {
                n_internal_invariant_violation += 1;
                eprintln!("InternalInvariantViolation: {smi}: {reason}");
                continue;
            }
        };

        // Cross-check apply() against assign() on the same input: apply()
        // re-derives its own kekulized molecule internally, so this also
        // catches drift between the two entry points, not just panics.
        if let Ok(applied) = apply_aromaticity_rdkit_parity_experimental(&raw) {
            let disagrees = applied
                .atoms()
                .any(|(idx, atom)| atom.aromatic != model.is_atom_aromatic(idx));
            if disagrees {
                n_apply_disagrees_with_assign += 1;
                eprintln!("APPLY/ASSIGN DISAGREE: {smi}");
            }
        }

        for (idx, _atom) in raw.atoms() {
            let row = json!({
                "kind": "atom",
                "smiles": smi,
                "atom_idx": idx.0,
                "rdkit_parity_atom_aromatic": model.is_atom_aromatic(idx),
            });
            println!("{row}");
        }
        for (bidx, bond) in raw.bonds() {
            let (a1, a2) = (
                bond.atom1.0.min(bond.atom2.0),
                bond.atom1.0.max(bond.atom2.0),
            );
            let row = json!({
                "kind": "bond",
                "smiles": smi,
                "bond_atoms": [a1, a2],
                "rdkit_parity_bond_aromatic": model.is_bond_aromatic(bidx),
            });
            println!("{row}");
        }
        n_ok += 1;
    }

    eprintln!(
        "processed {n_ok} molecules, {n_parse_fail} parse failures, \
         {n_kekulization_failed} KekulizationFailed, \
         {n_internal_invariant_violation} InternalInvariantViolation, \
         {n_apply_disagrees_with_assign} apply/assign disagreements"
    );
}
