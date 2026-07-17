//! Aromaticity-H0: classify every molecule in the 94-case canonical
//! round-trip idempotency corpus (`validation/aromaticity_h0_corpus.jsonl`,
//! built pre-Kekule-S0) into its root cause, post-Kekule-S0.
//!
//! Diagnostic only -- reads frozen corpus input, writes one classification
//! report. No production code touched.
//!
//! For each of the 94 molecules, re-checks canonical round-trip stability
//! under both `RdkitLike` (current default) and
//! `apply_aromaticity_rdkit_parity_experimental` (A1-1b-1's opt-in engine),
//! post-Kekule-S0, and classifies:
//!
//! - `fixed_by_kekule_s0`: now stable under both paths (was NOT stable
//!   under at least one, pre-fix) -- genuinely resolved by
//!   `apply_kekule`'s `stereo_neighbor_order` preservation.
//! - `bond_level_aromaticity_disagreement`: still unstable under at least
//!   one path; atom-level `.aromatic` flags AND `stereo_neighbor_order`
//!   are identical between the two paths, but at least one BOND's `order`
//!   differs (`RdkitLike` promotes a bond to `Aromatic` that the
//!   RDKit-parity engine leaves `Single`, or vice versa) -- a genuine
//!   aromaticity-model disagreement between the two engines (part of the
//!   already-known ~1.18% bond-level gap between current production and
//!   RDKit), which then triggers the *same* pre-existing canonicalization
//!   sensitivity as the next category.
//! - `canonicalizer_representation_sensitivity`: still unstable under at
//!   least one path; atom flags, stereo_neighbor_order, AND all bond
//!   orders are identical between the two paths -- the instability exists
//!   independent of which aromaticity engine is used, i.e. a pre-existing
//!   canonical-SMILES writer sensitivity on this specific
//!   fused/complex-ring molecule, not an aromaticity-assignment or
//!   stereo-metadata issue at all.
//! - `unexplained`: none of the above -- would need further tracing (goal
//!   is 0 of these).
//!
//! Run:
//! ```text
//! cargo run -p chematic-smarts --release --example aromaticity_h0_classify \
//!     -- validation/aromaticity_h0_corpus.jsonl \
//!     > validation/results/aromaticity_h0_classification.jsonl
//! ```

use std::fs;

use chematic_perception::{
    AromaticityAlgorithm, apply_aromaticity_ex, apply_aromaticity_rdkit_parity_experimental,
};
use chematic_smiles::canonical_smiles;
use serde_json::{Value, json};

fn round_trip_stable(applied: &chematic_core::Molecule) -> bool {
    let c1 = canonical_smiles(applied);
    match chematic_smiles::parse(&c1) {
        Ok(reparsed) => c1 == canonical_smiles(&reparsed),
        Err(_) => false,
    }
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "validation/aromaticity_h0_corpus.jsonl".to_string());
    let content =
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    let mut counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();

    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let case: Value = serde_json::from_str(line).expect("valid corpus JSON line");
        let smi = case["smiles"].as_str().unwrap().to_string();
        let case_id = case["case_id"].as_str().unwrap_or("").to_string();
        let bucket_before = case["bucket"].as_str().unwrap_or("").to_string();

        let raw = chematic_smiles::parse(&smi).expect("valid SMILES (already in frozen corpus)");
        let default_applied = apply_aromaticity_ex(&raw, AromaticityAlgorithm::RdkitLike);
        let Ok(experimental_applied) = apply_aromaticity_rdkit_parity_experimental(&raw) else {
            eprintln!("SKIP {case_id} ({smi}): KekulizationFailed post-fix (unexpected)");
            continue;
        };

        let default_stable = round_trip_stable(&default_applied);
        let experimental_stable = round_trip_stable(&experimental_applied);

        let classification = if default_stable && experimental_stable {
            "fixed_by_kekule_s0"
        } else {
            let atom_flags_match = raw.atoms().all(|(idx, _)| {
                default_applied.atom(idx).aromatic == experimental_applied.atom(idx).aromatic
            });
            let stereo_match = raw.atoms().all(|(idx, _)| {
                default_applied.stereo_neighbor_order(idx)
                    == experimental_applied.stereo_neighbor_order(idx)
            });
            let bond_orders_match = raw.bonds().all(|(bidx, _)| {
                default_applied.bond(bidx).order == experimental_applied.bond(bidx).order
            });

            if atom_flags_match && stereo_match && !bond_orders_match {
                "bond_level_aromaticity_disagreement"
            } else if atom_flags_match && stereo_match && bond_orders_match {
                "canonicalizer_representation_sensitivity"
            } else {
                "unexplained"
            }
        };

        *counts.entry(classification).or_insert(0) += 1;

        let row = json!({
            "case_id": case_id,
            "smiles": smi,
            "bucket_before_kekule_s0": bucket_before,
            "default_stable_now": default_stable,
            "experimental_stable_now": experimental_stable,
            "classification": classification,
        });
        println!("{row}");
    }

    eprintln!("\n=== classification counts ===");
    let mut sorted: Vec<_> = counts.into_iter().collect();
    sorted.sort_by_key(|(_, n)| std::cmp::Reverse(*n));
    let mut total = 0;
    for (k, n) in &sorted {
        eprintln!("{k}: {n}");
        total += n;
    }
    eprintln!("total classified: {total}");
}
