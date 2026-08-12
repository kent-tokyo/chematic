//! DrugScore — composite drug-likeness score.
//!
//! Product of four independent factors (each in 0–1):
//! - cLogP: mountain-shaped sigmoid, optimal 0–5
//! - logS: solubility, optimal > −5 log mol/L
//! - MW: molecular weight, optimal < 500 Da
//! - Toxicity: penalty per PAINS (×0.5) and Brenk (×0.75) alert

use crate::alerts::{brenk_matches, pains_matches};
use crate::descriptors::{logp_crippen, molecular_weight};
use crate::esol::esol_solubility;
use chematic_core::Molecule;

// ---------------------------------------------------------------------------
// Individual factor functions (0–1)
// ---------------------------------------------------------------------------

fn sigmoid(x: f64, k: f64) -> f64 {
    1.0 / (1.0 + (-k * x).exp())
}

// Mountain-shaped: rises at LogP = 0, falls at LogP = 5
fn f_logp(logp: f64) -> f64 {
    sigmoid(logp, 2.0) * (1.0 - sigmoid(logp - 5.0, 2.0))
}

// Right-ascending: logS > −5 is optimal
fn f_logs(logs: f64) -> f64 {
    sigmoid(logs + 5.0, 1.0)
}

// Right-descending: MW < 500 is optimal
fn f_mw(mw: f64) -> f64 {
    1.0 - sigmoid(mw - 500.0, 0.02)
}

// Toxicity penalty: each PAINS halves the score, each Brenk reduces by 25%
fn f_tox(mol: &Molecule) -> f64 {
    let n_pains = pains_matches(mol).len() as i32;
    let n_brenk = brenk_matches(mol).len() as i32;
    let score = 0.5f64.powi(n_pains) * 0.75f64.powi(n_brenk);
    score.max(0.05) // floor at 5% to avoid complete suppression
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Composite drug-likeness score (0–1, higher = more drug-like).
///
/// Combines four independent properties as a product of sigmoidal factors:
/// - **cLogP** (Crippen): optimal 0–5, penalty outside
/// - **logS** (ESOL): solubility optimal > −5 log mol/L
/// - **MW**: molecular weight optimal < 500 Da
/// - **Toxicity**: 0.5× per PAINS alert, 0.75× per Brenk alert
///
/// Comparable to OCL DrugScore (LogP × LogS × MW × toxicity product).
/// Use alongside [`qed`](crate::qed) for a second opinion on drug-likeness.
pub fn drug_score(mol: &Molecule) -> f64 {
    f_logp(logp_crippen(mol))
        * f_logs(esol_solubility(mol))
        * f_mw(molecular_weight(mol))
        * f_tox(mol)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> chematic_core::Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn drug_score_in_range_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
        let score = drug_score(&m);
        assert!(
            score > 0.0 && score <= 1.0,
            "drug_score must be in (0, 1]: {score}"
        );
    }

    #[test]
    fn drug_score_clean_molecule_scores_well() {
        // Ibuprofen: moderate LogP (~3.5), no PAINS, MW=206 — clean profile
        let m = mol("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O");
        let score = drug_score(&m);
        assert!(score > 0.4, "ibuprofen should score > 0.4, got {score:.3}");
    }

    #[test]
    fn drug_score_large_molecule_penalised() {
        // Very large molecule (MW >> 500) should score low due to MW penalty
        let heavy = mol("CCCCCCCCCCCCCCCCCCCCC"); // henicosane, MW ~296, LogP ~11
        let small = mol("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"); // ibuprofen
        let s_heavy = drug_score(&heavy);
        let s_small = drug_score(&small);
        assert!(s_heavy.is_finite() && s_small.is_finite());
        // High LogP penalises heavily
        assert!(
            s_heavy < s_small,
            "high-LogP molecule should score lower: {s_heavy:.3} vs {s_small:.3}"
        );
    }

    #[test]
    #[ignore = "slow real-budget VF2 calibration (~200ms release / ~3-5s debug \
                in isolation, but a real-budget VF2 search is a poor fit for \
                the default suite under concurrent CI load) -- run explicitly \
                with `cargo test -p chematic-chem --lib -- --ignored \
                drug_score_pathological_macrocycle_does_not_hang`. The \
                underlying budget-exhaustion-folds-to-flagged mechanism this \
                test exercises has a fast, deterministic default-suite test \
                in alerts.rs \
                (checked_matches_zero_budget_is_deterministic_budget_exhausted)."]
    fn drug_score_pathological_macrocycle_does_not_hang() {
        // Regression test for the exact repro molecule from PR #137's descriptor
        // census RFC (`docs/rfcs/descriptor_census_rfc.md` on `diag/descriptor-census`):
        // a symmetric bis-isoquinolinium macrocycle that, before this fix, drove
        // `drug_score` -> `alerts::pains_matches` -> `match_vf2::match_recursive`
        // into a multi-minute VF2 combinatorial blowup (confirmed via a stack
        // sampler showing zero forward progress ~15 recursion levels deep).
        //
        // The fix adds a hard VF2 visit-budget cap (`MatchConfig::max_visit_budget`)
        // to every PAINS/Brenk match site in `alerts.rs`, on top of the existing
        // existence-only `max_matches: Some(1)` short-circuit — so a pathological
        // symmetric target fails fast instead of hanging. A budget cutoff is
        // never silently read as "no match": it resolves to
        // `MatchOutcome::BudgetExhausted`, which folds into `pains_matches`'s
        // output as if it *had* matched (fail-safe: flag when uncertain, never
        // silently clean — see `alerts::checked_matches`/`fold_conservative`).
        // This test has no per-test timeout mechanism to hook into (none exists
        // elsewhere in this suite), so it asserts wall-clock elapsed time
        // directly; 30s is generous slack over the ~200ms (release) / ~5s (debug)
        // observed locally, while still catching a regression back to "hangs".
        let m =
            mol("C1=C\\c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCNc2cc[n+](c3ccccc23)Cc2ccc/1cc2");
        let start = std::time::Instant::now();
        let score = drug_score(&m);
        let elapsed = start.elapsed();
        assert!(score.is_finite(), "drug_score must return a finite value");
        assert!(
            elapsed < std::time::Duration::from_secs(30),
            "drug_score on the pathological macrocycle took {elapsed:?}, expected it to \
             fail fast within the VF2 visit budget"
        );
    }

    #[test]
    fn drug_score_factors_independent() {
        // Each factor function stays in [0, 1]
        for logp in [-5.0, 0.0, 2.5, 5.0, 10.0] {
            let f = f_logp(logp);
            assert!((0.0..=1.0).contains(&f), "f_logp({logp}) = {f}");
        }
        for logs in [-10.0, -5.0, -2.0, 0.0] {
            let f = f_logs(logs);
            assert!((0.0..=1.0).contains(&f), "f_logs({logs}) = {f}");
        }
        for mw in [100.0, 300.0, 500.0, 800.0] {
            let f = f_mw(mw);
            assert!((0.0..=1.0).contains(&f), "f_mw({mw}) = {f}");
        }
    }
}
