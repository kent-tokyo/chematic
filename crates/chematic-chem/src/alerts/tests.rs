use super::*;
use chematic_smiles::parse;

fn mol(s: &str) -> Molecule {
    parse(s).unwrap()
}

#[test]
fn test_pains_caffeine_passes() {
    // Caffeine has no PAINS alerts
    assert!(
        pains_passes(&mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C")),
        "caffeine should pass PAINS"
    );
}

#[test]
fn test_pains_rhodanine_fails() {
    // Rhodanine scaffold triggers ene_rhod_A pattern
    let hits = pains_matches(&mol("O=C1CSC(=S)N1"));
    assert!(!hits.is_empty(), "rhodanine should trigger PAINS alert");
}

#[test]
fn test_pains_azo_fails() {
    // Azo compound: R-N=N-R triggers azo_A
    let hits = pains_matches(&mol("c1ccc(cc1)N=Nc2ccccc2"));
    assert!(!hits.is_empty(), "azobenzene should trigger PAINS azo_A");
}

#[test]
fn test_pains_count_at_most_480() {
    let n = compiled_pains_patterns().len();
    assert!(n <= 480, "at most 480 patterns expected, got {n}");
    assert!(
        n >= 470,
        "at least 470 of 480 patterns should compile, got {n}"
    );
}

#[test]
fn test_brenk_count() {
    let n = compiled_brenk_patterns().len();
    assert!(n > 0, "at least 1 Brenk pattern expected");
    assert!(n <= 105, "at most 105 Brenk patterns expected, got {n}");
}

#[test]
fn test_brenk_phenol_fails() {
    // Phenol is a Brenk alert (phenol pattern)
    let hits = brenk_matches(&mol("Oc1ccccc1"));
    assert!(!hits.is_empty(), "phenol should trigger Brenk alert");
}

#[test]
fn test_brenk_thiol_fails() {
    // Thiol is a Brenk alert
    let hits = brenk_matches(&mol("Cc1ccccc1S"));
    assert!(!hits.is_empty(), "thiol should trigger Brenk alert");
}

#[test]
fn test_aldehyde_formaldehyde_matches() {
    // BUG FIX #1: Formaldehyde should match aldehyde pattern [CH1,CH2]=O
    let hits = brenk_matches(&mol("C=O"));
    assert!(
        hits.contains(&"aldehyde"),
        "formaldehyde should match aldehyde alert"
    );
}

#[test]
fn test_aldehyde_acetaldehyde_matches() {
    // BUG FIX #1 regression: acetaldehyde should still match
    let hits = brenk_matches(&mol("CC=O"));
    assert!(
        hits.contains(&"aldehyde"),
        "acetaldehyde should match aldehyde alert"
    );
}

#[test]
fn test_brenk_aniline_fails() {
    // Aniline is a Brenk alert
    let hits = brenk_matches(&mol("Nc1ccccc1"));
    assert!(!hits.is_empty(), "aniline should trigger Brenk alert");
}

#[test]
fn test_pains_catechol_fails() {
    // Catechol (ortho-dihydroxybenzene) → catechol_A pattern
    let hits = pains_matches(&mol("Oc1ccccc1O"));
    assert!(!hits.is_empty(), "catechol should trigger PAINS catechol_A");
}

#[test]
fn test_pains_mannich_base_fails() {
    // Mannich base: N-CH2-ArOH pattern → mannich_A
    let hits = pains_matches(&mol("NCc1ccccc1O"));
    assert!(
        !hits.is_empty(),
        "Mannich base should trigger PAINS mannich_A"
    );
}

/// Pins the bounded negative result for a formerly expensive symmetric
/// scaffold so the VF2 ordering improvement cannot silently regress.
///
/// `CC(C)(C)c1cc(C(=O)/C=C/c2ccsc2)cc(C(C)(C)C)c1O` is a 2,6-di-tert-
/// -butylphenol (BHT-class antioxidant) with a thienyl chalcone
/// substituent -- an ordinary drug-like structure, not an exotic
/// macrocycle. Its two tert-butyl groups give the `tert_butyl_B(1)`
/// PAINS pattern's search enough local symmetry (each tert-butyl's 3
/// methyls, each methyl's 3 H's, both tert-butyl groups mutually
/// symmetric under the ring) to blow up the old lowest-index VF2 order.
///
/// The correct answer is `NotFound` (`tert_butyl_B` genuinely does not
/// match). The constrained-query VF2 ordering now proves that negative
/// within the production budget, removing the former fail-safe false
/// positive without changing the conservative treatment of genuine
/// budget exhaustion.
///
/// This test intentionally exercises the real production budget. It is a
/// bounded regression for #139, so it runs in the default suite and keeps
/// the measured path visible to CI.
#[test]
fn test_pains_di_tert_butylphenol_resolves_negative_within_budget() {
    let smiles = "CC(C)(C)c1cc(C(=O)/C=C/c2ccsc2)cc(C(C)(C)C)c1O";
    let start = std::time::Instant::now();
    let checked = pains_matches_checked(&mol(smiles));
    let elapsed = start.elapsed();

    let outcome = checked
        .iter()
        .find(|(name, _)| *name == "tert_butyl_B(1)")
        .map(|(_, outcome)| *outcome);
    assert_eq!(
        outcome,
        Some(MatchOutcome::NotFound),
        "the constrained VF2 order should prove this known negative \
         within the production budget"
    );

    // A proven negative must not be surfaced as a conservative alert.
    assert!(
        !fold_conservative(checked.clone()).contains(&"tert_butyl_B(1)"),
        "a proven NotFound must not be surfaced as a PAINS alert"
    );

    assert!(
        elapsed < std::time::Duration::from_secs(10),
        "expected the constrained search to resolve within 10s -- took \
         {elapsed:?}"
    );
}

/// Fast, deterministic default-suite replacement for the `#[ignore]`d
/// real-budget calibration test above. Reuses the same real molecule and
/// pattern (`tert_butyl_B(1)` on the di-tert-butylphenol scaffold -- big
/// enough that this pattern doesn't hit `has_match_bounded`'s "query
/// bigger than target" early exit, the same way the slow test's `Found`
/// is never returned in that case either) but with `budget: 0` instead of
/// the real production budget: the very first VF2 visit is already
/// abandoned (see `match_vf2::match_recursive`'s early-exit branch), so
/// the outcome is deterministically `BudgetExhausted` and near-instant --
/// no backtracking search actually runs. This exercises the thing that
/// matters -- budget exhaustion folds to a conservative alert, never
/// silently to "no match" -- without waiting on the real (slow) budget.
#[test]
fn checked_matches_zero_budget_is_deterministic_budget_exhausted() {
    let smiles = "CC(C)(C)c1cc(C(=O)/C=C/c2ccsc2)cc(C(C)(C)C)c1O";
    let checked = checked_matches_with_budget(&mol(smiles), compiled_pains_patterns(), 0);

    let outcome = checked
        .iter()
        .find(|(name, _)| *name == "tert_butyl_B(1)")
        .map(|(_, outcome)| *outcome);
    assert_eq!(
        outcome,
        Some(MatchOutcome::BudgetExhausted),
        "budget=0 must resolve tert_butyl_B(1) to BudgetExhausted -- a \
         zero budget can never prove a negative (or a positive)"
    );
    assert!(
        fold_conservative(checked).contains(&"tert_butyl_B(1)"),
        "BudgetExhausted must fold into the flagged/alert side, never \
         silently dropped as if NotFound"
    );
}

/// Regression test for a real bug: `matches_detailed_impl` used to check
/// `ms.is_empty()` before `budget_exhausted`, so when the search found
/// *some* embeddings before the budget ran out, it returned those
/// partial embeddings as though they were the complete, exhaustive match
/// set -- silently misrepresenting an incomplete highlight set as a full
/// one.
///
/// `budget: 3` against benzene's 6-fold-symmetric `c` query
/// deterministically finds 2 (of 6 possible) embeddings before
/// exhausting -- empirically verified: budget 0-1 -> 0 matches
/// (exhausted), budget 2 -> 1 (exhausted), budget 3 -> 2 (exhausted), ...,
/// budget 7+ -> all 6 (not exhausted). Budget 3 is exactly the "found
/// some, not all" case the fix targets, deterministically and
/// near-instantly (benzene is tiny; no real backtracking search needed).
#[test]
fn matches_detailed_never_returns_partial_highlights_on_budget_exhaustion() {
    let benzene = mol("c1ccccc1");
    let query = parse_smarts("c").unwrap();
    let patterns: &'static [(&'static str, QueryMolecule)] =
        Box::leak(vec![("test_c", query)].into_boxed_slice());

    let result = matches_detailed_impl_with_budget(&benzene, patterns, 3);

    assert_eq!(
        result.len(),
        1,
        "the pattern must still be reported as flagged, not dropped, \
         when its search was cut off"
    );
    assert_eq!(result[0].0, "test_c");
    assert!(
        result[0].1.is_empty(),
        "a budget cutoff mid-enumeration must report an EMPTY highlight \
         set, never the partial embeddings collected before the cutoff \
         -- those are not the complete, exhaustive match"
    );
}
