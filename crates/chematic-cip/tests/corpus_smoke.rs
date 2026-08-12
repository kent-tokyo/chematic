//! Residual-corpus smoke test (Milestone 1's structural bar, not an accuracy check --
//! see `docs/rfcs/cip_accurate_rfc.md`): every one of the 155 real molecules chematic's
//! existing engine gets wrong vs. RDKit's CIP oracle must still build a digraph without
//! panicking. Milestone 1 has no ranking logic, so this can't check whether the
//! *label* is right -- only that the *structure* is buildable (or, for a pathological
//! input, that expansion fails loudly via `CipError::BudgetExceeded` rather than
//! panicking or silently truncating).

use chematic_cip::{CipBudget, CipDigraph, CipError};
use chematic_core::AtomIdx;

const CORPUS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cip_label_corpus.jsonl"
));

#[test]
fn all_residual_corpus_cases_build_without_panicking() {
    let mut checked = 0usize;
    let mut budget_exceeded = 0usize;

    for line in CORPUS.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let value: serde_json::Value = serde_json::from_str(line).unwrap_or_else(|e| {
            panic!("malformed JSONL line in cip_label_corpus.jsonl: {e}\n{line}")
        });

        // Skip the manifest line (no "smiles" field), don't assume it's always line 1.
        let Some(smiles) = value.get("smiles").and_then(|v| v.as_str()) else {
            continue;
        };
        let atom_idx = value
            .get("atom_idx")
            .and_then(|v| v.as_u64())
            .unwrap_or_else(|| panic!("missing atom_idx for {smiles}"))
            as u32;

        checked += 1;
        let mol = chematic_smiles::parse(smiles)
            .unwrap_or_else(|e| panic!("corpus SMILES failed to parse: {smiles}: {e:?}"));

        let mut digraph = CipDigraph::new(&mol, AtomIdx(atom_idx), CipBudget::default_budget())
            .unwrap_or_else(|e| {
                panic!("digraph root construction failed for {smiles} atom {atom_idx}: {e}")
            });

        match digraph.expand_all(digraph.root()) {
            Ok(()) => {}
            Err(CipError::BudgetExceeded { .. }) => budget_exceeded += 1,
        }
    }

    assert_eq!(
        checked, 155,
        "expected exactly 155 frozen residual cases (corpus manifest count)"
    );
    // A budget failure is an acceptable outcome (see module docs) but should be rare
    // for real molecules under the generous default budget -- flag if it isn't, since
    // that would suggest the budget or the termination rule needs a second look.
    assert!(
        budget_exceeded <= 5,
        "unexpectedly many budget-exceeded cases ({budget_exceeded}/155) -- investigate before trusting this as a smoke test"
    );
}
