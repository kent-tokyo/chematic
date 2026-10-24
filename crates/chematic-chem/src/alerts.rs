//! PAINS (Pan Assay Interference Compounds) structural alert filtering.
//
//! 480 SMARTS patterns from Baell & Holloway 2010 (J. Med. Chem. 53, 2719–2740)
//! as distributed in the RDKit FilterCatalog (pains_a/b/c.in).

use rustc_hash::FxHashMap;
use std::sync::OnceLock;

use chematic_core::{
    Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder, implicit_hcount,
};
use chematic_perception::find_sssr;
mod catalog;

use catalog::{BRENK_SMARTS, PAINS_SMARTS_A, PAINS_SMARTS_B};
use chematic_smarts::{
    MatchConfig, MatchOutcome, QueryMolecule, find_matches_with_rings_and_config_checked,
    has_match_bounded, parse_smarts,
};

/// Hard ceiling on VF2 recursive visits per pattern per molecule.
///
/// Existence-only matching (`has_match_bounded`, below) stops enumerating once
/// a match is *found*, but on a highly symmetric target (e.g. a symmetric
/// macrocycle) VF2 can still backtrack combinatorially just to find — or rule
/// out — that first match. This budget caps that search. Same value already
/// used for the same reason in `standardize::SaltCatalog::is_salt` and
/// `chematic-wasm`'s SMARTS entry points.
///
/// **Exhaustion is not "no match".** When the budget runs out mid-pattern,
/// `has_match_bounded` reports [`MatchOutcome::BudgetExhausted`] — a distinct
/// outcome from `NotFound` — and every function below folds that into the
/// *positive* side of its return value (conservative / fail-safe: an unclear
/// pattern is reported as if it alerted, never silently as clean). Callers
/// that need the real three-way distinction per pattern use the `_checked`
/// variants ([`pains_matches_checked`], [`brenk_matches_checked`]), which
/// `drug_score`'s `f_tox` explicitly does not need to since it just wants a
/// hit count and inherits this same fail-safe fold via `pains_matches`/
/// `brenk_matches`.
const ALERT_MAX_VISIT_BUDGET: u64 = 1_000_000;

/// Return a copy of `mol` with all implicit H atoms materialized as explicit
/// graph nodes.  PAINS patterns (and many RDKit SMARTS) are written for
/// molecules in explicit-H form; this helper enables correct matching.
fn add_explicit_hs(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();

    for (_, atom) in mol.atoms() {
        builder.add_atom(atom.clone());
    }
    for (_, bond) in mol.bonds() {
        builder
            .add_bond(bond.atom1, bond.atom2, bond.order)
            .unwrap();
    }

    let h_elem = Element::from_symbol("H").unwrap();
    let n = mol.atom_count();
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let h_count = implicit_hcount(mol, idx) as usize;
        for _ in 0..h_count {
            let h_idx = builder.add_atom(Atom::new(h_elem));
            builder.add_bond(idx, h_idx, BondOrder::Single).unwrap();
        }
    }

    builder.build()
}

static COMPILED_PAINS: OnceLock<Vec<(&'static str, QueryMolecule)>> = OnceLock::new();
static COMPILED_BRENK: OnceLock<Vec<(&'static str, QueryMolecule)>> = OnceLock::new();

fn compiled_pains_patterns() -> &'static [(&'static str, QueryMolecule)] {
    COMPILED_PAINS.get_or_init(|| {
        PAINS_SMARTS_A
            .iter()
            .chain(PAINS_SMARTS_B.iter())
            .filter_map(|&(name, smarts)| parse_smarts(smarts).ok().map(|q| (name, q)))
            .collect()
    })
}

fn compiled_brenk_patterns() -> &'static [(&'static str, QueryMolecule)] {
    COMPILED_BRENK.get_or_init(|| {
        BRENK_SMARTS
            .iter()
            .filter_map(|&(name, smarts)| parse_smarts(smarts).ok().map(|q| (name, q)))
            .collect()
    })
}

/// Runs `patterns` against `mol` with an existence-only, budget-bounded VF2
/// search, returning one `(name, MatchOutcome)` per pattern — never per
/// embedding (see `MatchOutcome::Found`'s doc: this reports *whether* a rule
/// matched at least once, not how many times or where).
///
/// This is the shared core behind [`pains_matches_checked`] /
/// [`brenk_matches_checked`] and, transitively, every coarser PAINS/Brenk
/// function in this module. Always uses the real production budget
/// ([`ALERT_MAX_VISIT_BUDGET`]) — see [`checked_matches_with_budget`] for the
/// budget-injectable version tests use to exercise the exhaustion path
/// deterministically and near-instantly, without needing a specific
/// pathological molecule or the real (slow) budget.
fn checked_matches(
    mol: &Molecule,
    patterns: &'static [(&'static str, QueryMolecule)],
) -> Vec<(&'static str, MatchOutcome)> {
    checked_matches_with_budget(mol, patterns, ALERT_MAX_VISIT_BUDGET)
}

/// Like [`checked_matches`], but with an injectable visit budget instead of
/// the hardcoded [`ALERT_MAX_VISIT_BUDGET`]. Production call sites always go
/// through [`checked_matches`] (real budget, unchanged behaviour); this
/// exists so tests can force `MatchOutcome::BudgetExhausted` on demand
/// instead of relying on a real pathological molecule actually exhausting
/// the production budget, which can take tens of seconds.
///
/// For patterns that reach the VF2 recursion path, `budget: 0`
/// deterministically produces `BudgetExhausted` on the first attempted
/// visit. Structural early exits in [`has_match_bounded`] — e.g. a query
/// pattern with more atoms than the target molecule — still return
/// `NotFound` without ever consuming budget, regardless of what `budget` is
/// set to.
fn checked_matches_with_budget(
    mol: &Molecule,
    patterns: &'static [(&'static str, QueryMolecule)],
    budget: u64,
) -> Vec<(&'static str, MatchOutcome)> {
    let mol_h = add_explicit_hs(mol);
    let rings = find_sssr(&mol_h);
    let config = MatchConfig {
        max_visit_budget: Some(budget),
        ..Default::default()
    };
    patterns
        .iter()
        .map(|(name, q)| (*name, has_match_bounded(q, &mol_h, &rings, &config)))
        .collect()
}

/// Folds a `checked_matches` result down to `Vec<name>` the same way the
/// coarser PAINS/Brenk functions below do: both `Found` and
/// `BudgetExhausted` count as "alert present" (fail-safe — an unresolved
/// pattern is reported as if it matched, never silently as clean). Only
/// `NotFound` is excluded.
fn fold_conservative(checked: Vec<(&'static str, MatchOutcome)>) -> Vec<&'static str> {
    checked
        .into_iter()
        .filter(|(_, outcome)| *outcome != MatchOutcome::NotFound)
        .map(|(name, _)| name)
        .collect()
}

/// Returns `(name, MatchOutcome)` for every PAINS pattern, distinguishing a
/// confirmed match (`Found`), a confirmed non-match (`NotFound`), and a VF2
/// visit-budget cutoff that resolved neither way (`BudgetExhausted`).
///
/// Use this over [`pains_matches`] when a caller must not conflate "no alert"
/// with "couldn't tell" — [`pains_matches`] necessarily folds the two
/// together into its plain `Vec<&str>` shape (documented there).
pub fn pains_matches_checked(mol: &Molecule) -> Vec<(&'static str, MatchOutcome)> {
    checked_matches(mol, compiled_pains_patterns())
}

/// Returns `(name, MatchOutcome)` for every Brenk pattern — see
/// [`pains_matches_checked`].
pub fn brenk_matches_checked(mol: &Molecule) -> Vec<(&'static str, MatchOutcome)> {
    checked_matches(mol, compiled_brenk_patterns())
}

/// Returns the names of all PAINS patterns that match `mol`, or whose result
/// could not be resolved within the VF2 visit budget (folded in as if
/// matched — see [`pains_matches_checked`] for the real three-way outcome).
///
/// An empty vector means every pattern was conclusively `NotFound`. For
/// *returned* names, this API does not distinguish a confirmed match from
/// `BudgetExhausted` — use [`pains_matches_checked`] for that distinction.
///
/// The molecule is expanded to explicit H form before matching because the
/// PAINS SMARTS patterns reference explicit hydrogen atoms (`[#1]`).
pub fn pains_matches(mol: &Molecule) -> Vec<&'static str> {
    fold_conservative(pains_matches_checked(mol))
}

/// Returns `true` only when every PAINS pattern was confirmed *not* to match.
/// A VF2 visit-budget cutoff on any pattern counts as "does not pass" (see
/// [`pains_matches`]).
pub fn pains_passes(mol: &Molecule) -> bool {
    pains_matches_checked(mol)
        .iter()
        .all(|(_, outcome)| *outcome == MatchOutcome::NotFound)
}

/// Returns the names of all Brenk structural alerts that match `mol`, or
/// whose result could not be resolved within the VF2 visit budget (folded in
/// as if matched — see [`pains_matches`] for the same policy on PAINS).
///
/// An empty vector means every pattern was conclusively `NotFound`. For
/// *returned* names, this API does not distinguish a confirmed match from
/// `BudgetExhausted` — use [`brenk_matches_checked`] for that distinction.
///
/// Brenk et al. 2008 structural alerts for drug design (J. Med. Chem. 51, 5149–5171).
/// These patterns identify problematic functional groups and structural features.
pub fn brenk_matches(mol: &Molecule) -> Vec<&'static str> {
    fold_conservative(brenk_matches_checked(mol))
}

/// Returns `true` only when every Brenk pattern was confirmed *not* to match.
pub fn brenk_passes(mol: &Molecule) -> bool {
    brenk_matches_checked(mol)
        .iter()
        .all(|(_, outcome)| *outcome == MatchOutcome::NotFound)
}

/// Returns `(passes, alert_names)` in a single explicit-H + SSSR + pattern pass.
///
/// Equivalent to calling `pains_passes(mol)` and `pains_matches(mol)` separately,
/// but avoids duplicating the expensive explicit-H conversion and 480-pattern scan.
pub fn pains_passes_and_matches(mol: &Molecule) -> (bool, Vec<&'static str>) {
    let checked = pains_matches_checked(mol);
    let passes = checked
        .iter()
        .all(|(_, outcome)| *outcome == MatchOutcome::NotFound);
    (passes, fold_conservative(checked))
}

/// Returns `(passes, alert_names)` in a single explicit-H + SSSR + pattern pass.
///
/// Equivalent to calling `brenk_passes(mol)` and `brenk_matches(mol)` separately,
/// but avoids duplicating the expensive explicit-H conversion and ~300-pattern scan.
pub fn brenk_passes_and_matches(mol: &Molecule) -> (bool, Vec<&'static str>) {
    let checked = brenk_matches_checked(mol);
    let passes = checked
        .iter()
        .all(|(_, outcome)| *outcome == MatchOutcome::NotFound);
    (passes, fold_conservative(checked))
}

// ---------------------------------------------------------------------------
// Detailed alert matches — returns atom indices for highlighting
// ---------------------------------------------------------------------------

/// Extract unique sorted heavy-atom indices from a set of SMARTS match mappings.
fn heavy_atoms_from_matches(
    matches: &[FxHashMap<usize, AtomIdx>],
    mol_h: &Molecule,
) -> Vec<AtomIdx> {
    use std::collections::BTreeSet;
    let mut set = BTreeSet::new();
    for mapping in matches {
        for &atom_idx in mapping.values() {
            if mol_h.atom(atom_idx).element.atomic_number() != 1 {
                set.insert(atom_idx);
            }
        }
    }
    set.into_iter().collect()
}

/// Shared implementation: scan `patterns` against `mol` and return matching
/// alert names with their heavy-atom indices. Always uses the real
/// production budget ([`ALERT_MAX_VISIT_BUDGET`]) — see
/// [`matches_detailed_impl_with_budget`] for the budget-injectable version
/// tests use to exercise the "found some, then cut off" case deterministically
/// and near-instantly.
///
/// **Budget exhaustion note:** an entry with an empty `atom_indices` vec can
/// mean one of two things — see the public doc comments on
/// [`pains_matches_detailed`]/[`brenk_matches_detailed`] for the caller-facing
/// contract. In particular, when the search is cut off before enumeration
/// completes, any embeddings collected *before* the cutoff are discarded
/// rather than returned as a silently-incomplete highlight set.
fn matches_detailed_impl(
    mol: &Molecule,
    patterns: &'static [(&'static str, QueryMolecule)],
) -> Vec<(&'static str, Vec<AtomIdx>)> {
    matches_detailed_impl_with_budget(mol, patterns, ALERT_MAX_VISIT_BUDGET)
}

/// Like [`matches_detailed_impl`], but with an injectable visit budget. Lets
/// tests deterministically force the "some embeddings were collected, then
/// the budget ran out mid-enumeration" case (a non-zero but insufficient
/// budget) without needing a real pathological molecule and the real (slow)
/// production budget.
fn matches_detailed_impl_with_budget(
    mol: &Molecule,
    patterns: &'static [(&'static str, QueryMolecule)],
    budget: u64,
) -> Vec<(&'static str, Vec<AtomIdx>)> {
    let mol_h = add_explicit_hs(mol);
    let rings = find_sssr(&mol_h);
    // This path enumerates every embedding (needed for highlighting), so it
    // can't be existence-only — but it still gets the same visit-budget
    // safety net as the existence-only variants above, so a pathological
    // symmetric target fails fast instead of hanging.
    let config = MatchConfig {
        max_visit_budget: Some(budget),
        ..Default::default()
    };
    patterns
        .iter()
        .filter_map(|(name, q)| {
            let (ms, budget_exhausted) =
                find_matches_with_rings_and_config_checked(q, &mol_h, &rings, &config);

            if budget_exhausted {
                // Enumeration is incomplete -- whatever embeddings were
                // collected before the cutoff are NOT the full match set, so
                // they must never be returned as though they were exhaustive.
                // Same fail-safe policy as the existence-only functions
                // above: a budget cutoff is reported as present (flagged),
                // never silently dropped like a confirmed non-match.
                return Some((*name, Vec::new()));
            }

            if ms.is_empty() {
                return None;
            }

            Some((*name, heavy_atoms_from_matches(&ms, &mol_h)))
        })
        .collect()
}

/// Like [`pains_matches`] but also returns the matched heavy-atom indices.
///
/// Returns `Vec<(name, atom_indices)>` — one entry per triggered alert. The
/// `atom_indices` identify the heavy atoms in `mol` that are part of the
/// problematic substructure, suitable for use with
/// `chematic_depict::render_svg_highlighted`.
///
/// **An entry with an empty `atom_indices` vec means that alert rule's
/// search was cut off by the VF2 visit budget before enumeration could
/// complete** — not a normal match with zero atoms. Non-empty indices are
/// only ever returned when enumeration ran to completion. Callers must not
/// interpret an empty vec as "matched, but no atoms to highlight"; treat it
/// the same as any other flagged-but-unresolved alert (see
/// [`pains_matches_checked`] for the equivalent distinction on the coarser
/// existence-only API).
pub fn pains_matches_detailed(mol: &Molecule) -> Vec<(&'static str, Vec<AtomIdx>)> {
    matches_detailed_impl(mol, compiled_pains_patterns())
}

/// Like [`brenk_matches`] but also returns the matched heavy-atom indices.
///
/// Returns `Vec<(name, atom_indices)>` — one entry per triggered alert.
///
/// **An entry with an empty `atom_indices` vec means that alert rule's
/// search was cut off by the VF2 visit budget before enumeration could
/// complete** — not a normal match with zero atoms. Non-empty indices are
/// only ever returned when enumeration ran to completion. See
/// [`pains_matches_detailed`] for the full explanation.
pub fn brenk_matches_detailed(mol: &Molecule) -> Vec<(&'static str, Vec<AtomIdx>)> {
    matches_detailed_impl(mol, compiled_brenk_patterns())
}

#[cfg(test)]
#[path = "alerts/tests.rs"]
mod tests;
