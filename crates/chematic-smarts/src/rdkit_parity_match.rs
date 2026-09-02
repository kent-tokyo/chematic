//! Opt-in RDKit-parity SMARTS matching mode.
//!
//! This module is a **deliberate, near-total duplication** of
//! [`crate::match_vf2`]'s VF2 recursive matcher, not a refactor of it. The
//! duplication is intentional: it keeps `match_vf2.rs` at **zero diff**, so
//! "the default matcher is byte-identical before and after this change" is
//! trivially provable (there is no change to prove anything about) rather
//! than something that has to be argued from a shared-code refactor. The
//! only place this module's evaluator actually diverges from
//! `match_vf2.rs`'s is `AtomPrimitive::RingCount` (`[RN]`), which consults
//! [`crate::rdkit_ring_model`] instead of a plain SSSR count — see that
//! module's doc comment for the full root-cause/design rationale (RDKit's
//! `symmetrizeSSSR`) and for why every *other* ring-shaped primitive
//! (`[R]`/`[R0]`, `[rN]`, `[kN]`, `[xN]`, ring-bond `@`/`!@`) is left
//! wired to the identical plain-SSSR formula `match_vf2.rs` already uses.
//!
//! Everything else in this file — chirality/isotope handling, valence,
//! hybridization, recursive-SMARTS anchoring, the visit-budget/
//! `MatchOutcome` contract — is copied verbatim in behavior from
//! `match_vf2.rs` so this mode's non-ring-count results are identical to
//! the default matcher's.

use rustc_hash::{FxHashMap, FxHashSet};

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};
use chematic_perception::RingSet;

use crate::match_vf2::{MatchConfig, MatchOutcome};
use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};
use crate::rdkit_ring_model::{
    RdkitParityError, RdkitParityRingModel, RdkitRingModelBudget, build_rdkit_parity_ring_model,
};

/// Configuration for [`find_matches_rdkit_parity`] / [`has_match_rdkit_parity_bounded`].
#[derive(Debug, Clone, Default)]
pub struct RdkitParityConfig {
    /// Everything [`MatchConfig`] already covers (chirality/isotope
    /// enforcement, `max_matches`, `uniquify`, VF2 `max_visit_budget`) —
    /// unaffected by, and orthogonal to, this mode's ring-count model.
    pub base: MatchConfig,
    /// Resource bound on the RDKit-parity ring-count model's candidate
    /// search. See [`RdkitRingModelBudget`].
    pub ring_model_budget: RdkitRingModelBudget,
    /// When `true`, the target molecule is re-perceived with
    /// `chematic_perception::apply_aromaticity_rdkit_parity_experimental`
    /// before matching (addresses the bridgehead-N ring-fusion
    /// over-aromatization residual named "SMARTS-A0" in
    /// `docs/rdkit_compat.md`). **Default `false`** and deliberately kept
    /// as a separate flag from the ring-count model above: conflating the
    /// two would make any mismatch unattributable to either mechanism. On
    /// failure (a known limitation of that engine on a small class of
    /// bridgehead-heteroatom-fused rings), the whole call returns
    /// `Err(RdkitParityError::Aromaticity(..))` — never a silent fallback
    /// to the default Hückel engine's flags.
    pub use_rdkit_parity_aromaticity: bool,
}

/// Find all non-overlapping (injective) embeddings of `query` in `mol` using
/// the opt-in RDKit-parity ring-count model for `[RN]`.
///
/// Returns `(matches, budget_exhausted)` exactly like
/// [`crate::find_matches_with_rings_and_config_checked`] — `budget_exhausted`
/// is the VF2 state-space search budget (`config.base.max_visit_budget`),
/// never conflated with the separate ring-model budget below (that one
/// surfaces as `Err`, not as a flag on a successful result).
pub fn find_matches_rdkit_parity(
    query: &QueryMolecule,
    mol: &Molecule,
    config: &RdkitParityConfig,
) -> Result<(Vec<FxHashMap<usize, AtomIdx>>, bool), RdkitParityError> {
    if query.atoms.is_empty() {
        return Ok((vec![], false));
    }
    if query.atoms.len() > mol.atom_count() {
        return Ok((vec![], false));
    }

    let mol_owned;
    let mol_ref: &Molecule = if config.use_rdkit_parity_aromaticity {
        mol_owned = chematic_perception::apply_aromaticity_rdkit_parity_experimental(mol)
            .map_err(RdkitParityError::Aromaticity)?;
        &mol_owned
    } else {
        mol
    };

    let rings = chematic_perception::find_sssr(mol_ref);
    let ring_model = if query_uses_ring_count(query) {
        Some(build_rdkit_parity_ring_model(
            mol_ref,
            &rings,
            &config.ring_model_budget,
        )?)
    } else {
        None
    };

    let ctx = EvalCtx {
        mol: mol_ref,
        rings: &rings,
        ring_model: ring_model.as_ref(),
        config: &config.base,
        visit_budget: std::cell::Cell::new(config.base.max_visit_budget.unwrap_or(u64::MAX)),
        budget_exhausted: std::cell::Cell::new(false),
        min_ring_size_by_atom: std::cell::RefCell::new(None),
    };
    let mut mapping: FxHashMap<usize, AtomIdx> = FxHashMap::default();
    let mut results: Vec<FxHashMap<usize, AtomIdx>> = Vec::new();
    match_recursive(
        query,
        &ctx,
        &mut mapping,
        &mut results,
        config.base.max_matches,
    );

    if config.base.uniquify {
        let mut seen = FxHashSet::default();
        results.retain(|m| {
            let mut key: Vec<u32> = m.values().map(|idx| idx.0).collect();
            key.sort_unstable();
            seen.insert(key)
        });
    }

    Ok((results, ctx.budget_exhausted.get()))
}

/// Existence-only search, mirroring [`crate::has_match_bounded`]'s 3-way
/// [`MatchOutcome`] contract on top of the RDKit-parity ring-count model.
pub fn has_match_rdkit_parity_bounded(
    query: &QueryMolecule,
    mol: &Molecule,
    config: &RdkitParityConfig,
) -> Result<MatchOutcome, RdkitParityError> {
    let mut one_match_cfg = config.clone();
    one_match_cfg.base.max_matches = Some(1);
    let (results, budget_exhausted) = find_matches_rdkit_parity(query, mol, &one_match_cfg)?;
    Ok(if !results.is_empty() {
        MatchOutcome::Found
    } else if budget_exhausted {
        MatchOutcome::BudgetExhausted
    } else {
        MatchOutcome::NotFound
    })
}

/// Scan a query (including nested recursive `$(...)` sub-queries) for any
/// use of `AtomPrimitive::RingCount` (`[RN]`) — the only predicate this
/// mode's ring model actually changes. Used to skip building the model
/// entirely (and its candidate-search budget) for queries that don't need
/// it, matching this crate's existing "don't pay for what you don't use"
/// precedent (`EvalCtx::min_ring_size_by_atom` in `match_vf2.rs`).
fn query_uses_ring_count(query: &QueryMolecule) -> bool {
    query
        .atoms
        .iter()
        .any(|a| atom_query_uses_ring_count(&a.query))
}

fn atom_query_uses_ring_count(q: &AtomQuery) -> bool {
    match q {
        AtomQuery::Primitive(AtomPrimitive::RingCount(_)) => true,
        AtomQuery::Primitive(AtomPrimitive::Recursive(sub)) => query_uses_ring_count(sub),
        AtomQuery::Primitive(_) => false,
        AtomQuery::And(a, b) | AtomQuery::Or(a, b) => {
            atom_query_uses_ring_count(a) || atom_query_uses_ring_count(b)
        }
        AtomQuery::Not(a) => atom_query_uses_ring_count(a),
    }
}

// ---------------------------------------------------------------------------
// Evaluation context -- same shape as match_vf2::EvalCtx, plus the ring model.
// ---------------------------------------------------------------------------

struct EvalCtx<'a> {
    mol: &'a Molecule,
    rings: &'a RingSet,
    ring_model: Option<&'a RdkitParityRingModel>,
    config: &'a MatchConfig,
    visit_budget: std::cell::Cell<u64>,
    budget_exhausted: std::cell::Cell<bool>,
    min_ring_size_by_atom: std::cell::RefCell<Option<Vec<Option<u8>>>>,
}

impl EvalCtx<'_> {
    fn min_ring_size(&self, idx: AtomIdx) -> Option<u8> {
        let mut cache = self.min_ring_size_by_atom.borrow_mut();
        let table = cache.get_or_insert_with(|| {
            let mut table = vec![None; self.mol.atom_count()];
            for ring in self.rings.rings() {
                let size = ring.len() as u8;
                for &atom in ring {
                    let slot = &mut table[atom.0 as usize];
                    *slot = Some(slot.map_or(size, |current: u8| current.min(size)));
                }
            }
            table
        });
        table[idx.0 as usize]
    }
}

// ---------------------------------------------------------------------------
// Recursive VF2 search -- identical control flow to match_vf2::match_recursive.
// ---------------------------------------------------------------------------

fn next_unmapped(mapping: &FxHashMap<usize, AtomIdx>, query_len: usize) -> usize {
    (0..query_len).find(|i| !mapping.contains_key(i)).unwrap()
}

fn match_recursive(
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
    mapping: &mut FxHashMap<usize, AtomIdx>,
    results: &mut Vec<FxHashMap<usize, AtomIdx>>,
    max: Option<usize>,
) {
    if max.is_some_and(|m| results.len() >= m) {
        return;
    }
    let remaining = ctx.visit_budget.get();
    if remaining == 0 {
        ctx.budget_exhausted.set(true);
        return;
    }
    ctx.visit_budget.set(remaining - 1);

    if mapping.len() == query.atoms.len() {
        results.push(mapping.clone());
        return;
    }

    let q_next = next_unmapped(mapping, query.atoms.len());
    let used_targets: FxHashSet<AtomIdx> = mapping.values().copied().collect();

    for t in 0..ctx.mol.atom_count() {
        if max.is_some_and(|m| results.len() >= m) {
            break;
        }
        let t_idx = AtomIdx(t as u32);
        if used_targets.contains(&t_idx) {
            continue;
        }
        if !eval_atom_query(&query.atoms[q_next].query, t_idx, ctx) {
            continue;
        }
        if !bonds_compatible(q_next, t_idx, mapping, query, ctx) {
            continue;
        }
        mapping.insert(q_next, t_idx);
        match_recursive(query, ctx, mapping, results, max);
        mapping.remove(&q_next);
    }
}

fn bonds_compatible(
    q: usize,
    t: AtomIdx,
    mapping: &FxHashMap<usize, AtomIdx>,
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
) -> bool {
    for &(bond_idx, q_nb) in &query.adj[q] {
        if let Some(&t_nb) = mapping.get(&q_nb) {
            match ctx.mol.bond_between(t, t_nb) {
                None => return false,
                Some((_bidx, bond_entry)) => {
                    let bq = &query.bonds[bond_idx].query;
                    if !eval_bond_query(bq, bond_entry.order, t, t_nb, ctx) {
                        return false;
                    }
                }
            }
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Atom query evaluation -- identical to match_vf2, except RingCount.
// ---------------------------------------------------------------------------

fn eval_atom_query(q: &AtomQuery, idx: AtomIdx, ctx: &EvalCtx<'_>) -> bool {
    match q {
        AtomQuery::Primitive(p) => eval_atom_primitive(p, idx, ctx),
        AtomQuery::And(a, b) => eval_atom_query(a, idx, ctx) && eval_atom_query(b, idx, ctx),
        AtomQuery::Or(a, b) => eval_atom_query(a, idx, ctx) || eval_atom_query(b, idx, ctx),
        AtomQuery::Not(a) => !eval_atom_query(a, idx, ctx),
    }
}

fn eval_atom_primitive(p: &AtomPrimitive, idx: AtomIdx, ctx: &EvalCtx<'_>) -> bool {
    let atom = ctx.mol.atom(idx);
    match p {
        AtomPrimitive::AtomicNum(n) => atom.element.atomic_number() == *n,
        AtomPrimitive::Symbol(s) => atom.element.symbol() == s.as_str(),
        AtomPrimitive::Aromatic(a) => atom.aromatic == *a,
        AtomPrimitive::Charge(c) => atom.charge == *c,
        AtomPrimitive::HCount(h) => eval_hcount(idx, ctx, *h),
        AtomPrimitive::ImplicitHCount(h) => implicit_hcount(ctx.mol, idx) == *h,
        AtomPrimitive::Degree(d) => ctx.mol.neighbors(idx).count() as u8 == *d,
        // [R]/[!R] -- provably invariant to which SSSR basis is used (see
        // `rdkit_ring_model`'s module doc comment) -- left on plain SSSR.
        AtomPrimitive::RingMembership(r) => ctx.rings.contains_atom(idx) == *r,
        // [kN] -- already ~100% (99.98%) against RDKit on plain SSSR
        // (SMARTS-R1); deliberately left unchanged here, see this module's
        // doc comment.
        AtomPrimitive::RingSize(n) => ctx
            .rings
            .rings()
            .iter()
            .any(|ring| ring.len() == *n as usize && ring.contains(&idx)),
        // [rN] -- same rationale as RingSize above.
        AtomPrimitive::MinRingSize(n) => ctx.min_ring_size(idx) == Some(*n),
        AtomPrimitive::Wildcard => true,
        AtomPrimitive::Recursive(sub_query) => has_match_anchored(sub_query, idx, ctx),
        AtomPrimitive::Valence(v) => eval_valence(idx, ctx, *v),
        // [xN] -- provably invariant, same rationale as RingMembership above.
        AtomPrimitive::RingBondCount(x) => eval_ring_bond_count(idx, ctx, *x),
        AtomPrimitive::TotalConnectivity(x) => {
            ctx.mol.neighbors(idx).count() as u8 + implicit_hcount(ctx.mol, idx) == *x
        }
        // [RN], N >= 1 -- the one primitive this mode actually changes.
        // Falls back to the plain-SSSR count if the query somehow reaches
        // here without the model having been built (shouldn't happen:
        // `query_uses_ring_count` scans for this exact primitive before
        // the model is constructed) -- documented fallback, not a silent
        // divergence, since it's identical to the default matcher's own
        // formula in that (unreachable in practice) case.
        AtomPrimitive::RingCount(n) => match ctx.ring_model {
            Some(model) => model.ring_count(idx) == *n,
            None => {
                ctx.rings
                    .rings()
                    .iter()
                    .filter(|r| r.contains(&idx))
                    .count() as u8
                    == *n
            }
        },
        AtomPrimitive::Hybridization(h) => eval_hybridization(idx, ctx, *h),
        AtomPrimitive::Isotope(mass) => {
            !ctx.config.use_isotopes || ctx.mol.atom(idx).isotope == Some(*mass)
        }
        AtomPrimitive::Chirality(kind) => eval_chirality(idx, ctx, *kind),
    }
}

fn eval_hcount(idx: AtomIdx, ctx: &EvalCtx<'_>, h: u8) -> bool {
    let explicit_h = ctx
        .mol
        .neighbors(idx)
        .filter(|(nb, _)| ctx.mol.atom(*nb).element.atomic_number() == 1)
        .count() as u8;
    explicit_h + implicit_hcount(ctx.mol, idx) == h
}

fn eval_valence(idx: AtomIdx, ctx: &EvalCtx<'_>, v: u8) -> bool {
    let bond_sum: u8 = ctx
        .mol
        .neighbors(idx)
        .map(|(_, bid)| bond_order_int(ctx.mol.bond(bid).order))
        .sum();
    bond_sum + implicit_hcount(ctx.mol, idx) == v
}

fn eval_ring_bond_count(idx: AtomIdx, ctx: &EvalCtx<'_>, x: u8) -> bool {
    let count = ctx
        .mol
        .neighbors(idx)
        .filter(|(nb, _)| {
            ctx.rings
                .rings()
                .iter()
                .any(|ring| ring.contains(&idx) && ring.contains(nb))
        })
        .count() as u8;
    count == x
}

fn eval_hybridization(idx: AtomIdx, ctx: &EvalCtx<'_>, h: u8) -> bool {
    let atom = ctx.mol.atom(idx);
    let hyb = if atom.aromatic {
        2u8
    } else {
        let mut has_triple = false;
        let mut has_double = false;
        for (_, bid) in ctx.mol.neighbors(idx) {
            match ctx.mol.bond(bid).order {
                BondOrder::Triple => {
                    has_triple = true;
                    break;
                }
                BondOrder::Double => has_double = true,
                _ => {}
            }
        }
        if has_triple {
            1
        } else if has_double {
            2
        } else {
            3
        }
    };
    hyb == h
}

fn eval_chirality(idx: AtomIdx, ctx: &EvalCtx<'_>, kind: u8) -> bool {
    if !ctx.config.use_chirality {
        return true;
    }
    use chematic_core::Chirality;
    let c = ctx.mol.atom(idx).chirality;
    match kind {
        1 => c == Chirality::CounterClockwise,
        2 => c == Chirality::Clockwise,
        _ => c != Chirality::None,
    }
}

// ---------------------------------------------------------------------------
// Anchored match helpers (for recursive SMARTS) -- identical to match_vf2.
// ---------------------------------------------------------------------------

fn has_match_anchored(query: &QueryMolecule, anchor: AtomIdx, ctx: &EvalCtx<'_>) -> bool {
    if query.atoms.is_empty() {
        return false;
    }
    if query.atoms.len() > ctx.mol.atom_count() {
        return false;
    }
    if !eval_atom_query(&query.atoms[0].query, anchor, ctx) {
        return false;
    }
    let mut mapping = FxHashMap::default();
    mapping.insert(0usize, anchor);
    if query.atoms.len() == 1 {
        return true;
    }
    has_match_recursive(query, ctx, &mut mapping)
}

fn has_match_recursive(
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
    mapping: &mut FxHashMap<usize, AtomIdx>,
) -> bool {
    let remaining = ctx.visit_budget.get();
    if remaining == 0 {
        ctx.budget_exhausted.set(true);
        return false;
    }
    ctx.visit_budget.set(remaining - 1);

    if mapping.len() == query.atoms.len() {
        return true;
    }

    let q_next = next_unmapped(mapping, query.atoms.len());
    let used_targets: FxHashSet<AtomIdx> = mapping.values().copied().collect();

    for t in 0..ctx.mol.atom_count() {
        let t_idx = AtomIdx(t as u32);
        if used_targets.contains(&t_idx) {
            continue;
        }
        if !eval_atom_query(&query.atoms[q_next].query, t_idx, ctx) {
            continue;
        }
        if !bonds_compatible(q_next, t_idx, mapping, query, ctx) {
            continue;
        }
        mapping.insert(q_next, t_idx);
        if has_match_recursive(query, ctx, mapping) {
            mapping.remove(&q_next);
            return true;
        }
        mapping.remove(&q_next);
    }
    false
}

// ---------------------------------------------------------------------------
// Bond query evaluation -- identical to match_vf2.
// ---------------------------------------------------------------------------

fn eval_bond_query(
    q: &BondQuery,
    order: BondOrder,
    a: AtomIdx,
    b: AtomIdx,
    ctx: &EvalCtx<'_>,
) -> bool {
    match q {
        BondQuery::Primitive(p) => eval_bond_primitive(p, order, a, b, ctx),
        BondQuery::And(x, y) => {
            eval_bond_query(x, order, a, b, ctx) && eval_bond_query(y, order, a, b, ctx)
        }
        BondQuery::Or(x, y) => {
            eval_bond_query(x, order, a, b, ctx) || eval_bond_query(y, order, a, b, ctx)
        }
        BondQuery::Not(x) => !eval_bond_query(x, order, a, b, ctx),
        BondQuery::Any => true,
    }
}

fn bond_order_int(order: BondOrder) -> u8 {
    match order {
        BondOrder::Zero => 0,
        BondOrder::Single
        | BondOrder::Up
        | BondOrder::Down
        | BondOrder::Aromatic
        | BondOrder::Dative
        | BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
    }
}

fn eval_bond_primitive(
    p: &BondPrimitive,
    order: BondOrder,
    a: AtomIdx,
    b: AtomIdx,
    ctx: &EvalCtx<'_>,
) -> bool {
    match p {
        BondPrimitive::Single => {
            matches!(
                order,
                BondOrder::Single
                    | BondOrder::Up
                    | BondOrder::Down
                    | BondOrder::QuerySingleOrDouble
                    | BondOrder::QuerySingleOrAromatic
            )
        }
        BondPrimitive::Double => matches!(
            order,
            BondOrder::Double | BondOrder::QuerySingleOrDouble | BondOrder::QueryDoubleOrAromatic
        ),
        BondPrimitive::Triple => matches!(order, BondOrder::Triple),
        BondPrimitive::Aromatic => matches!(
            order,
            BondOrder::Aromatic
                | BondOrder::QuerySingleOrAromatic
                | BondOrder::QueryDoubleOrAromatic
        ),
        BondPrimitive::Any => true,
        // Ring-bond `@`/`!@` -- provably invariant, same rationale as
        // RingMembership/RingBondCount above.
        BondPrimitive::Ring => ctx
            .rings
            .rings()
            .iter()
            .any(|ring| ring.contains(&a) && ring.contains(&b)),
        BondPrimitive::Up => matches!(order, BondOrder::Up),
        BondPrimitive::Down => matches!(order, BondOrder::Down),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_smarts;
    use chematic_smiles::parse;

    fn atom_sets(matches: &[FxHashMap<usize, AtomIdx>]) -> Vec<Vec<u32>> {
        let mut sets: Vec<Vec<u32>> = matches
            .iter()
            .map(|m| {
                let mut v: Vec<u32> = m.values().map(|a| a.0).collect();
                v.sort_unstable();
                v
            })
            .collect();
        sets.sort();
        sets
    }

    #[test]
    fn matches_default_on_ordinary_query() {
        // A query with no [RN] at all: opt-in mode must produce the exact
        // same match set as the default matcher (no ring model even built).
        let mol = parse("c1ccc2ccccc2c1").unwrap(); // naphthalene
        let query = parse_smarts("c1ccccc1").unwrap();
        let default = crate::find_matches(&query, &mol);
        let (parity, exhausted) =
            find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
        assert!(!exhausted);
        assert_eq!(atom_sets(&default), atom_sets(&parity));
    }

    #[test]
    fn r2_matches_extra_ring_bridgeheads_on_adamantane() {
        // Adamantane: RDKit ground truth (rdkit==2026.03.3, live oracle,
        // atom indices 0-9 in the same left-to-right SMILES parse order as
        // chematic's) gives exactly 4 atoms (indices 1,3,5,7) a ring count
        // of 3, the rest (0,2,4,6,8,9) a ring count of 2.
        //
        // chematic's raw SSSR (cycle_rank 3, no augmentation) happens to
        // pick a *different but equally valid* basis where atom 1 alone
        // already sits in all 3 basis rings -- a concrete illustration of
        // the "genuine SSSR-basis-cardinality disagreement" `docs/
        // rdkit_compat.md`'s SMARTS-R0/R2 sections describe: different
        // (both valid) basis choices produce different naive [R3] answers
        // even before any symmetrization. See `dbg_print_sssr` below for
        // the raw per-atom counts this depends on.
        let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap();
        let query = parse_smarts("[R3]").unwrap();
        let (parity, _) =
            find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
        let mut parity_atoms: Vec<u32> = parity.iter().map(|m| m[&0].0).collect();
        parity_atoms.sort_unstable();
        assert_eq!(parity_atoms, vec![1, 3, 5, 7], "RDKit-parity [R3] atom set");

        // Default matcher (plain SSSR) only finds atom 1 -- the one atom
        // chematic's own basis choice happens to already give count 3 to,
        // without needing the 4th "extra" ring at all. This is exactly the
        // partial/coincidental overlap the mode difference is meant to fix
        // (the other 3 -- atoms 3, 5, 7 -- are missed on plain SSSR).
        let default = crate::find_matches(&query, &mol);
        let mut default_atoms: Vec<u32> = default.iter().map(|m| m[&0].0).collect();
        default_atoms.sort_unstable();
        assert_eq!(default_atoms, vec![1], "default (plain-SSSR) [R3] atom set");
    }

    #[test]
    fn ring_model_budget_exceeded_is_typed_not_silent() {
        let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap(); // adamantane
        let query = parse_smarts("[R3]").unwrap();
        let config = RdkitParityConfig {
            ring_model_budget: RdkitRingModelBudget { max_candidates: 0 },
            ..RdkitParityConfig::default()
        };
        let result = find_matches_rdkit_parity(&query, &mol, &config);
        assert!(matches!(
            result,
            Err(RdkitParityError::RingModelBudgetExceeded { .. })
        ));
    }

    #[test]
    fn non_ring_count_query_unaffected_by_zero_ring_budget() {
        // A query that never touches [RN] must not even attempt to build
        // the ring model, so a starved budget must not affect it.
        let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap(); // adamantane
        let query = parse_smarts("[R]").unwrap();
        let config = RdkitParityConfig {
            ring_model_budget: RdkitRingModelBudget { max_candidates: 0 },
            ..RdkitParityConfig::default()
        };
        let (parity, _) = find_matches_rdkit_parity(&query, &mol, &config).unwrap();
        assert_eq!(parity.len(), 10); // every atom in adamantane is in a ring
    }

    #[test]
    fn has_match_rdkit_parity_bounded_three_way_outcome() {
        let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap();
        let query = parse_smarts("[R3]").unwrap();
        let found =
            has_match_rdkit_parity_bounded(&query, &mol, &RdkitParityConfig::default()).unwrap();
        assert_eq!(found, MatchOutcome::Found);

        let query_none = parse_smarts("[R5]").unwrap();
        let not_found =
            has_match_rdkit_parity_bounded(&query_none, &mol, &RdkitParityConfig::default())
                .unwrap();
        assert_eq!(not_found, MatchOutcome::NotFound);
    }

    // -- SMARTS-A0 bridgehead-N bucket: this pipeline vs. RdkitLike re-perception --
    //
    // This reproducer guards the bridgehead-N ring-fusion aromaticity case.
    // RDKit ground truth is that only atoms 2-7, the benzo ring, are
    // aromatic. Both the direct parser and explicit RdkitLike re-perception
    // must preserve that result.

    #[test]
    fn smarts_a0_does_not_fire_on_direct_parse_no_reperception() {
        // Neither `find_matches` nor `find_matches_rdkit_parity` ever calls
        // any aromaticity re-perception on their own -- they match whatever
        // flags the input molecule already carries. The reproducer's mixed
        // aromatic/Kekule SMILES already encodes the *correct* (RDKit-
        // matching) flags directly from parsing, so SMARTS-A0's precondition
        // (a `RdkitLike` re-perception pass) is simply never reached on this
        // pipeline as used by this crate's default entry points.
        let mol = parse("C1=Cc2ccccc2C2=NCCCN12").unwrap();
        let c_query = parse_smarts("c").unwrap();
        let default = crate::find_matches(&c_query, &mol);
        let mut atoms: Vec<u32> = default.iter().map(|m| m[&0].0).collect();
        atoms.sort_unstable();
        // RDKit ground truth (rdkit==2026.03.3, live oracle): exactly atoms
        // 2-7 (the benzo ring), never atoms 0, 1, or 13.
        assert_eq!(atoms, vec![2, 3, 4, 5, 6, 7]);
    }

    #[test]
    fn smarts_a0_remains_rdkit_compatible_after_explicit_reperception() {
        // Callers may explicitly re-perceive with
        // `AromaticityAlgorithm::RdkitLike` before handing the molecule to
        // this crate. The fixed aromaticity path must retain the RDKit
        // ground-truth benzo ring rather than extend it into the fused ring.
        let mol = parse("C1=Cc2ccccc2C2=NCCCN12").unwrap();
        let kekulized = chematic_core::kekulize(&mol)
            .map(|k| chematic_core::apply_kekule(&mol, &k))
            .unwrap_or(mol);
        let reperceived = chematic_perception::apply_aromaticity_ex(
            &kekulized,
            chematic_perception::AromaticityAlgorithm::RdkitLike,
        );
        let c_query = parse_smarts("c").unwrap();
        let default = crate::find_matches(&c_query, &reperceived);
        let mut atoms: Vec<u32> = default.iter().map(|m| m[&0].0).collect();
        atoms.sort_unstable();
        assert_eq!(atoms, vec![2, 3, 4, 5, 6, 7]);
    }

    #[test]
    fn use_rdkit_parity_aromaticity_flag_avoids_the_rdkit_like_bug_on_this_reproducer() {
        // This crate's OTHER opt-in flag, `use_rdkit_parity_aromaticity`
        // (a from-scratch port of RDKit's actual aromaticity algorithm, see
        // `chematic_perception::apply_aromaticity_rdkit_parity_experimental`'s
        // own doc comment) is a materially different engine from the
        // heuristic `AromaticityAlgorithm::RdkitLike` used above. On this
        // specific bare-core reproducer it does NOT reproduce SMARTS-A0's
        // over-extension -- checked directly, not assumed, since the two
        // engines are unrelated code paths and this crate must not imply
        // one fixes the other without evidence.
        let mol = parse("C1=Cc2ccccc2C2=NCCCN12").unwrap();
        let c_query = parse_smarts("c").unwrap();
        let config = RdkitParityConfig {
            use_rdkit_parity_aromaticity: true,
            ..RdkitParityConfig::default()
        };
        let (parity, _) = find_matches_rdkit_parity(&c_query, &mol, &config).unwrap();
        let mut atoms: Vec<u32> = parity.iter().map(|m| m[&0].0).collect();
        atoms.sort_unstable();
        assert_eq!(
            atoms,
            vec![2, 3, 4, 5, 6, 7],
            "use_rdkit_parity_aromaticity=true should match RDKit's real answer \
             (benzo ring only) on this reproducer, not RdkitLike's over-extension"
        );
    }
}
