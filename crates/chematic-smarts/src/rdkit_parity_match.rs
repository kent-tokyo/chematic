//! Opt-in RDKit-parity SMARTS matching mode.
//!
//! This module is a **deliberate, near-total duplication** of
//! [`crate::match_vf2`]'s VF2 recursive matcher, not a refactor of it. The
//! duplication is intentional: it keeps `match_vf2.rs` at **zero diff**, so
//! "the default matcher is byte-identical before and after this change" is
//! trivially provable (there is no change to prove anything about) rather
//! than something that has to be argued from a shared-code refactor. The
//! evaluator diverges from `match_vf2.rs` for `AtomPrimitive::RingCount`
//! (`[RN]`) and `AtomPrimitive::RingSize` (`[kN]`). See
//! [`crate::rdkit_ring_model`] for the ring-count model. Ring-size membership
//! uses a separately bounded symmetrized-SSSR selection. The narrow
//! Fe--[C-] cleanup below also makes `[xN]` and ring-bond `@`/`!@` exclude
//! a converted dative bond. Other ring primitives retain the plain-SSSR
//! formula from `match_vf2.rs`.
//!
//! Everything else in this file — chirality/isotope handling, valence,
//! hybridization, recursive-SMARTS anchoring, the visit-budget/
//! `MatchOutcome` contract — is copied verbatim in behavior from
//! `match_vf2.rs` so this mode's other results are identical to
//! the default matcher's.

use rustc_hash::{FxHashMap, FxHashSet};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Element, Molecule, implicit_hcount};
use chematic_perception::{RingSet, SymmetrizedSssrStatus};

use crate::match_vf2::hetero_neighbor_count;
use crate::match_vf2::{MatchConfig, MatchOutcome};
use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};
use crate::rdkit_ring_model::{
    RdkitParityError, RdkitParityRingModel, RdkitRingModelBudget, build_rdkit_parity_ring_model,
    build_shared_symmetrized_ring_model,
};

/// Configuration for [`find_matches_rdkit_parity`] / [`has_match_rdkit_parity_bounded`].
#[derive(Debug, Clone, Default)]
pub struct RdkitParityConfig {
    /// Everything [`MatchConfig`] already covers (chirality/isotope
    /// enforcement, `max_matches`, `uniquify`, VF2 `max_visit_budget`) —
    /// unaffected by, and orthogonal to, this mode's ring-count model.
    pub base: MatchConfig,
    /// Resource bound on the RDKit-parity ring-count and ring-size candidate
    /// searches. See [`RdkitRingModelBudget`].
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
    /// Use perception's alternate bounded symmetrized-SSSR selector for
    /// `[RN]` queries. This is an experimental comparison lane; the default
    /// remains the SMARTS-specific selector because the two models differ on
    /// bridged-cage fixtures.
    pub use_shared_symmetrized_sssr: bool,
}

/// Find all non-overlapping (injective) embeddings of `query` in `mol` using
/// the opt-in RDKit-parity ring models for `[RN]` and `[kN]`.
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

    let view;
    let mol_ref: &Molecule = if config.use_rdkit_parity_aromaticity {
        view = chematic_perception::apply_aromaticity_rdkit_parity_shared(mol);
        view.as_ref()
            .as_ref()
            .map_err(|e| RdkitParityError::Aromaticity(e.clone()))?
    } else {
        mol
    };

    // RDKit sanitization converts a hypervalent carbanion--metal single
    // bond into a dative bond before ring perception. Keep this narrowly
    // established Fe case inside the opt-in matcher: the caller's molecule,
    // native SMARTS, and the general SMILES parser remain unchanged.
    let organometallic_view = rdkit_parity_iron_carbanion_view(mol_ref);
    let mol_ref = organometallic_view.as_ref().unwrap_or(mol_ref);

    let rings = chematic_perception::find_sssr(mol_ref);
    let uses_ring_size = query_uses_ring_size(query);
    let symmetrized_rings = if uses_ring_size {
        let result = chematic_perception::find_symmetrized_sssr_with_diagnostics_bounded(
            mol_ref,
            Some(config.ring_model_budget.max_candidates),
        );
        if result.status() == SymmetrizedSssrStatus::CapExhausted {
            return Err(RdkitParityError::RingModelBudgetExceeded {
                candidates_examined: result.candidates_examined(),
                cap: config.ring_model_budget.max_candidates,
            });
        }
        Some(result.into_ring_set())
    } else {
        None
    };
    let size_model = if uses_ring_size {
        Some(build_rdkit_parity_ring_model(
            mol_ref,
            &rings,
            &config.ring_model_budget,
        )?)
    } else {
        None
    };
    let ring_model = if query_uses_ring_count(query) {
        let model = if config.use_shared_symmetrized_sssr {
            build_shared_symmetrized_ring_model(mol_ref, &rings, &config.ring_model_budget)?
        } else {
            build_rdkit_parity_ring_model(mol_ref, &rings, &config.ring_model_budget)?
        };
        let aromatic_cations = mol_ref
            .atoms()
            .filter(|(_, atom)| atom.aromatic && atom.charge > 0)
            .count();
        if aromatic_cations >= 2 && model.extra_ring_count() > 0 {
            return Err(RdkitParityError::RingModelAmbiguous {
                aromatic_cations,
                extra_rings: model.extra_ring_count(),
            });
        }
        Some(model)
    } else {
        None
    };

    let ctx = EvalCtx {
        mol: mol_ref,
        rings: &rings,
        symmetrized_rings: symmetrized_rings.as_ref(),
        size_model: size_model.as_ref(),
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

/// Reproduce the verified Fe--[C-] valence-four cleanup case without
/// pretending to implement RDKit's full organometallic sanitizer. RDKit
/// 2026.03.6 converts this single bond to dative before its ring search.
/// We only change a uniquely identified Fe bond; multiple metal choices
/// require canonical ranking and are deliberately left for separate work.
fn rdkit_parity_iron_carbanion_view(mol: &Molecule) -> Option<Molecule> {
    let mut convert: Vec<BondIdx> = Vec::new();
    for (idx, atom) in mol.atoms() {
        if atom.element != Element::C
            || atom.charge != -1
            || atom.aromatic
            || atom.hydrogen_count != Some(0)
        {
            continue;
        }
        // Keep this narrow: four explicit single bonds, no H, and exactly
        // one of those singles to Fe. Other valence/metal combinations need
        // their own oracle adjudication before changing the private view.
        let bonds: Vec<_> = mol.neighbors(idx).collect();
        if bonds.len() != 4
            || bonds
                .iter()
                .any(|(_, bond)| mol.bond(*bond).order != BondOrder::Single)
        {
            continue;
        }
        let metal_bonds: Vec<_> = bonds
            .iter()
            .filter(|(neighbor, _)| mol.atom(*neighbor).element == Element::FE)
            .map(|(_, bond)| bond)
            .collect();
        if metal_bonds.len() == 1 {
            convert.push(*metal_bonds[0]);
        }
    }
    if convert.is_empty() {
        return None;
    }
    let mut view = mol.clone();
    for bond in convert {
        view.set_bond_order(bond, BondOrder::Dative);
    }
    Some(view)
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
        // `[R0]` is ring-membership negation, not an exact positive ring
        // count. It is invariant under adding RDKit's extra SSSR rings and
        // must not trigger the ambiguity gate by itself.
        AtomQuery::Primitive(AtomPrimitive::RingCount(n)) => *n > 0,
        AtomQuery::Primitive(AtomPrimitive::Recursive(sub)) => query_uses_ring_count(sub),
        AtomQuery::Primitive(_) => false,
        AtomQuery::And(a, b) | AtomQuery::Or(a, b) => {
            atom_query_uses_ring_count(a) || atom_query_uses_ring_count(b)
        }
        AtomQuery::Not(a) => atom_query_uses_ring_count(a),
    }
}

fn query_uses_ring_size(query: &QueryMolecule) -> bool {
    query
        .atoms
        .iter()
        .any(|atom| atom_query_uses_ring_size(&atom.query))
}

fn atom_query_uses_ring_size(q: &AtomQuery) -> bool {
    match q {
        AtomQuery::Primitive(AtomPrimitive::RingSize(_)) => true,
        AtomQuery::Primitive(AtomPrimitive::Recursive(sub)) => query_uses_ring_size(sub),
        AtomQuery::Primitive(_) => false,
        AtomQuery::And(a, b) | AtomQuery::Or(a, b) => {
            atom_query_uses_ring_size(a) || atom_query_uses_ring_size(b)
        }
        AtomQuery::Not(a) => atom_query_uses_ring_size(a),
    }
}

// ---------------------------------------------------------------------------
// Evaluation context -- same shape as match_vf2::EvalCtx, plus the ring model.
// ---------------------------------------------------------------------------

struct EvalCtx<'a> {
    mol: &'a Molecule,
    rings: &'a RingSet,
    symmetrized_rings: Option<&'a RingSet>,
    size_model: Option<&'a RdkitParityRingModel>,
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
                    let qbond = &query.bonds[bond_idx];
                    // Whether the target bond runs in the query bond's own
                    // atom1 -> atom2 direction (dative `->` / `<-`).
                    let image_of_atom1 = if qbond.atom1 == q { t } else { t_nb };
                    let forward = bond_entry.atom1 == image_of_atom1;
                    if !eval_bond_query(&qbond.query, bond_entry.order, t, t_nb, forward, ctx) {
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
        // [kN] means membership in any selected ring of size N. RDKit's
        // default ring information is symmetrized; plain SSSR can miss a
        // symmetry-equivalent ring containing this atom. This opt-in lane
        // uses the bounded symmetrized set and refuses when it cannot be
        // computed completely.
        AtomPrimitive::RingSize(n) => {
            ctx.symmetrized_rings
                .unwrap_or(ctx.rings)
                .rings()
                .iter()
                .any(|ring| ring.len() == *n as usize && ring.contains(&idx))
                || ctx
                    .size_model
                    .is_some_and(|model| model.has_extra_ring_of_size(idx, *n as usize))
        }
        // [rN] -- same rationale as RingSize above.
        AtomPrimitive::MinRingSize(n) => ctx.min_ring_size(idx) == Some(*n),
        AtomPrimitive::Wildcard => true,
        AtomPrimitive::Recursive(sub_query) => has_match_anchored(sub_query, idx, ctx),
        AtomPrimitive::Valence(v) => eval_valence(idx, ctx, *v),
        // [xN] -- use the same ring basis, but never count a dative edge.
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
            !ctx.config.use_isotopes || ctx.mol.atom(idx).isotope.unwrap_or(0) == *mass
        }
        AtomPrimitive::Chirality(kind) => eval_chirality(idx, ctx, *kind),
        AtomPrimitive::HeteroNeighborCount(n) => {
            hetero_neighbor_count(ctx.mol, idx, false) == usize::from(*n)
        }
        AtomPrimitive::AliphaticHeteroNeighborCount(n) => {
            hetero_neighbor_count(ctx.mol, idx, true) == usize::from(*n)
        }
        AtomPrimitive::HeavyDegree(n) => {
            ctx.mol
                .neighbors(idx)
                .filter(|(nb, _)| ctx.mol.atom(*nb).element.atomic_number() != 1)
                .count()
                == usize::from(*n)
        }
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
    crate::match_vf2::total_valence(ctx.mol, idx) == v
}

fn eval_ring_bond_count(idx: AtomIdx, ctx: &EvalCtx<'_>, x: u8) -> bool {
    let count = ctx
        .mol
        .neighbors(idx)
        .filter(|(nb, bond)| {
            ctx.mol.bond(*bond).order != BondOrder::Dative
                && ctx
                    .rings
                    .rings()
                    .iter()
                    .any(|ring| ring.contains(&idx) && ring.contains(nb))
        })
        .count() as u8;
    count == x
}

/// RDKit's hybridization (`^n`), see [`crate::hybridization`].
fn eval_hybridization(idx: AtomIdx, ctx: &EvalCtx<'_>, h: u8) -> bool {
    crate::hybridization::rdkit_hybridization(ctx.mol, idx) == Some(h)
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
    forward: bool,
    ctx: &EvalCtx<'_>,
) -> bool {
    match q {
        BondQuery::Primitive(p) => eval_bond_primitive(p, order, a, b, forward, ctx),
        BondQuery::And(x, y) => {
            eval_bond_query(x, order, a, b, forward, ctx)
                && eval_bond_query(y, order, a, b, forward, ctx)
        }
        BondQuery::Or(x, y) => {
            eval_bond_query(x, order, a, b, forward, ctx)
                || eval_bond_query(y, order, a, b, forward, ctx)
        }
        BondQuery::Not(x) => !eval_bond_query(x, order, a, b, forward, ctx),
        // Unspecified SMARTS bond: single or aromatic (RDKit semantics).
        BondQuery::Any => matches!(
            order,
            BondOrder::Single
                | BondOrder::Up
                | BondOrder::Down
                | BondOrder::Aromatic
                | BondOrder::QuerySingleOrAromatic
        ),
    }
}

fn eval_bond_primitive(
    p: &BondPrimitive,
    order: BondOrder,
    a: AtomIdx,
    b: AtomIdx,
    forward: bool,
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
        // Dative bonds are excluded from the ring model, even when both
        // endpoints happen to lie in another ring of that model.
        BondPrimitive::Ring => {
            order != BondOrder::Dative
                && ctx
                    .rings
                    .rings()
                    .iter()
                    .any(|ring| ring.contains(&a) && ring.contains(&b))
        }
        // `/` and `\\` match a single or aromatic bond without constraining
        // cis/trans, as in RDKit (which ignores bond stereo when matching):
        // `C/C` matches ethane, and the answer for `F/C=C/F` does not depend
        // on how the target was written.
        BondPrimitive::Up | BondPrimitive::Down => matches!(
            order,
            BondOrder::Single
                | BondOrder::Up
                | BondOrder::Down
                | BondOrder::Aromatic
                | BondOrder::QuerySingleOrDouble
                | BondOrder::QuerySingleOrAromatic
        ),
        // Dative bond: the target's donor (its `atom1`) is the image of the
        // query bond's atom1 for `->`, of its atom2 for `<-`.
        BondPrimitive::DativeForward => order == BondOrder::Dative && forward,
        BondPrimitive::DativeBackward => order == BondOrder::Dative && !forward,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_smarts;
    use chematic_smiles::parse;

    #[test]
    fn organometallic_ring_closure_is_single_in_raw_graph() {
        let mol = parse("CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10").unwrap();
        let (bond_idx, bond) = mol.bond_between(AtomIdx(4), AtomIdx(9)).unwrap();
        assert_eq!(bond.order, BondOrder::Single);
        assert_eq!(mol.atom(AtomIdx(4)).charge, -1);
        assert_eq!(mol.atom(AtomIdx(9)).element.atomic_number(), 26);
        assert_eq!(
            mol.neighbors(AtomIdx(4))
                .map(|(_, idx)| mol.bond(idx).order.order_int() as usize)
                .sum::<usize>(),
            4
        );
        assert_eq!(mol.bond(bond_idx).order, BondOrder::Single);

        let view = rdkit_parity_iron_carbanion_view(&mol).unwrap();
        assert_eq!(view.bond(bond_idx).order, BondOrder::Dative);
        assert_eq!(mol.bond(bond_idx).order, BondOrder::Single);
        let (ordinary_fe_bond, _) = mol.bond_between(AtomIdx(9), AtomIdx(10)).unwrap();
        assert_eq!(view.bond(ordinary_fe_bond).order, BondOrder::Single);
        let rings = chematic_perception::find_sssr(&view);
        assert_eq!(
            rings
                .rings()
                .iter()
                .filter(|ring| ring.contains(&AtomIdx(4)))
                .count(),
            1
        );

        for (pattern, should_contain) in [
            ("[R1]", true),
            ("[R2]", false),
            ("[x2]", true),
            ("[x3]", false),
        ] {
            let query = parse_smarts(pattern).unwrap();
            let (matches, exhausted) =
                find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
            assert!(!exhausted);
            assert_eq!(
                matches
                    .iter()
                    .any(|mapping| mapping.values().any(|idx| *idx == AtomIdx(4))),
                should_contain,
                "{pattern}"
            );
        }

        for (pattern, should_contain) in [("*@*", false), ("*!@*", true)] {
            let query = parse_smarts(pattern).unwrap();
            let (matches, exhausted) =
                find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
            assert!(!exhausted);
            assert_eq!(
                atom_sets(&matches).contains(&vec![4, 9]),
                should_contain,
                "{pattern}"
            );
        }
    }

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
    fn k6_uses_symmetry_equivalent_ring_in_opt_in_mode() {
        // In this bridged fragment, atom 46 belongs to the six-membered
        // ring selected by RDKit 2026.03.6, but not to chematic's plain
        // SSSR basis. The default matcher must retain its existing answer.
        let mol = parse("C=CC[C@H](NC(=O)[C@@H]1C[C@@H](CCCc2cccc3ccccc23)c2c(Cl)nc(NCc3cccc(OC)c3)c(=O)n21)B1OC2CC3CC(C3(C)C)[C@@]2(C)O1").unwrap();
        let query = parse_smarts("[k6]").unwrap();
        let default = crate::find_matches(&query, &mol);
        let (parity, exhausted) =
            find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
        assert!(!exhausted);
        assert!(!default.iter().any(|m| m[&0] == AtomIdx(46)));
        assert!(parity.iter().any(|m| m[&0] == AtomIdx(46)));
    }

    #[test]
    fn k6_corpus_3498_uses_mixed_size_replacement() {
        // RDKit 2026.03.6 includes atom 22 in the replacement six-ring
        // 1-2-20-21-22-23. Perception's symmetrized set keeps only the
        // four-ring through atom 22, while the SMARTS-specific bounded
        // selector independently recovers the six-ring.
        let mol = parse("C=C1[C@H]2Oc3cc(C(C)(C)CCCCCC)cc(O)c3[C@H]2[C@H]2C[C@@H]1C2(C)C").unwrap();
        let query = parse_smarts("[k6]").unwrap();
        let (parity, _) =
            find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
        let mut atoms: Vec<u32> = parity.iter().map(|m| m[&0].0).collect();
        atoms.sort_unstable();
        assert_eq!(atoms, vec![1, 2, 4, 5, 6, 16, 17, 19, 20, 21, 22, 23, 24]);
    }

    #[test]
    fn k6_ring_model_budget_exceeded_is_typed() {
        let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap();
        let query = parse_smarts("[k6]").unwrap();
        let config = RdkitParityConfig {
            ring_model_budget: RdkitRingModelBudget { max_candidates: 0 },
            ..RdkitParityConfig::default()
        };
        assert!(matches!(
            find_matches_rdkit_parity(&query, &mol, &config),
            Err(RdkitParityError::RingModelBudgetExceeded { cap: 0, .. })
        ));
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
            Err(RdkitParityError::RingModelBudgetExceeded {
                candidates_examined: 1,
                cap: 0,
            })
        ));
    }

    #[test]
    fn shared_ring_model_propagates_measured_budget_exhaustion() {
        let mol = parse("C12C3C4C1C5C4C3C25").unwrap();
        let query = parse_smarts("[R3]").unwrap();
        let config = RdkitParityConfig {
            use_shared_symmetrized_sssr: true,
            ring_model_budget: RdkitRingModelBudget { max_candidates: 0 },
            ..RdkitParityConfig::default()
        };
        let result = find_matches_rdkit_parity(&query, &mol, &config);
        assert!(matches!(
            result,
            Err(RdkitParityError::RingModelBudgetExceeded {
                candidates_examined: 1,
                cap: 0,
            })
        ));
    }

    #[test]
    fn ambiguous_aromatic_cation_ring_model_fails_closed() {
        // This highly charged fused polyaromatic scaffold was the sole
        // SMARTS parity regression in the 5,021-molecule A4 census: the
        // generated replacement rings produced many false-positive [R3]
        // matches where RDKit returned only atoms 3 and 4.  Refuse the
        // ambiguous model instead of returning a plausible but wrong set.
        let mol =
            parse("c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1")
                .unwrap();
        let query = parse_smarts("[R3]").unwrap();
        let result = find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default());
        assert!(matches!(
            result,
            Err(RdkitParityError::RingModelAmbiguous {
                aromatic_cations: 2..,
                extra_rings: 1..
            })
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
