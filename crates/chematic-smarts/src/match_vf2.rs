//! VF2 subgraph isomorphism: find all embeddings of a `QueryMolecule` in a target `Molecule`.
//!
//! The classic VF2 algorithm explores a state-space search tree.  At each step
//! it picks the next unmapped query atom and tries to extend the current partial
//! mapping with every compatible target atom.  Compatibility is checked at two
//! levels:
//!
//! 1. **Atom compatibility** — the target atom must satisfy the query atom's
//!    `AtomQuery` expression.
//! 2. **Bond compatibility** — for every already-mapped query neighbour of the
//!    candidate query atom, the corresponding target atoms must be bonded in the
//!    target molecule, and that bond must satisfy the query bond's `BondQuery`
//!    expression.

use std::collections::HashMap;

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};
use chematic_perception::find_sssr;
use chematic_perception::RingSet;

use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

// ---------------------------------------------------------------------------
// Evaluation context (precomputed per `find_matches` call)
// ---------------------------------------------------------------------------

/// Per-call evaluation context: the target molecule and its precomputed ring set.
///
/// Computing the SSSR is expensive; this struct ensures it is done once per
/// `find_matches` call, not once per atom evaluation.
struct EvalCtx<'a> {
    mol: &'a Molecule,
    rings: RingSet,
}

impl<'a> EvalCtx<'a> {
    fn new(mol: &'a Molecule) -> Self {
        Self { mol, rings: find_sssr(mol) }
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Find all non-overlapping (injective) embeddings of `query` in `mol`.
///
/// Returns a `Vec` of mappings, each mapping a query atom index to a target
/// `AtomIdx`.  Each individual mapping is injective (no two query atoms map to
/// the same target atom), but the same target atom may appear in different
/// mappings.
pub fn find_matches(
    query: &QueryMolecule,
    mol: &Molecule,
) -> Vec<HashMap<usize, AtomIdx>> {
    if query.atoms.is_empty() {
        return vec![];
    }

    let ctx = EvalCtx::new(mol);
    let mut mapping: HashMap<usize, AtomIdx> = HashMap::new();
    let mut results: Vec<HashMap<usize, AtomIdx>> = Vec::new();

    match_recursive(query, &ctx, &mut mapping, &mut results);
    results
}

// ---------------------------------------------------------------------------
// Recursive VF2 search
// ---------------------------------------------------------------------------

fn match_recursive(
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
    mapping: &mut HashMap<usize, AtomIdx>,
    results: &mut Vec<HashMap<usize, AtomIdx>>,
) {
    // Base case: all query atoms have been mapped.
    if mapping.len() == query.atoms.len() {
        results.push(mapping.clone());
        return;
    }

    // Pick the next unmapped query atom (smallest index not yet in mapping).
    let q_next = (0..query.atoms.len())
        .find(|i| !mapping.contains_key(i))
        .unwrap(); // safe: mapping.len() < query.atoms.len()

    // Collect the set of target atoms already used in this mapping so we can
    // enforce injectivity.
    let used_targets: std::collections::HashSet<AtomIdx> = mapping.values().copied().collect();

    // Try each target atom as a candidate for q_next.
    for t in 0..ctx.mol.atom_count() {
        let t_idx = AtomIdx(t as u32);

        // 1. Injectivity: target atom must not already be mapped.
        if used_targets.contains(&t_idx) {
            continue;
        }

        // 2. Atom query must match.
        if !eval_atom_query(&query.atoms[q_next].query, t_idx, ctx) {
            continue;
        }

        // 3. Bond constraints from already-mapped neighbours of q_next.
        if !bonds_compatible(q_next, t_idx, mapping, query, ctx) {
            continue;
        }

        // Extend the mapping and recurse.
        mapping.insert(q_next, t_idx);
        match_recursive(query, ctx, mapping, results);
        mapping.remove(&q_next);
    }
}

// ---------------------------------------------------------------------------
// Bond compatibility check
// ---------------------------------------------------------------------------

/// For every already-mapped query neighbour of `q`, verify that the
/// corresponding target atoms are bonded and satisfy the query bond condition.
fn bonds_compatible(
    q: usize,
    t: AtomIdx,
    mapping: &HashMap<usize, AtomIdx>,
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
) -> bool {
    for &(bond_idx, q_nb) in &query.adj[q] {
        // Only check neighbours that are already mapped.
        if let Some(&t_nb) = mapping.get(&q_nb) {
            // The target must have a bond between t and t_nb.
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
// Atom query evaluation
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
        AtomPrimitive::HCount(h) => implicit_hcount(ctx.mol, idx) == *h,
        AtomPrimitive::Degree(d) => ctx.mol.neighbors(idx).count() as u8 == *d,
        AtomPrimitive::RingMembership(r) => ctx.rings.contains_atom(idx) == *r,
        AtomPrimitive::RingSize(n) => ctx.rings.rings().iter().any(|ring| {
            ring.len() == *n as usize && ring.contains(&idx)
        }),
        AtomPrimitive::Wildcard => true,
        AtomPrimitive::Recursive(sub_query) => {
            // The target atom `idx` must be the root of at least one embedding
            // of `sub_query` in the same target molecule.
            has_match_anchored(sub_query, idx, ctx)
        }
        AtomPrimitive::Valence(v) => {
            // Total valence = sum of explicit bond orders + implicit H count.
            let bond_sum: u8 = ctx.mol.neighbors(idx).map(|(_, bid)| {
                match ctx.mol.bond(bid).order {
                    BondOrder::Single | BondOrder::Up | BondOrder::Down => 1u8,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    BondOrder::Quadruple => 4,
                    BondOrder::Aromatic => 1, // aromatic bond counted as 1
                }
            }).sum();
            bond_sum + implicit_hcount(ctx.mol, idx) == *v
        }
        AtomPrimitive::RingBondCount(x) => {
            // Count bonds where both endpoints share at least one SSSR ring.
            let count = ctx.mol.neighbors(idx).filter(|(nb, _)| {
                ctx.rings.rings().iter().any(|ring| ring.contains(&idx) && ring.contains(nb))
            }).count() as u8;
            count == *x
        }
        AtomPrimitive::TotalConnectivity(x) => {
            // Total connectivity = heavy-atom degree + implicit H count.
            let deg = ctx.mol.neighbors(idx).count() as u8;
            deg + implicit_hcount(ctx.mol, idx) == *x
        }
        AtomPrimitive::RingCount(n) => {
            // Count how many SSSR rings contain this atom.
            let count = ctx.rings.rings().iter().filter(|ring| ring.contains(&idx)).count() as u8;
            count == *n
        }
        AtomPrimitive::Hybridization(h) => {
            // Inferred hybridization:
            //   aromatic atom → sp2 (2)
            //   atom with any triple bond → sp (1)
            //   atom with any double bond → sp2 (2)
            //   otherwise → sp3 (3)
            let atom = ctx.mol.atom(idx);
            let hyb = if atom.aromatic {
                2u8
            } else if ctx.mol.neighbors(idx).any(|(_, bid)| {
                matches!(ctx.mol.bond(bid).order, BondOrder::Triple)
            }) {
                1
            } else if ctx.mol.neighbors(idx).any(|(_, bid)| {
                matches!(ctx.mol.bond(bid).order, BondOrder::Double)
            }) {
                2
            } else {
                3
            };
            hyb == *h
        }
    }
}

// ---------------------------------------------------------------------------
// Anchored match helpers (for recursive SMARTS)
// ---------------------------------------------------------------------------

/// Returns `true` if there exists at least one embedding of `query` in `ctx.mol`
/// with query atom 0 forced to map to `anchor`.
fn has_match_anchored(query: &QueryMolecule, anchor: AtomIdx, ctx: &EvalCtx<'_>) -> bool {
    if query.atoms.is_empty() {
        return false;
    }
    // Quick check: anchor must satisfy the first query atom.
    if !eval_atom_query(&query.atoms[0].query, anchor, ctx) {
        return false;
    }
    // Seed the mapping with query atom 0 → anchor.
    let mut mapping = HashMap::new();
    mapping.insert(0usize, anchor);
    // Single-atom query — the anchor already satisfies it.
    if query.atoms.len() == 1 {
        return true;
    }
    has_match_recursive(query, ctx, &mut mapping)
}

/// Depth-first search for a complete embedding, starting from a partial
/// `mapping`.  Returns as soon as the first complete match is found.
fn has_match_recursive(
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
    mapping: &mut HashMap<usize, AtomIdx>,
) -> bool {
    // Base case: all query atoms mapped.
    if mapping.len() == query.atoms.len() {
        return true;
    }

    // Pick the next unmapped query atom.
    let q_next = (0..query.atoms.len())
        .find(|i| !mapping.contains_key(i))
        .unwrap();

    let used_targets: std::collections::HashSet<AtomIdx> = mapping.values().copied().collect();

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
// Bond query evaluation
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
        // Implicit "any bond" — matches any bond order.
        BondQuery::Any => true,
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
            matches!(order, BondOrder::Single | BondOrder::Up | BondOrder::Down)
        }
        BondPrimitive::Double => matches!(order, BondOrder::Double),
        BondPrimitive::Triple => matches!(order, BondOrder::Triple),
        BondPrimitive::Aromatic => matches!(order, BondOrder::Aromatic),
        BondPrimitive::Any => true,
        BondPrimitive::Ring => {
            // A bond is a "ring bond" if both its endpoints share at least one common ring.
            ctx.rings.rings().iter().any(|ring| ring.contains(&a) && ring.contains(&b))
        }
    }
}
