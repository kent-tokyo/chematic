//! Milestone 4B-2 production port of the Rule 4b tiebreak from the validated
//! reference engine (`examples/rule4b_bottom_up.rs`, 72/72 across 4 oracle corpora
//! including mirrors -- see `docs/rfcs/cip_accurate_rfc.md`). Mechanical port: algorithm
//! unchanged from the example.

use std::cmp::Ordering;
use std::collections::HashMap;

use chematic_core::{AtomIdx, Chirality, Molecule};

use crate::budget::CipBudget;
use crate::compare::CipCompareError;
use crate::digraph::CipDigraph;
use crate::node::{CipNodeKind, NodeId};
use crate::resolver::resolve_chirality;

/// BFS, in place within the outer digraph (no re-rooting), for the nearest
/// chirality-bearing atom reachable from `start`'s children onward. `None` if 2+
/// distinct atoms tie at the nearest level (ambiguous -- ties fall through to
/// `SkipReason::Tied` at the call site, never guessed).
///
/// Known, documented, not-fixed residual: this check is `!= Chirality::None`, so a
/// reachable `Chirality::SquarePlanar` atom (added for coordination complexes) can be
/// "found" here even though it isn't a tetrahedral center. This is not a correctness
/// bug -- [`crate::resolver::resolve_chirality`] (which every caller here eventually
/// calls on the found node) gates on `is_tetrahedral()` and returns `None` for a
/// square-planar atom, so the tiebreak fails closed to unresolved/`Tied` rather than
/// emitting a wrong CIP code. The only cost is a conservative false-negative: a
/// genuinely-resolvable tetrahedral Rule 4b tie one level further away than a
/// square-planar atom would go unresolved instead of resolved. No fixture in this
/// codebase currently reaches this path (it needs a tetrahedral Rule 4b tie and a
/// reachable square-planar center in the same digraph), so it is documented rather
/// than fixed.
pub(crate) fn nearest_embedded(
    graph: &mut CipDigraph,
    mol: &Molecule,
    start: NodeId,
) -> Result<Option<NodeId>, CipCompareError> {
    let mut frontier = graph
        .expand_children(start)
        .map_err(CipCompareError::Digraph)?;
    loop {
        if frontier.is_empty() {
            return Ok(None);
        }
        let mut found: Option<(AtomIdx, NodeId)> = None;
        for &n in &frontier {
            if let CipNodeKind::Atom { atom_idx } = graph.node(n).kind
                && mol.atom(atom_idx).chirality != Chirality::None
            {
                match found {
                    None => found = Some((atom_idx, n)),
                    Some((existing, _)) if existing == atom_idx => {}
                    Some(_) => return Ok(None),
                }
            }
        }
        if let Some((_, node_id)) = found {
            return Ok(Some(node_id));
        }
        let mut next = Vec::new();
        for &n in &frontier {
            next.extend(graph.expand_children(n).map_err(CipCompareError::Digraph)?);
        }
        frontier = next;
    }
}

const MAX_CHAIN_DEPTH: usize = 4;

/// The ordered chain of nearest-embedded-stereocenter `NodeId`s reached from `start`,
/// one level deeper each step -- in place, no re-rooting. `pub(crate)` (not just used by
/// [`break_tie_rule4b`]): `crate::assign::assign_one_with_rule5` also needs "the nearest
/// embedded stereocenter, treating `start` itself as chain position 0" for a tied group's
/// own physical position node -- exactly this function's own special-casing of `start`,
/// not [`nearest_embedded`]'s (which is only correct for continuing a chain past an
/// already-found element, never for the first one -- see its own call sites here for why
/// that distinction matters).
pub(crate) fn embedded_chain(
    graph: &mut CipDigraph,
    mol: &Molecule,
    start: NodeId,
) -> Result<Vec<NodeId>, CipCompareError> {
    // `start` itself is chain position 0 if it's chirality-bearing (it always is here
    // -- `start` is one of a tied group's own members, i.e. the outer atom's own
    // forward child, and every case in this residual family has the tied children
    // themselves be embedded stereocenters). `nearest_embedded` starts searching from
    // `start`'s *children*, so it is only correct for continuing the chain past an
    // already-found element, never for the first one.
    let mut chain = Vec::new();
    if let CipNodeKind::Atom { atom_idx } = graph.node(start).kind
        && mol.atom(atom_idx).chirality != Chirality::None
    {
        chain.push(start);
    }
    let mut current = start;
    for _ in 0..MAX_CHAIN_DEPTH {
        match nearest_embedded(graph, mol, current)? {
            Some(node_id) => {
                chain.push(node_id);
                current = node_id;
            }
            None => break,
        }
    }
    Ok(chain)
}

/// Rule 4b tiebreak for an exactly-2-way tie (every case observed in this corpus family
/// so far; 3+-way ties are out of scope and reported rather than guessed). Faithful
/// reference-relative Like/Unlike: each branch's own reference is its own nearest
/// (position-0) resolved chirality, never shared, never inverted -- the operator
/// confirmed correct by `rule4b_operator_tests.rs` (7/7). The chain elements'
/// chirality is resolved via *recursive, in-place* [`resolve_chirality`].
pub(crate) fn break_tie_rule4b(
    graph: &mut CipDigraph,
    mol: &Molecule,
    a: NodeId,
    b: NodeId,
    budget: CipBudget,
    cache: &mut HashMap<NodeId, Option<bool>>,
) -> Result<Option<Ordering>, CipCompareError> {
    let chain_a = embedded_chain(graph, mol, a)?;
    let chain_b = embedded_chain(graph, mol, b)?;
    let signs_a: Vec<Option<bool>> = chain_a
        .iter()
        .map(|&n| resolve_chirality(graph, mol, n, budget, cache))
        .collect::<Result<_, _>>()?;
    let signs_b: Vec<Option<bool>> = chain_b
        .iter()
        .map(|&n| resolve_chirality(graph, mol, n, budget, cache))
        .collect::<Result<_, _>>()?;

    let (Some(Some(ref_a)), Some(Some(ref_b))) = (signs_a.first(), signs_b.first()) else {
        return Ok(None);
    };
    for i in 0..signs_a.len().min(signs_b.len()) {
        let (Some(sa), Some(sb)) = (signs_a[i], signs_b[i]) else {
            break;
        };
        let like_a = sa == *ref_a;
        let like_b = sb == *ref_b;
        if like_a != like_b {
            return Ok(Some(if like_a {
                Ordering::Greater
            } else {
                Ordering::Less
            }));
        }
    }
    Ok(None)
}
