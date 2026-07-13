//! A minimal, tetrahedral-only R/S assignment, built on the new comparator, for the
//! sole purpose of producing labels the Milestone 2 corpus report can diff against
//! RDKit's oracle. **Not** wired into `chematic_chem::assign_cip()` or any public
//! surface beyond this crate -- explicitly experimental. No E/Z, no allene, no Rule 5 /
//! pseudoasymmetry (Milestone 4's scope); this module only touches atoms with a
//! `Chirality` annotation and exactly 4 resolvable substituent positions.
//!
//! # Positions come from `stereo_neighbor_order`, ranks come from the new comparator
//!
//! [`crate::digraph::CipDigraph`]'s root children are built by iterating
//! `Molecule::neighbors()` -- raw adjacency order, which reflects bond-*creation*
//! time, not SMILES textual encounter order. That is exactly the wrong order for
//! interpreting a stereocenter's `@`/`@@` marker (its meaning is defined relative to
//! encounter order) -- precisely the bug `d0e726b` fixed in the older engine
//! (`crates/chematic-chem/src/cip.rs`) by switching to
//! `Molecule::stereo_neighbor_order`. This module sources the four substituent
//! *positions* from `stereo_neighbor_order` (mapping `STEREO_H_SENTINEL` to the
//! digraph's `ImplicitHydrogen` child), and only the *priority ranking* of those
//! positions from [`crate::compare::rank_children`] -- reusing the exact swap-counting
//! parity algorithm already correct in `assign_tetrahedral` (mirrored below, not
//! redesigned), so this module doesn't reintroduce the order bug on its first day.

use std::collections::HashMap;

use chematic_core::{AtomIdx, Chirality, CipCode, Molecule, STEREO_H_SENTINEL};

use crate::CipError;
use crate::budget::CipBudget;
use crate::compare::{CipCompareError, CompareContext, rank_children};
use crate::digraph::CipDigraph;
use crate::node::{CipNodeKind, NodeId};

/// Why a candidate stereocenter got no assignment -- distinct from "assigned but
/// mismatched," so a caller (the corpus report, in particular) can tell "we don't know"
/// from "we got it wrong."
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkipReason {
    /// Not exactly 4 resolvable substituent positions (not a plain tetrahedral center,
    /// or `stereo_neighbor_order` unavailable).
    NotFourSubstituents,
    /// Two or more substituents are mutually tied under Rules 1a/1b/2 alone -- a
    /// genuine tie CIP can't resolve without Rule 3+ (out of scope this milestone), not
    /// a guess.
    Tied,
    /// The underlying digraph or comparator exceeded its budget for this atom.
    BudgetExceeded,
}

/// Result of the experimental tetrahedral-only assignment pass.
#[derive(Debug, Clone, Default)]
pub struct AccurateCipAssignment {
    pub assignments: Vec<(AtomIdx, CipCode)>,
    pub skipped: Vec<(AtomIdx, SkipReason)>,
}

/// Assign R/S to every tetrahedral stereocenter in `mol` that Rules 1a/1b/2 alone can
/// resolve. See module docs for scope and the positions-vs-ranks distinction.
pub fn assign_cip_accurate_experimental(
    mol: &Molecule,
    budget: CipBudget,
) -> Result<AccurateCipAssignment, CipCompareError> {
    let mut result = AccurateCipAssignment::default();

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        if atom.chirality == Chirality::None {
            continue;
        }

        let Some(stereo_order) = mol.stereo_neighbor_order(idx) else {
            result.skipped.push((idx, SkipReason::NotFourSubstituents));
            continue;
        };
        if stereo_order.len() != 4 {
            result.skipped.push((idx, SkipReason::NotFourSubstituents));
            continue;
        }

        match assign_one(mol, idx, atom.chirality, stereo_order, budget) {
            Ok(Some(code)) => result.assignments.push((idx, code)),
            Ok(None) => result.skipped.push((idx, SkipReason::NotFourSubstituents)),
            Err(SkipReason::Tied) => result.skipped.push((idx, SkipReason::Tied)),
            Err(SkipReason::BudgetExceeded) => {
                result.skipped.push((idx, SkipReason::BudgetExceeded))
            }
            Err(other) => result.skipped.push((idx, other)),
        }
    }

    Ok(result)
}

fn assign_one(
    mol: &Molecule,
    idx: AtomIdx,
    chirality: Chirality,
    stereo_order: &[u32],
    budget: CipBudget,
) -> Result<Option<CipCode>, SkipReason> {
    let mut graph = CipDigraph::new(mol, idx, budget).map_err(map_digraph_err)?;
    let root = graph.root();
    let root_children = graph.expand_children(root).map_err(map_digraph_err)?;
    if root_children.len() != 4 {
        return Ok(None);
    }

    let Some(position_nodes) = position_node_ids(&graph, &root_children, stereo_order) else {
        return Ok(None);
    };

    let mut ctx = CompareContext::new();
    let groups = rank_children(&mut graph, &root_children, &mut ctx).map_err(map_compare_err)?;
    if groups.iter().any(|g| g.len() > 1) {
        return Err(SkipReason::Tied);
    }

    // Highest-priority group first (rank_children's own convention) -> rank N down to
    // rank 1, matching assign_tetrahedral's swap-counting convention below.
    let n = groups.len() as u8;
    let mut rank_of: HashMap<NodeId, u8> = HashMap::new();
    for (group_idx, group) in groups.iter().enumerate() {
        rank_of.insert(group[0], n - group_idx as u8);
    }

    let ranks: Vec<u8> = position_nodes.iter().map(|node| rank_of[node]).collect();

    // Mirrors crates/chematic-chem/src/cip.rs::assign_tetrahedral's parity computation
    // verbatim (already correct there, fixed in d0e726b) -- only the source of
    // `ranks` differs (the new recursive comparator, not the old shell-pooling one).
    let Some(lowest_pos) = ranks.iter().position(|&r| r == 1) else {
        return Ok(None);
    };
    let parity_odd = lowest_pos % 2 == 1;
    let smiles_cw = chirality == Chirality::Clockwise;
    let cw_from_lowest = smiles_cw ^ parity_odd;

    let remaining_ranks: Vec<u8> = (0..4usize)
        .filter(|&i| i != lowest_pos)
        .map(|i| ranks[i])
        .collect();
    let Some(remaining_swaps_odd) = swap_parity(&remaining_ranks) else {
        return Ok(None);
    };

    let is_r = cw_from_lowest ^ remaining_swaps_odd;
    Ok(Some(if is_r { CipCode::R } else { CipCode::S }))
}

fn map_digraph_err(e: CipError) -> SkipReason {
    let CipError::BudgetExceeded { .. } = e;
    SkipReason::BudgetExceeded
}

fn map_compare_err(e: CipCompareError) -> SkipReason {
    match e {
        CipCompareError::BudgetExceeded { .. } | CipCompareError::Digraph(_) => {
            SkipReason::BudgetExceeded
        }
        CipCompareError::InvalidDigraph(_) => SkipReason::NotFourSubstituents,
    }
}

/// Map each `stereo_neighbor_order` position to the digraph node representing it. A
/// tetrahedral stereocenter's own substituents are always single-bonded (a `@`/`@@`
/// marker only appears on genuinely tetrahedral centers), so the root's direct
/// children are always `Atom`/`ImplicitHydrogen` kinds here -- never a duplicate.
fn position_node_ids(
    graph: &CipDigraph,
    root_children: &[NodeId],
    stereo_order: &[u32],
) -> Option<Vec<NodeId>> {
    let mut result = Vec::with_capacity(stereo_order.len());
    for &pos_val in stereo_order {
        let node_id = if pos_val == STEREO_H_SENTINEL {
            root_children
                .iter()
                .copied()
                .find(|&id| matches!(graph.node(id).kind, CipNodeKind::ImplicitHydrogen))?
        } else {
            let atom_idx = AtomIdx(pos_val);
            root_children.iter().copied().find(|&id| {
                matches!(graph.node(id).kind, CipNodeKind::Atom { atom_idx: a } if a == atom_idx)
            })?
        };
        result.push(node_id);
    }
    Some(result)
}

/// Count swaps needed to bring `remaining_ranks` (3 elements, each in `{2,3,4}`) into
/// ascending order `[2,3,4]`. Identical to `assign_tetrahedral`'s own helper.
fn swap_parity(remaining_ranks: &[u8]) -> Option<bool> {
    let mut r = remaining_ranks.to_vec();
    let target = [2u8, 3, 4];
    let mut swaps = 0usize;
    for i in 0..3 {
        if r[i] != target[i] {
            let j_rel = r[i + 1..].iter().position(|&x| x == target[i])?;
            r.swap(i, j_rel + i + 1);
            swaps += 1;
        }
    }
    Some(swaps % 2 == 1)
}
