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
//!
//! # Physical ligands vs. duplicate nodes
//!
//! The digraph root's children are not always exactly the stereocenter's 4 physical
//! neighbors: a multiple bond *at* the stereocenter itself (e.g. a P=N phosphazene
//! center) adds one or more [`CipNodeKind::MultipleBondDuplicate`] siblings alongside the
//! real neighbor -- 5 root children for one double bond, not 4. `stereo_neighbor_order`
//! only ever names the 4 real physical neighbors (never a duplicate), so
//! [`position_node_ids`] already only ever resolves to real `Atom`/`ImplicitHydrogen`
//! nodes. What must NOT happen is treating a duplicate as if it were competing for one of
//! those 4 slots: [`assign_one`] ranks the *entire* root-children set (duplicates
//! included, since a duplicate's presence is real information for ranking a real
//! neighbor's own priority), but only ever treats a tie as unresolvable when it's between
//! two of the 4 *physical* positions -- a duplicate tying with anything doesn't block
//! assignment. Ranks are then dense-remapped to `1..=4` before the swap-parity step,
//! since a duplicate can occupy a rank slot between two physical positions.

use std::collections::{HashMap, HashSet};

use chematic_core::{AtomIdx, Chirality, CipCode, Molecule, STEREO_H_SENTINEL};

use crate::CipError;
use crate::budget::CipBudget;
use crate::compare::{CipCompareError, CompareContext, rank_children};
use crate::digraph::CipDigraph;
use crate::mancude::{MancudeContext, prepare_kekule_form};
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
///
/// Computes `mol`'s Kekulé-form clone and [`MancudeContext`] **once**, before the
/// per-atom loop, and shares both across every stereocenter's digraph -- never
/// recomputed per atom or per subtree expansion (a whole-molecule quantity). If Kekulé
/// form can't be computed at all for `mol` (rare -- e.g. a non-bipartite aromatic system
/// `chematic_core::kekulization::kekulize` can't resolve), falls back to the plain,
/// pre-Milestone-3B-1 digraph path (`CipDigraph::new` on the original aromatic-notation
/// `mol`, no `MancudeContext`) for that molecule rather than failing the whole
/// assignment -- exactly today's behavior for such a molecule, since it never had a
/// MANCUDE-fractional path to lose.
pub fn assign_cip_accurate_experimental(
    mol: &Molecule,
    budget: CipBudget,
) -> Result<AccurateCipAssignment, CipCompareError> {
    let kekule = prepare_kekule_form(mol).ok();
    assign_all(mol, budget, kekule.as_ref())
}

/// Identical to [`assign_cip_accurate_experimental`], but never attaches a
/// [`MancudeContext`] -- reproduces exactly the pre-Milestone-3B-1b digraph
/// construction (plain `CipDigraph::new`, aromatic bonds contribute no
/// `MultipleBondDuplicate` nodes). Exists as a stable reference point for regression
/// tooling and tests that need to classify a stereocenter's wrong-vs-tied outcome
/// independent of whatever the live, MANCUDE-aware engine currently does -- see
/// `tests/common/mod.rs::is_bucket_misclassified`'s module docs for why that
/// independence matters (gating a *structural corpus scope* on the live engine's current
/// correctness makes the scope shrink every time the engine improves).
pub fn assign_cip_accurate_experimental_without_mancude(
    mol: &Molecule,
    budget: CipBudget,
) -> Result<AccurateCipAssignment, CipCompareError> {
    assign_all(mol, budget, None)
}

fn assign_all(
    mol: &Molecule,
    budget: CipBudget,
    kekule: Option<&(Molecule, MancudeContext)>,
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

        match assign_one(mol, idx, atom.chirality, stereo_order, budget, kekule) {
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
    kekule: Option<&(Molecule, MancudeContext)>,
) -> Result<Option<CipCode>, SkipReason> {
    // `apply_kekule` preserves `AtomIdx` values exactly (verified Milestone 3B-0), so
    // `idx`/`stereo_order` (sourced from the original `mol`, above) name the same
    // physical atoms in `kekule_mol` -- no remapping needed either way.
    let mut graph = match kekule {
        Some((kekule_mol, ctx)) => {
            CipDigraph::new_with_mancude(kekule_mol, idx, budget, ctx).map_err(map_digraph_err)?
        }
        None => CipDigraph::new(mol, idx, budget).map_err(map_digraph_err)?,
    };
    let root = graph.root();
    let root_children = graph.expand_children(root).map_err(map_digraph_err)?;

    let Some(position_nodes) = position_node_ids(&graph, &root_children, stereo_order) else {
        return Ok(None);
    };

    let mut ctx = CompareContext::new();
    let groups = rank_children(&mut graph, &root_children, &mut ctx).map_err(map_compare_err)?;

    // A tie only blocks resolution when two of the 4 *physical* positions land in the
    // same group -- a duplicate node tying with anything (another duplicate, or even a
    // physical position) doesn't compete for a stereo_neighbor_order slot. See module
    // docs ("Physical ligands vs. duplicate nodes").
    let position_set: HashSet<NodeId> = position_nodes.iter().copied().collect();
    for group in &groups {
        if group.iter().filter(|n| position_set.contains(n)).count() > 1 {
            return Err(SkipReason::Tied);
        }
    }

    // Rank every node in every group (not just each group's first member) -- a duplicate
    // can share a physical position's group, and that position's rank_of lookup below
    // must still resolve. Highest-priority group first (rank_children's own convention)
    // -> rank N down to rank 1, matching assign_tetrahedral's swap-counting convention.
    let n = groups.len() as u8;
    let mut rank_of: HashMap<NodeId, u8> = HashMap::new();
    for (group_idx, group) in groups.iter().enumerate() {
        for &node in group {
            rank_of.insert(node, n - group_idx as u8);
        }
    }

    let raw_ranks: Vec<u8> = position_nodes.iter().map(|node| rank_of[node]).collect();
    // Dense-remap to 1..=4: a duplicate sibling (e.g. from a double bond at the
    // stereocenter itself) can occupy a rank slot between two physical positions, so
    // their raw ranks aren't necessarily {1,2,3,4} contiguously.
    let mut distinct_ranks = raw_ranks.clone();
    distinct_ranks.sort_unstable();
    distinct_ranks.dedup();
    let ranks: Vec<u8> = raw_ranks
        .iter()
        .map(|&r| distinct_ranks.iter().position(|&x| x == r).unwrap() as u8 + 1)
        .collect();

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
