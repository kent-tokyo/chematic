//! A minimal, tetrahedral-only R/S assignment, built on the new comparator, for the
//! sole purpose of producing labels the Milestone 2 corpus report can diff against
//! RDKit's oracle. **Not** wired into `chematic_chem::assign_cip()` or any public
//! surface beyond this crate -- explicitly experimental. No E/Z, no allene; this module
//! only touches atoms with a `Chirality` annotation and exactly 4 resolvable substituent
//! positions.
//!
//! # Rule 5 (pseudoasymmetry) -- Milestone 4A, deliberately narrow scope
//!
//! [`assign_cip_accurate_experimental`] runs a second pass ([`apply_rule5_pass`]) after
//! the Rules-1a/1b/2 pass ([`assign_all`]) completes, refining only atoms the first pass
//! left [`SkipReason::Tied`]. It resolves exactly one shape: a tie between precisely 2 of
//! the 4 physical positions, where each tied branch's *nearest* embedded, already-
//! Pass-1-resolved stereocenter (unambiguous: only one atom at that minimum depth) has a
//! code that differs from the other branch's (one R, one S) -- CIP's textbook
//! pseudoasymmetric-center pattern. "Nearest," not "only one in the whole subtree": in a
//! monocyclic ring, both ring-direction branches from the pseudoasymmetric center
//! eventually wrap around and reach *every* embedded stereocenter on the ring, just in
//! opposite order, so distinguishing "which one is closer per branch" is required, not
//! "does this branch contain only one at all" (verified empirically on both target
//! rows -- an earlier "exactly one in the whole subtree" version of this check
//! wrongly disqualified both). The resolved atom is labeled
//! [`chematic_core::CipCode::LowerR`]/[`chematic_core::CipCode::LowerS`], not `R`/`S`.
//!
//! **What this deliberately does not attempt**: a three-armed, locally-symmetric cage
//! family (verified present in the validation corpus -- flipping one arm's stereo tag
//! reclassifies all three embedded centers at once) has no seed for this pairwise
//! provisional-map approach to refine from -- every relevant neighbor is *also* tied in
//! Pass 1. That family needs symmetry/automorphism-aware joint resolution, a different
//! architecture, and is out of scope here (tracked as Milestone 4A-2 in
//! `docs/cip_accurate_rfc.md`). [`apply_rule5_pass`] is a structural no-op on it: the tie
//! detection requires *exactly* 2 physical positions in *exactly* one tied group, which
//! the cage family's 3-way (or worse) ties never satisfy.
//!
//! **Auxiliary vs. molecular descriptor**: real CIP Rule 4c/5 compares an embedded
//! center's *auxiliary* descriptor -- computed within the digraph rooted at the outer
//! stereocenter -- not necessarily its own independently-computed molecular R/S. This
//! pass uses the molecular descriptor (`provisional`, built from Pass 1's whole-molecule
//! results) as a stand-in. That is verified, not assumed, correct for this milestone's
//! 2-row target: both target molecules' embedded reference centers are already
//! independently and uniquely resolved by Pass 1 with no ring-duplicate/phantom-atom
//! complication in their own ranking. It is not a general auxiliary-descriptor
//! implementation and should not be trusted as one outside this scope.
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
    let pass1 = assign_all(mol, budget, kekule.as_ref())?;
    Ok(apply_rule5_pass(mol, budget, kekule.as_ref(), pass1))
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

/// Milestone 4A's Rule 5 refinement -- see module docs for scope. Only ever touches
/// atoms Pass 1 (`assign_all`, above) left [`SkipReason::Tied`]; every other atom
/// (resolved or skipped for any other reason) is carried through unchanged. Not called
/// by [`assign_cip_accurate_experimental_without_mancude`] -- that function's whole
/// purpose is to stay a frozen, Rule-5-independent reference point (see its own doc
/// comment), so Rule 5 is deliberately only wired into the main, live entry point.
fn apply_rule5_pass(
    mol: &Molecule,
    budget: CipBudget,
    kekule: Option<&(Molecule, MancudeContext)>,
    pass1: AccurateCipAssignment,
) -> AccurateCipAssignment {
    let provisional: HashMap<AtomIdx, CipCode> = pass1.assignments.iter().copied().collect();

    let mut assignments = pass1.assignments;
    let mut skipped = Vec::with_capacity(pass1.skipped.len());

    for (idx, reason) in pass1.skipped {
        if reason != SkipReason::Tied {
            skipped.push((idx, reason));
            continue;
        }

        let atom = mol.atom(idx);
        let Some(stereo_order) = mol.stereo_neighbor_order(idx) else {
            skipped.push((idx, reason));
            continue;
        };
        if stereo_order.len() != 4 {
            skipped.push((idx, reason));
            continue;
        }

        match assign_one_with_rule5(
            mol,
            idx,
            atom.chirality,
            stereo_order,
            budget,
            kekule,
            &provisional,
        ) {
            Ok(Some(code)) => assignments.push((idx, code)),
            _ => skipped.push((idx, reason)),
        }
    }

    AccurateCipAssignment {
        assignments,
        skipped,
    }
}

/// Retry a Pass-1-tied atom with Rule 5. Identical digraph/ranking setup to
/// [`assign_one`]; the only new logic is locating a single, exactly-2-physical-position
/// tied group and, if each side embeds exactly one already-resolved stereocenter with a
/// differing `provisional` code, breaking the tie by that code (R precedes S) before
/// handing the (now fully split) group partition to the same
/// [`resolve_is_r_from_groups`] parity math Pass 1 uses. Any shape this narrow detector
/// doesn't recognize -- 0 or 2+ tied groups, a tied group of size != 2, an
/// unresolved/ambiguous embedded center, or matching codes on both sides -- returns
/// `Err(SkipReason::Tied)` unchanged, exactly Pass 1's own outcome.
fn assign_one_with_rule5(
    mol: &Molecule,
    idx: AtomIdx,
    chirality: Chirality,
    stereo_order: &[u32],
    budget: CipBudget,
    kekule: Option<&(Molecule, MancudeContext)>,
    provisional: &HashMap<AtomIdx, CipCode>,
) -> Result<Option<CipCode>, SkipReason> {
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

    let position_set: HashSet<NodeId> = position_nodes.iter().copied().collect();
    let mut tied_group_idx = None;
    for (gi, group) in groups.iter().enumerate() {
        let physical_count = group.iter().filter(|n| position_set.contains(n)).count();
        if physical_count > 1 {
            // A 3+-way tie in one group, or a second independently-tied group, is
            // outside this milestone's narrow scope -- stay tied rather than guess.
            if physical_count != 2 || tied_group_idx.is_some() {
                return Err(SkipReason::Tied);
            }
            tied_group_idx = Some(gi);
        }
    }
    let Some(tied_group_idx) = tied_group_idx else {
        // Pass 2 only ever runs on atoms Pass 1 already found tied, so this shouldn't
        // happen -- stay conservative rather than assume Pass 1/Pass 2 agree.
        return Err(SkipReason::Tied);
    };

    let tied_group = &groups[tied_group_idx];
    let physical_in_tied: Vec<NodeId> = tied_group
        .iter()
        .copied()
        .filter(|n| position_set.contains(n))
        .collect();
    let (pos_a, pos_b) = (physical_in_tied[0], physical_in_tied[1]);

    let Some(atom_a) =
        nearest_embedded_stereocenter(&mut graph, mol, pos_a).map_err(map_digraph_err)?
    else {
        return Err(SkipReason::Tied);
    };
    let Some(atom_b) =
        nearest_embedded_stereocenter(&mut graph, mol, pos_b).map_err(map_digraph_err)?
    else {
        return Err(SkipReason::Tied);
    };
    let (Some(&code_a), Some(&code_b)) = (provisional.get(&atom_a), provisional.get(&atom_b))
    else {
        return Err(SkipReason::Tied);
    };
    if code_a == code_b {
        // Both R or both S -- not a distinguishing pseudoasymmetric pair in this scope;
        // stay tied rather than guess.
        return Err(SkipReason::Tied);
    }

    // Rule 5: R precedes S. Verified against this milestone's own 2-row corpus target
    // (both currently-known cases resolve to lowercase `r`), not assumed from the
    // textbook statement alone -- see module docs.
    let (higher, lower) = if code_a == CipCode::R {
        (pos_a, pos_b)
    } else {
        (pos_b, pos_a)
    };

    // Split the one tied group into two singleton-physical-position groups (higher
    // priority first), leaving every other group's relative order untouched. Any
    // non-physical (duplicate) sibling that was in the tied group rides along with
    // `lower` -- harmless, since `resolve_is_r_from_groups` only ever looks up ranks
    // for the 4 physical `position_nodes`, never for duplicates directly.
    let mut resolved_groups: Vec<Vec<NodeId>> = Vec::with_capacity(groups.len() + 1);
    for (gi, group) in groups.iter().enumerate() {
        if gi == tied_group_idx {
            resolved_groups.push(vec![higher]);
            let mut lower_group = vec![lower];
            lower_group.extend(group.iter().copied().filter(|&n| n != higher && n != lower));
            resolved_groups.push(lower_group);
        } else {
            resolved_groups.push(group.clone());
        }
    }

    let Some(is_r) = resolve_is_r_from_groups(&resolved_groups, &position_nodes, chirality) else {
        return Ok(None);
    };
    Ok(Some(if is_r {
        CipCode::LowerR
    } else {
        CipCode::LowerS
    }))
}

/// Breadth-first search from `node`, level by level, for the *nearest* atom with its own
/// `Chirality` set (checked on `mol`, the original un-kekulized molecule -- chirality is
/// atom-level SMILES stereo info, unaffected by kekulization). Returns `None` if the
/// nearest level containing a chirality-bearing atom contains more than one (ambiguous --
/// which one is "nearest" is undefined), or if no such atom exists anywhere in the
/// subtree.
///
/// Deliberately "nearest," not "the only one in the whole subtree": in a monocyclic
/// ring, walking either ring-direction branch from the pseudoasymmetric center
/// eventually wraps all the way around and reaches *every* embedded stereocenter on the
/// ring (just in opposite order per branch) -- so a whole-subtree "exactly one" count
/// would find 2+ in both branches and always disqualify, which is wrong (verified
/// against this milestone's own 2-row target, which are both single-ring cases). The
/// three-armed cage family this milestone excludes (see module docs) instead fails via
/// *ambiguity at the nearest level* (that family's branches reach 2+ equally-near
/// embedded stereocenters), which still safely falls through to `SkipReason::Tied`
/// rather than producing a wrong lowercase label.
fn nearest_embedded_stereocenter(
    graph: &mut CipDigraph,
    mol: &Molecule,
    node: NodeId,
) -> Result<Option<AtomIdx>, CipError> {
    let mut frontier = vec![node];
    while !frontier.is_empty() {
        let mut found_this_level: Option<AtomIdx> = None;
        let mut next_frontier = Vec::new();
        for &current in &frontier {
            if let CipNodeKind::Atom { atom_idx } = graph.node(current).kind
                && mol.atom(atom_idx).chirality != Chirality::None
            {
                match found_this_level {
                    None => found_this_level = Some(atom_idx),
                    Some(existing) if existing == atom_idx => {}
                    Some(_) => return Ok(None),
                }
            }
            next_frontier.extend(graph.expand_children(current)?);
        }
        if let Some(atom_idx) = found_this_level {
            return Ok(Some(atom_idx));
        }
        frontier = next_frontier;
    }
    Ok(None)
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

    let Some(is_r) = resolve_is_r_from_groups(&groups, &position_nodes, chirality) else {
        return Ok(None);
    };
    Ok(Some(if is_r { CipCode::R } else { CipCode::S }))
}

/// Given a fully-ranked group partition (highest-priority group first, matching
/// [`rank_children`]'s convention) and the 4 physical positions in `stereo_order`
/// order, compute whether the center is *rectus* (R). Shared by [`assign_one`] (Rules
/// 1a/1b/2 only) and [`assign_one_with_rule5`] (which passes a `groups` partition with
/// one tied group already split by a Rule 5 tiebreak) -- the swap-counting parity math
/// itself doesn't know or care which rule produced the ordering.
fn resolve_is_r_from_groups(
    groups: &[Vec<NodeId>],
    position_nodes: &[NodeId],
    chirality: Chirality,
) -> Option<bool> {
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
    let lowest_pos = ranks.iter().position(|&r| r == 1)?;
    let parity_odd = lowest_pos % 2 == 1;
    let smiles_cw = chirality == Chirality::Clockwise;
    let cw_from_lowest = smiles_cw ^ parity_odd;

    let remaining_ranks: Vec<u8> = (0..4usize)
        .filter(|&i| i != lowest_pos)
        .map(|i| ranks[i])
        .collect();
    let remaining_swaps_odd = swap_parity(&remaining_ranks)?;

    Some(cw_from_lowest ^ remaining_swaps_odd)
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
