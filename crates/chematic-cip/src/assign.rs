//! A minimal, tetrahedral-only R/S assignment, built on the new comparator, for the
//! sole purpose of producing labels the Milestone 2 corpus report can diff against
//! RDKit's oracle. **Not** wired into `chematic_chem::assign_cip()` or any public
//! surface beyond this crate -- explicitly experimental. No E/Z, no allene; this module
//! only touches atoms with a `Chirality` annotation and exactly 4 resolvable substituent
//! positions.
//!
//! # Rule 5 (pseudoasymmetry) -- Milestone 4A (2-row scope), generalized in
//! Milestone 4A-2 to a real bottom-up auxiliary descriptor
//!
//! [`assign_cip_accurate_experimental`] runs a pass ([`apply_rule5_pass`]) after the
//! Rules-1a/1b/2 pass ([`assign_all`]) and the Rule 4b pass (`crate::resolver`) complete,
//! refining only atoms still left [`SkipReason::Tied`]. It resolves one shape: a tie
//! between precisely 2 of the 4 physical positions, where each tied branch's *nearest*
//! embedded stereocenter has an **auxiliary** R/S sign (see below) that differs from the
//! other branch's (one R, one S) -- CIP's textbook pseudoasymmetric-center pattern. The
//! resolved atom is labeled [`chematic_core::CipCode::LowerR`]/
//! [`chematic_core::CipCode::LowerS`], not `R`/`S`.
//!
//! **Auxiliary, not molecular, descriptor** (Milestone 4A-2's fix to Milestone 4A's own
//! documented limitation). Real CIP Rule 4c/5 compares an embedded center's *auxiliary*
//! descriptor -- computed bottom-up, within the very same digraph rooted at the outer
//! stereocenter under examination, per Hanson, Musacchio, Mayfield et al. 2018 (*J. Chem.
//! Inf. Model.* 58(9), 1755-1765): "the descriptor for an auxiliary center does not
//! depend upon any center between it and the root... the priority of a ligand leading
//! back to the digraph root will always be ranked by Rule 1a." Milestone 4A shipped a
//! *molecular*-descriptor stand-in instead (a `provisional: HashMap<AtomIdx, CipCode>`
//! built from Pass 1's whole-molecule results), verified correct only for its own 2-row
//! target where the embedded reference happened to already be independently resolved.
//! That stand-in structurally cannot handle an embedded reference that is *itself*
//! Pass-1/Rule-4b tied -- exactly the three-armed, locally-symmetric adamantane-cage
//! family (`[C]12C[CH]3C[CH](C[CH](C3)C1)C2`-shaped, verified present in the validation
//! corpus, 15 rows across 5 distinct molecules) Milestone 4A deferred as "Milestone
//! 4A-2, needs symmetry/automorphism-aware joint resolution, a different architecture."
//!
//! That framing turned out to overstate what was needed. Diagnosing the cage family
//! directly (see `docs/rfcs/cip_accurate_rfc.md`'s Milestone 4A-2 entry) found **no genuine
//! cross-atom cycle**: within one digraph rooted at the outer atom, every embedded
//! reference is strictly deeper (Hanson's own bottom-up postulate above), so
//! `crate::resolver::resolve_chirality` -- already built, and already validated 72/72 for
//! Rule 4b -- computes a correct-shaped `Option<bool>` auxiliary sign for the tied
//! branches' nearest embedded stereocenters directly, with no modification, on every one
//! of the 15 rows. [`assign_one_with_rule5`] now calls `resolve_chirality`
//! (`crate::rule4b::nearest_embedded` to locate the reference node, matching Rule 4b's own
//! mechanism) instead of looking `provisional` up, so an embedded reference that is
//! itself tied gets resolved in place, recursively, rather than looked up in a map that
//! was never populated for it. No SCC/fixed-point solver was built or needed -- see the
//! RFC entry for the falsified alternative hypothesis and the concrete evidence.
//!
//! **Still out of scope**: a tied pair whose *nearest* embedded references share the same
//! auxiliary sign (both R or both S -- not a distinguishing mirror-image pair) does not
//! attempt a deeper chain comparison the way Rule 4b's `break_tie_rule4b` does; no row in
//! this project's validation corpus exercises that shape for Rule 5, so extending to it
//! is deferred rather than guessed. See `docs/rfcs/cip_accurate_rfc.md` for the full writeup.
//!
//! **Element-level guard: phosphorus stays tied.** The same code-path fix above, as a
//! side effect, also reaches 2 cyclophosphazene phosphorus stereocenters
//! (`docs/rfcs/cip_accurate_rfc.md` Milestone 4C-1) that were previously `SkipReason::Tied`
//! for the identical chain-length-1 Rule 4b degeneracy the carbon cage family has.
//! Unlike the 15 carbon rows, Milestone 4C-1 found **neither** RDKit CIP engine has a
//! representation-stable answer for that phosphorus molecule -- both flip under a
//! chemically-neutral Kekule respelling -- so there is no oracle a resolved phosphorus
//! label could be checked against. [`assign_one_with_rule5`] therefore never emits a
//! resolved label for a **phosphorus** stereocenter; the element is checked once,
//! cheaply, before any digraph work, and falls back to `Err(SkipReason::Tied)` -- exactly
//! this function's own existing convention for every other shape it doesn't recognize
//! (see below). This is an element-level guard, not a molecule-specific one: it excludes
//! every phosphorus stereocenter that reaches this path, not one specific SMILES. Note
//! this project has, as of this writing, zero verified examples of this path being
//! correct for *any* non-carbon element -- phosphorus is simply the only non-carbon
//! element the validation corpus happens to exercise here, so a broader "unverified for
//! any non-carbon element" framing may be more honest than "unverified for phosphorus
//! specifically"; that broader guard is *not* implemented here, only flagged (see
//! `docs/rfcs/cip_accurate_rfc.md`'s Milestone 4C-1 entry) -- narrowing to exactly the
//! element asked about avoids expanding scope
//! unilaterally.
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

use chematic_core::{AtomIdx, Chirality, CipCode, Element, Molecule, STEREO_H_SENTINEL};

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
    // Rule 4b (Milestone 4B-2) must run before Rule 5 -- CIP rule order. Rule 5's own
    // `provisional` map (built inside `apply_rule5_pass` from whatever
    // `AccurateCipAssignment` it's handed) sees Rule 4b's newly-resolved atoms for
    // free; `apply_rule5_pass` needs no changes for this. See `crate::resolver`.
    let pass2 = crate::resolver::apply_rule4b_pass(mol, budget, kekule.as_ref(), pass1);
    Ok(apply_rule5_pass(mol, budget, kekule.as_ref(), pass2))
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
        if !atom.chirality.is_tetrahedral() {
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

/// Rule 5 refinement -- see module docs for scope. Only ever touches atoms Pass 1
/// (`assign_all`, above) and the Rule 4b pass (`crate::resolver::apply_rule4b_pass`) left
/// [`SkipReason::Tied`]; every other atom (resolved or skipped for any other reason) is
/// carried through unchanged. Not called by
/// [`assign_cip_accurate_experimental_without_mancude`] -- that function's whole purpose
/// is to stay a frozen, Rule-5-independent reference point (see its own doc comment), so
/// Rule 5 is deliberately only wired into the main, live entry point.
fn apply_rule5_pass(
    mol: &Molecule,
    budget: CipBudget,
    kekule: Option<&(Molecule, MancudeContext)>,
    pass1: AccurateCipAssignment,
) -> AccurateCipAssignment {
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

        match assign_one_with_rule5(mol, idx, atom.chirality, stereo_order, budget, kekule) {
            Ok(Some(code)) => assignments.push((idx, code)),
            _ => skipped.push((idx, reason)),
        }
    }

    AccurateCipAssignment {
        assignments,
        skipped,
    }
}

/// Retry a Pass-1/Rule-4b-tied atom with Rule 5. Identical digraph/ranking setup to
/// [`assign_one`]; the new logic is locating a single, exactly-2-physical-position tied
/// group and, if each side's *nearest* embedded stereocenter resolves (via
/// `crate::resolver::resolve_chirality`, the same bottom-up auxiliary-descriptor
/// mechanism Rule 4b validated 72/72) to a differing auxiliary R/S sign, breaking the tie
/// by that sign (R precedes S) before handing the (now fully split) group partition to
/// the same [`resolve_is_r_from_groups`] parity math Pass 1 uses. Any shape this detector
/// doesn't recognize -- 0 or 2+ tied groups, a tied group of size != 2, an
/// unresolved/ambiguous embedded reference (`crate::rule4b::nearest_embedded` returning
/// `None`), or matching auxiliary signs on both sides -- returns `Err(SkipReason::Tied)`
/// unchanged, exactly Pass 1's own outcome (never a guess).
fn assign_one_with_rule5(
    mol: &Molecule,
    idx: AtomIdx,
    chirality: Chirality,
    stereo_order: &[u32],
    budget: CipBudget,
    kekule: Option<&(Molecule, MancudeContext)>,
) -> Result<Option<CipCode>, SkipReason> {
    // Element-level guard (see module docs, "Element-level guard: phosphorus stays
    // tied"): this auxiliary-resolution path is oracle-verified only for carbon
    // stereocenters (all 15 corpus rows it resolves are carbon). The only other element
    // it has ever been observed to reach is phosphorus, on a cyclophosphazene molecule
    // whose RDKit oracle is itself representation-unstable (Milestone 4C-1) -- there is
    // no stable answer to check a resolved phosphorus label against, so never emit one.
    // `idx`/`mol` name the same atom identity regardless of Kekule respelling (`mol` is
    // always the original, pre-Kekule molecule here -- see `assign_one`'s own doc
    // comment), so this check is stable across resonance respellings by construction.
    if mol.atom(idx).element == Element::P {
        return Err(SkipReason::Tied);
    }

    let mut graph = match kekule {
        Some((kekule_mol, ctx)) => {
            CipDigraph::new_with_mancude(kekule_mol, idx, budget, ctx).map_err(map_digraph_err)?
        }
        None => CipDigraph::new(mol, idx, budget).map_err(map_digraph_err)?,
    };
    let root = graph.root();
    let root_children = graph.expand_children(root).map_err(map_digraph_err)?;

    let Some(position_nodes) = position_node_ids(&graph, stereo_order, &root_children) else {
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

    // Auxiliary descriptor, not molecular: each branch's own nearest embedded
    // stereocenter's R/S sign, computed bottom-up *within this same digraph*
    // (`crate::resolver::resolve_chirality`) -- not a lookup into a whole-molecule
    // "already resolved" map, which cannot succeed when the embedded atom is itself
    // Pass-1/Rule-4b tied (exactly the three-armed cage family's shape; see module
    // docs). `embedded_chain` (not `nearest_embedded` directly -- see that function's
    // own doc comment) is the same reference-location mechanism Rule 4b validated
    // 72/72, reused rather than reimplemented: `pos_a`/`pos_b` are the tied group's own
    // physical position nodes, so `pos_a`/`pos_b` themselves are chain position 0 when
    // they're directly chirality-bearing (a monocyclic ring's immediate ring
    // neighbor), and `embedded_chain` searches onward only when they aren't (a plain
    // CH2 bridge, e.g. the three-armed cage family). Only chain position 0 (the
    // nearest reference) is used here -- see module docs for why a deeper chain
    // comparison is out of scope.
    let chain_a = crate::rule4b::embedded_chain(&mut graph, mol, pos_a).map_err(map_compare_err)?;
    let chain_b = crate::rule4b::embedded_chain(&mut graph, mol, pos_b).map_err(map_compare_err)?;
    let (Some(&embedded_a), Some(&embedded_b)) = (chain_a.first(), chain_b.first()) else {
        return Err(SkipReason::Tied);
    };

    let mut cache = HashMap::new();
    let sign_a =
        crate::resolver::resolve_chirality(&mut graph, mol, embedded_a, budget, &mut cache)
            .map_err(map_compare_err)?;
    let sign_b =
        crate::resolver::resolve_chirality(&mut graph, mol, embedded_b, budget, &mut cache)
            .map_err(map_compare_err)?;
    let (Some(is_r_a), Some(is_r_b)) = (sign_a, sign_b) else {
        // The embedded reference itself has no resolvable auxiliary sign (a deeper,
        // still-unresolvable tie) -- stay tied rather than guess.
        return Err(SkipReason::Tied);
    };
    if is_r_a == is_r_b {
        // Both R or both S -- not a distinguishing mirror-image pair in this scope
        // (would need a deeper chain comparison, unexercised by this project's
        // validation corpus for Rule 5 -- see module docs); stay tied rather than guess.
        return Err(SkipReason::Tied);
    }

    // Rule 5: R precedes S. Verified against this project's full 17-row pseudoasymmetric
    // corpus (2 Milestone-4A rows + 15 Milestone-4A-2 cage-family rows, all lowercase
    // `r`/`s` per the RDKit oracle) -- not assumed from the textbook statement alone.
    let (higher, lower) = if is_r_a {
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

    let Some(position_nodes) = position_node_ids(&graph, stereo_order, &root_children) else {
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
/// 1a/1b/2 only), [`assign_one_with_rule5`] (which passes a `groups` partition with
/// one tied group already split by a Rule 5 tiebreak), and (Milestone 4B-2)
/// `crate::resolver::assign_one_with_rule4b`/`resolve_chirality` -- the swap-counting
/// parity math itself doesn't know or care which rule produced the ordering.
pub(crate) fn resolve_is_r_from_groups(
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

pub(crate) fn map_digraph_err(e: CipError) -> SkipReason {
    let CipError::BudgetExceeded { .. } = e;
    SkipReason::BudgetExceeded
}

pub(crate) fn map_compare_err(e: CipCompareError) -> SkipReason {
    match e {
        CipCompareError::BudgetExceeded { .. } | CipCompareError::Digraph(_) => {
            SkipReason::BudgetExceeded
        }
        CipCompareError::InvalidDigraph(_) => SkipReason::NotFourSubstituents,
    }
}

/// Map each `stereo_neighbor_order` position to the digraph node representing it.
///
/// For a tetrahedral stereocenter's own root children (the only case before Milestone
/// 4B-2), substituents are always single-bonded, so `candidates` only ever contains
/// `Atom`/`ImplicitHydrogen` kinds -- never a duplicate. Milestone 4B-2's Rule 4b
/// resolver additionally calls this for an *embedded* stereocenter's own 4 positions
/// (its forward children plus the back-to-root parent node,
/// `crate::resolver::resolve_chirality`), where a physical position can be a
/// `RingDuplicate` (the stereo_neighbor_order-named atom is reached via a ring
/// closure at that embedded position, not as a fresh real node) -- hence the
/// additional match arm below. Provably inert for every pre-existing caller
/// ([`assign_one`]/[`assign_one_with_rule5`]), which only ever pass a root's own
/// direct children.
pub(crate) fn position_node_ids(
    graph: &CipDigraph,
    stereo_order: &[u32],
    candidates: &[NodeId],
) -> Option<Vec<NodeId>> {
    let mut result = Vec::with_capacity(stereo_order.len());
    for &pos_val in stereo_order {
        let node_id = if pos_val == STEREO_H_SENTINEL {
            candidates
                .iter()
                .copied()
                .find(|&id| matches!(graph.node(id).kind, CipNodeKind::ImplicitHydrogen))?
        } else {
            let atom_idx = AtomIdx(pos_val);
            candidates
                .iter()
                .copied()
                .find(|&id| match graph.node(id).kind {
                    CipNodeKind::Atom { atom_idx: a } => a == atom_idx,
                    CipNodeKind::RingDuplicate { closure_atom, .. } => closure_atom == atom_idx,
                    _ => false,
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
