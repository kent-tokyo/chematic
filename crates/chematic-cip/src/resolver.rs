//! Milestone 4B-2 production port of the validated Rule 4b reference engine
//! (`examples/rule4b_bottom_up.rs`, 72/72 across 4 oracle corpora including mirrors --
//! see `docs/rfcs/cip_accurate_rfc.md`). Mechanical port: [`resolve_chirality`] is
//! unchanged from the example. [`assign_one_with_rule4b`]/[`apply_rule4b_pass`] are
//! new, fitting the validated algorithm into `assign.rs`'s existing pass-based
//! architecture (mirroring [`crate::assign`]'s own `assign_one_with_rule5`/
//! `apply_rule5_pass` plumbing shape, but the tie-detection loop itself mirrors the
//! reference engine's `resolve_outer_root` -- looping over *every* group and
//! attempting a Rule 4b tiebreak on each clean 2-way tie found, not
//! `assign_one_with_rule5`'s narrower "exactly one tied group" restriction, since
//! that restriction was never part of what was validated 72/72).

use std::cmp::Ordering;
use std::collections::HashMap;

use chematic_core::{AtomIdx, Chirality, CipCode, Molecule};

use crate::assign::{
    AccurateCipAssignment, SkipReason, map_compare_err, map_digraph_err, position_node_ids,
    resolve_is_r_from_groups,
};
use crate::budget::CipBudget;
use crate::compare::{CipCompareError, CompareContext, rank_children};
use crate::digraph::CipDigraph;
use crate::mancude::MancudeContext;
use crate::node::{CipNodeKind, NodeId};
use crate::rule4b::break_tie_rule4b;

/// The core recursive, memoized, in-place resolver. `node_id` must be a real `Atom`
/// node (not the true digraph root -- the root has no "back to root" ligand and is
/// handled by [`assign_one_with_rule4b`] directly, mirroring the reference engine's
/// own `resolve_outer_root`/`resolve_chirality` split).
pub(crate) fn resolve_chirality(
    graph: &mut CipDigraph,
    mol: &Molecule,
    node_id: NodeId,
    budget: CipBudget,
    cache: &mut HashMap<NodeId, Option<bool>>,
) -> Result<Option<bool>, CipCompareError> {
    if let Some(&cached) = cache.get(&node_id) {
        return Ok(cached);
    }

    let atom_idx = match graph.node(node_id).kind {
        CipNodeKind::Atom { atom_idx } => atom_idx,
        _ => {
            cache.insert(node_id, None);
            return Ok(None);
        }
    };
    if !mol.atom(atom_idx).chirality.is_tetrahedral() {
        cache.insert(node_id, None);
        return Ok(None);
    }
    let Some(parent_id) = graph.node(node_id).parent else {
        // True digraph root -- use assign_one_with_rule4b instead.
        cache.insert(node_id, None);
        return Ok(None);
    };

    let forward_children = graph
        .expand_children(node_id)
        .map_err(CipCompareError::Digraph)?;
    let mut ctx = CompareContext::new();
    let mut groups = rank_children(graph, &forward_children, &mut ctx)?;

    // Rule 4b tiebreak: only a clean 2-way tie is handled (scope match with every case
    // seen in this residual family). 3+-way ties are reported, not guessed.
    for gi in 0..groups.len() {
        if groups[gi].len() == 2 {
            let (a, b) = (groups[gi][0], groups[gi][1]);
            if let Some(ord) = break_tie_rule4b(graph, mol, a, b, budget, cache)? {
                let (higher, lower) = match ord {
                    Ordering::Greater => (a, b),
                    Ordering::Less => (b, a),
                    Ordering::Equal => (a, b),
                };
                if ord != Ordering::Equal {
                    groups[gi] = vec![lower];
                    groups.insert(gi, vec![higher]);
                }
            }
        } else if groups[gi].len() > 2 {
            cache.insert(node_id, None);
            return Ok(None); // out of scope this round
        }
    }

    // Insert the back-to-root ligand at its Rule-1a-only rank among `groups`. Compare
    // sphere-1-to-sphere-1: the back ligand's sphere 1 is `ancestor` itself; the
    // forward ligand's sphere 1 is `rep` itself -- NOT `rep`'s children (that would be
    // off by one sphere).
    let mut insert_at = groups.len();
    for (gi, group) in groups.iter().enumerate() {
        let rep = group[0];
        let back_frontier = vec![crate::auxiliary::BackItem::Ascending {
            ancestor: parent_id,
            came_from: node_id,
        }];
        let ord = crate::auxiliary::compare_rule1a_only(graph, back_frontier, vec![rep])?;
        if ord == Ordering::Greater {
            insert_at = gi;
            break;
        }
    }
    let mut final_groups = groups.clone();
    final_groups.insert(insert_at, vec![parent_id]);

    let stereo_order = match mol.stereo_neighbor_order(atom_idx) {
        Some(o) => o,
        None => {
            cache.insert(node_id, None);
            return Ok(None);
        }
    };
    let mut candidates = forward_children.clone();
    candidates.push(parent_id);
    let position_nodes = match position_node_ids(graph, stereo_order, &candidates) {
        Some(p) => p,
        None => {
            cache.insert(node_id, None);
            return Ok(None);
        }
    };
    let result =
        resolve_is_r_from_groups(&final_groups, &position_nodes, mol.atom(atom_idx).chirality);
    cache.insert(node_id, result);
    Ok(result)
}

/// Retry a Pass-1-tied atom with Rule 4b. Builds one fresh digraph rooted at the
/// outer atom itself (`idx`) -- matching `assign_one`/`assign_one_with_rule5`'s
/// existing per-outer-atom pattern -- and resolves every embedded stereocenter *in
/// place* within that same digraph via [`resolve_chirality`], never re-rooting a
/// second digraph at any embedded atom. Loops over every group produced by
/// Rules-1a/1b/2 (mirroring the reference engine's `resolve_outer_root`), attempting
/// the Rule 4b tiebreak on each clean 2-way tie found -- not restricted to exactly one
/// tied group the way [`crate::assign::assign_one_with_rule5`] is, since that
/// restriction isn't part of what the reference engine validated 72/72.
pub(crate) fn assign_one_with_rule4b(
    mol: &Molecule,
    idx: AtomIdx,
    chirality: Chirality,
    stereo_order: &[u32],
    budget: CipBudget,
    kekule: Option<&(Molecule, MancudeContext)>,
) -> Result<Option<CipCode>, SkipReason> {
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
    let mut groups =
        rank_children(&mut graph, &root_children, &mut ctx).map_err(map_compare_err)?;

    let mut cache = HashMap::new();
    for gi in 0..groups.len() {
        if groups[gi].len() == 2 {
            let (a, b) = (groups[gi][0], groups[gi][1]);
            if let Some(ord) = break_tie_rule4b(&mut graph, mol, a, b, budget, &mut cache)
                .map_err(map_compare_err)?
            {
                let (higher, lower) = match ord {
                    Ordering::Greater => (a, b),
                    _ => (b, a),
                };
                if ord != Ordering::Equal {
                    groups[gi] = vec![lower];
                    groups.insert(gi, vec![higher]);
                }
            }
        } else if groups[gi].len() > 2 {
            return Err(SkipReason::Tied);
        }
    }

    let Some(is_r) = resolve_is_r_from_groups(&groups, &position_nodes, chirality) else {
        return Err(SkipReason::Tied);
    };
    Ok(Some(if is_r { CipCode::R } else { CipCode::S }))
}

/// Milestone 4B-2's Rule 4b refinement -- see module docs for scope. Only ever touches
/// atoms Pass 1 (`assign_all`) left [`SkipReason::Tied`]; every other atom is carried
/// through unchanged. Mirrors [`crate::assign::apply_rule5_pass`]'s exact shape and,
/// crucially, runs *before* it in [`crate::assign::assign_cip_accurate_experimental`]
/// (Rule 4b precedes Rule 5 in CIP rule order) -- Rule 5's own `provisional` map
/// (built from whatever `AccurateCipAssignment` it's handed) sees Rule 4b's
/// newly-resolved atoms for free, no changes needed inside `apply_rule5_pass` itself.
pub(crate) fn apply_rule4b_pass(
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

        match assign_one_with_rule4b(mol, idx, atom.chirality, stereo_order, budget, kekule) {
            Ok(Some(code)) => assignments.push((idx, code)),
            _ => skipped.push((idx, reason)),
        }
    }

    AccurateCipAssignment {
        assignments,
        skipped,
    }
}
