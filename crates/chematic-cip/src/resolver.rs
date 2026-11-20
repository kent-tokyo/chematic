//! Resolve auxiliary tetrahedral descriptors and refine constitutional ties.

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
use crate::mancude::PreparedKekule;
use crate::node::{CipNodeKind, NodeId};
use crate::rule4b::{StereoRuleLimit, compare_tied_ligands};

/// The core recursive, memoized, in-place resolver. `node_id` must be a real `Atom`
/// node (not the true digraph root -- the root has no "back to root" ligand and is
/// handled by [`assign_one_with_rule4b`] directly, mirroring the reference engine's
/// own `resolve_outer_root`/`resolve_chirality` split).
pub(crate) fn resolve_chirality(
    graph: &mut CipDigraph,
    mol: &Molecule,
    node_id: NodeId,
    cache: &mut HashMap<NodeId, Option<CipCode>>,
) -> Result<Option<CipCode>, CipCompareError> {
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
    let groups = rank_children(graph, &forward_children, &mut ctx)?;
    let mut groups = crate::rule3::refine(graph, groups, &mut ctx)?;

    // Refine the forward ligands before inserting the back-to-root ligand.
    let mut pseudo_count = 0;
    for gi in (0..groups.len()).rev() {
        if let [a, b] = groups[gi].as_slice() {
            let (a, b) = (*a, *b);
            let result = compare_tied_ligands(graph, mol, a, b, cache, StereoRuleLimit::Rule5)?;
            if let Some((ord, is_pseudo)) = result {
                let (higher, lower) = if ord == Ordering::Greater {
                    (a, b)
                } else {
                    (b, a)
                };
                groups[gi] = vec![lower];
                groups.insert(gi, vec![higher]);
                pseudo_count += usize::from(is_pseudo);
            }
        } else if groups[gi].len() > 2 {
            cache.insert(node_id, None);
            return Ok(None);
        }
    }

    // Insert the back-to-root ligand at its Rule-1a-only rank among `groups`. Compare
    // sphere-1-to-sphere-1: the back ligand's sphere 1 is `ancestor` itself; the
    // forward ligand's sphere 1 is `rep` itself -- NOT `rep`'s children (that would be
    // off by one sphere).
    let mut insert_at = groups.len();
    for (gi, group) in groups.iter().enumerate() {
        let rep = group[0];
        let back = crate::auxiliary::BackItem::Ascending {
            ancestor: parent_id,
            came_from: node_id,
        };
        let ord = crate::auxiliary::compare_rule1a_only(graph, back, rep, &mut ctx)?;
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
    let chirality = mol.atom(atom_idx).chirality;
    let result = if stereo_order.len() == 3 {
        // A saturated ring centre has a lowest-priority phantom ligand.
        // Use the same parity for auxiliary and top-level descriptors.
        if crate::assign::is_saturated_ring_lone_pair_centre(mol, atom_idx) {
            crate::assign::lone_pair_is_r(&final_groups, &position_nodes, chirality)
                .ok()
                .flatten()
        } else {
            None
        }
    } else {
        resolve_is_r_from_groups(&final_groups, &position_nodes, chirality)
    };
    let result = result.map(|is_r| match (is_r, pseudo_count == 1) {
        (true, false) => CipCode::R,
        (false, false) => CipCode::S,
        (true, true) => CipCode::LowerR,
        (false, true) => CipCode::LowerS,
    });
    cache.insert(node_id, result);
    Ok(result)
}

/// Try the reflection-invariant auxiliary rules on a constitutional tie.
pub(crate) fn assign_one_with_rule4b(
    mol: &Molecule,
    idx: AtomIdx,
    chirality: Chirality,
    stereo_order: &[u32],
    budget: CipBudget,
    kekule: Option<&PreparedKekule<'_>>,
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
    let groups = rank_children(&mut graph, &root_children, &mut ctx).map_err(map_compare_err)?;
    let mut groups = crate::rule3::refine(&mut graph, groups, &mut ctx).map_err(map_compare_err)?;

    let mut cache = HashMap::new();
    for gi in (0..groups.len()).rev() {
        if groups[gi].len() == 2 {
            let (a, b) = (groups[gi][0], groups[gi][1]);
            if let Some((ord, _)) =
                compare_tied_ligands(&mut graph, mol, a, b, &mut cache, StereoRuleLimit::Rule4c)
                    .map_err(map_compare_err)?
            {
                let (higher, lower) = match ord {
                    Ordering::Greater => (a, b),
                    _ => (b, a),
                };
                groups[gi] = vec![lower];
                groups.insert(gi, vec![higher]);
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

/// Refine centers with a constitutional tie. Other results pass through unchanged.
pub(crate) fn apply_rule4b_pass(
    mol: &Molecule,
    budget: CipBudget,
    kekule: Option<&PreparedKekule<'_>>,
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
