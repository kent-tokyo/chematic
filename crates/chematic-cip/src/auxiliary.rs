//! Milestone 4B-2 production port of the "back to root" ligand from the validated
//! reference engine (`examples/rule4b_bottom_up.rs`, 72/72 across 4 oracle corpora
//! including mirrors -- see `docs/rfcs/cip_accurate_rfc.md`). Mechanical port: algorithm
//! unchanged from the example.
//!
//! Per Hanson, Musacchio, Mayfield et al. 2018 (*J. Chem. Inf. Model.* 58(9),
//! 1755-1765): "the priority of a ligand leading back to the digraph root will always
//! be ranked by Rule 1a, with no need to consider auxiliary centers... the path back
//! to the root is always unique in connectivity and atomic numbers." [`BackItem`]
//! implements this as a synthetic frontier that walks *up* through the existing parent
//! chain, reusing already-built off-path subtrees as-is (no rebuilding, no fresh digraph
//! rooted at the embedded atom -- that re-rooted architecture is what produced the
//! pair-antisymmetry bug this milestone's reference engine fixed).
//! [`compare_rule1a_only`] compares this frontier against a real forward ligand under
//! Rule 1a with the same hierarchical, branch-by-branch exploration as
//! `compare.rs` (issue #634: an earlier version pooled each sphere into one sorted
//! multiset, which is not how a hierarchical digraph is explored and mis-ranked the
//! back ligand of row 4480's embedded atoms).

use std::cmp::Ordering;

use crate::compare::{CipCompareError, CompareContext, rank_children, ranked_child_ids};
use crate::digraph::CipDigraph;
use crate::node::NodeId;
use crate::rational::{AtomicNumberKey, cmp_atomic_number_key};

/// A position in the back-to-root ligand's growing frontier.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum BackItem {
    /// A node already present in the outer digraph -- normal forward expansion from
    /// here on via `expand_children`, no special handling.
    Reused(NodeId),
    /// The walk continues upward: `ancestor` is this sphere's atom; `came_from` is the
    /// specific child of `ancestor` to exclude when expanding (the one leading back
    /// down toward the atom under resolution -- never re-enter its subtree).
    Ascending { ancestor: NodeId, came_from: NodeId },
}

fn back_item_atomic_number(graph: &CipDigraph, item: BackItem) -> AtomicNumberKey {
    let node = match item {
        BackItem::Reused(n) => n,
        BackItem::Ascending { ancestor, .. } => ancestor,
    };
    graph.node(node).atomic_number
}

/// One position of a sphere in [`compare_rule1a_only`]'s level walk.
#[derive(Clone, Copy)]
enum Slot {
    Back(BackItem),
    Forward(NodeId),
    /// Padding for a position whose counterpart has more substituents.
    Phantom,
}

fn slot_atomic_number(graph: &CipDigraph, slot: Slot) -> AtomicNumberKey {
    match slot {
        Slot::Back(item) => back_item_atomic_number(graph, item),
        Slot::Forward(n) => graph.node(n).atomic_number,
        Slot::Phantom => AtomicNumberKey::Integral(0),
    }
}

/// The sphere-(k+1) items below `item`, in rank order.
///
/// A reused node's children are ranked exactly as the forward comparison ranks
/// them ([`ranked_child_ids`]). An ascending item's children are its off-path
/// children, ranked the same way, plus the next ascending item (its own parent),
/// placed at its Rule 1a rank among them (after any group it ties with).
fn ranked_back_children(
    graph: &mut CipDigraph,
    item: BackItem,
    ctx: &mut CompareContext,
) -> Result<Vec<BackItem>, CipCompareError> {
    match item {
        BackItem::Reused(n) => Ok(ranked_child_ids(graph, n, ctx)?
            .into_iter()
            .map(BackItem::Reused)
            .collect()),
        BackItem::Ascending {
            ancestor,
            came_from,
        } => {
            let off_path: Vec<NodeId> = graph
                .expand_children(ancestor)
                .map_err(CipCompareError::Digraph)?
                .into_iter()
                .filter(|&c| c != came_from)
                .collect();
            let groups = rank_children(graph, &off_path, ctx)?;
            let up = graph
                .node(ancestor)
                .parent
                .map(|grandparent| BackItem::Ascending {
                    ancestor: grandparent,
                    came_from: ancestor,
                });
            let mut insert_at = groups.len();
            if let Some(up) = up {
                for (gi, group) in groups.iter().enumerate() {
                    if compare_rule1a_only(graph, up, group[0], ctx)? == Ordering::Greater {
                        insert_at = gi;
                        break;
                    }
                }
            }
            let mut out = Vec::with_capacity(off_path.len() + 1);
            for (gi, group) in groups.iter().enumerate() {
                if gi == insert_at {
                    out.extend(up);
                }
                out.extend(group.iter().copied().map(BackItem::Reused));
            }
            if insert_at == groups.len() {
                out.extend(up);
            }
            Ok(out)
        }
    }
}

/// Rule 1a comparison of the back-to-root ligand (`lhs`) against a forward
/// ligand (`rhs`); `Greater` means the back-to-root ligand outranks `rhs`.
///
/// Hierarchical, like [`crate::compare::compare_ligands`]: sphere by sphere,
/// each sphere compared position by position, where a sphere lists the
/// previous sphere's items' substituents set by set in the order of those
/// items' ranks (each set ranked, and padded with phantom atoms to the
/// counterpart's size). Atoms of one sphere are *not* pooled into a single
/// multiset: that loses which branch an atom belongs to, and the branch explored
/// first must decide first (IUPAC 2013 P-92.1.4). For example, from the
/// embedded atom of `CO[C@@H]1[C@@H](N)[C@@H](OC)[C@@H](O)[C@H]1O` issue #634
/// row 4480 compares a back ligand `C(OC)(C(N)...)` against a forward ligand
/// `C(O)(C(OC)...)`: pooled, the third sphere is {N,C,C,H} against {O,C,H,H}
/// and the forward ligand wins; hierarchically, the highest-ranked branches
/// (the oxygens) are compared first, {C} against {H}, and the back ligand wins.
pub(crate) fn compare_rule1a_only(
    graph: &mut CipDigraph,
    lhs: BackItem,
    rhs: NodeId,
    ctx: &mut CompareContext,
) -> Result<Ordering, CipCompareError> {
    let mut left = vec![Slot::Back(lhs)];
    let mut right = vec![Slot::Forward(rhs)];
    loop {
        ctx.recursive_calls += 1;
        if ctx.recursive_calls > ctx.max_recursive_calls {
            return Err(CipCompareError::BudgetExceeded {
                expanded_nodes: graph.nodes().len(),
                recursive_calls: ctx.recursive_calls,
            });
        }
        debug_assert_eq!(left.len(), right.len());
        for (&l, &r) in left.iter().zip(&right) {
            let ord =
                cmp_atomic_number_key(slot_atomic_number(graph, l), slot_atomic_number(graph, r));
            if ord != Ordering::Equal {
                return Ok(ord);
            }
        }
        if left.is_empty() {
            return Ok(Ordering::Equal);
        }
        let mut next_left = Vec::new();
        let mut next_right = Vec::new();
        for (&l, &r) in left.iter().zip(&right) {
            let lc: Vec<Slot> = match l {
                Slot::Back(item) => ranked_back_children(graph, item, ctx)?
                    .into_iter()
                    .map(Slot::Back)
                    .collect(),
                Slot::Forward(n) => ranked_child_ids(graph, n, ctx)?
                    .into_iter()
                    .map(Slot::Forward)
                    .collect(),
                Slot::Phantom => Vec::new(),
            };
            let rc: Vec<Slot> = match r {
                Slot::Back(item) => ranked_back_children(graph, item, ctx)?
                    .into_iter()
                    .map(Slot::Back)
                    .collect(),
                Slot::Forward(n) => ranked_child_ids(graph, n, ctx)?
                    .into_iter()
                    .map(Slot::Forward)
                    .collect(),
                Slot::Phantom => Vec::new(),
            };
            for j in 0..lc.len().max(rc.len()) {
                next_left.push(lc.get(j).copied().unwrap_or(Slot::Phantom));
                next_right.push(rc.get(j).copied().unwrap_or(Slot::Phantom));
            }
        }
        left = next_left;
        right = next_right;
    }
}
