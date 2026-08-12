//! Milestone 4B-2 production port of the "back to root" ligand from the validated
//! reference engine (`examples/rule4b_bottom_up.rs`, 72/72 across 4 oracle corpora
//! including mirrors -- see `docs/rfcs/cip_accurate_rfc.md`). Mechanical port: algorithm
//! unchanged from the example.
//!
//! Per Hanson, Musacchio, Mayfield et al. 2018 (*J. Chem. Inf. Model.* 58(9),
//! 1755-1765): "the priority of a ligand leading back to the digraph root will always
//! be ranked by Rule 1a, with no need to consider auxiliary centers... the path back
//! to the root is always unique in connectivity and atomic numbers." [`BackItem`]/
//! [`expand_back_item`] implement this as a synthetic frontier that walks *up* through
//! the existing parent chain, reusing already-built off-path subtrees as-is (no
//! rebuilding, no fresh digraph rooted at the embedded atom -- that re-rooted
//! architecture is what produced the pair-antisymmetry bug this milestone's reference
//! engine fixed). [`compare_rule1a_only`] compares this frontier against a real
//! forward ligand via plain sorted-atomic-number shell comparison -- shell-pooling is
//! unsound for the full rule cascade (see `compare.rs` module docs) but is exactly
//! Rule 1a's own classical definition, so it is sound restricted to this one rule.

use std::cmp::Ordering;

use crate::compare::CipCompareError;
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

pub(crate) fn expand_back_item(
    graph: &mut CipDigraph,
    item: BackItem,
) -> Result<Vec<BackItem>, CipCompareError> {
    match item {
        BackItem::Reused(n) => Ok(graph
            .expand_children(n)
            .map_err(CipCompareError::Digraph)?
            .into_iter()
            .map(BackItem::Reused)
            .collect()),
        BackItem::Ascending {
            ancestor,
            came_from,
        } => {
            let mut result: Vec<BackItem> = graph
                .expand_children(ancestor)
                .map_err(CipCompareError::Digraph)?
                .into_iter()
                .filter(|&c| c != came_from)
                .map(BackItem::Reused)
                .collect();
            if let Some(grandparent) = graph.node(ancestor).parent {
                result.push(BackItem::Ascending {
                    ancestor: grandparent,
                    came_from: ancestor,
                });
            }
            Ok(result)
        }
    }
}

/// Rule 1a only, shell-pooled (sound for Rule 1a alone -- it is the classical
/// sum/sorted-multiset-of-atomic-numbers rule; shell-pooling is only unsound once
/// branch-provenance-aware rules like 1b/4b/5 are mixed in, which this function never
/// invokes). `Greater` means the back-to-root ligand (`lhs`) outranks `rhs`.
pub(crate) fn compare_rule1a_only(
    graph: &mut CipDigraph,
    mut lhs: Vec<BackItem>,
    mut rhs: Vec<NodeId>,
) -> Result<Ordering, CipCompareError> {
    loop {
        let mut lhs_keys: Vec<_> = lhs
            .iter()
            .map(|&i| back_item_atomic_number(graph, i))
            .collect();
        let mut rhs_keys: Vec<_> = rhs.iter().map(|&n| graph.node(n).atomic_number).collect();
        lhs_keys.sort_unstable_by(|a, b| cmp_atomic_number_key(*b, *a));
        rhs_keys.sort_unstable_by(|a, b| cmp_atomic_number_key(*b, *a));
        let n = lhs_keys.len().max(rhs_keys.len());
        for i in 0..n {
            let l = lhs_keys.get(i).copied();
            let r = rhs_keys.get(i).copied();
            let ord = match (l, r) {
                (Some(a), Some(b)) => cmp_atomic_number_key(a, b),
                (Some(_), None) => Ordering::Greater,
                (None, Some(_)) => Ordering::Less,
                (None, None) => Ordering::Equal,
            };
            if ord != Ordering::Equal {
                return Ok(ord);
            }
        }
        if lhs.is_empty() && rhs.is_empty() {
            return Ok(Ordering::Equal);
        }
        let mut next_lhs = Vec::new();
        for item in lhs {
            next_lhs.extend(expand_back_item(graph, item)?);
        }
        let mut next_rhs = Vec::new();
        for node in rhs {
            next_rhs.extend(
                graph
                    .expand_children(node)
                    .map_err(CipCompareError::Digraph)?,
            );
        }
        lhs = next_lhs;
        rhs = next_rhs;
    }
}
