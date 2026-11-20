//! Hierarchical sphere profiles for rules that refine constitutional ties.

use rustc_hash::FxHashMap;

use crate::compare::{CipCompareError, CompareContext, rank_children};
use crate::digraph::CipDigraph;
use crate::node::NodeId;

pub(crate) type Branch = Vec<(NodeId, Vec<Vec<NodeId>>)>;
pub(crate) type Profile = Vec<Vec<u8>>;

pub(crate) fn collect_branches(
    graph: &mut CipDigraph,
    starts: impl Iterator<Item = NodeId>,
    context: &mut CompareContext,
) -> Result<Branch, CipCompareError> {
    let mut pending: Vec<_> = starts.collect();
    let mut branches = Vec::new();
    while let Some(node) = pending.pop() {
        let children = graph
            .expand_child_ids(node)
            .map_err(CipCompareError::Digraph)?;
        let ranks = rank_children(graph, &children, context)?;
        pending.extend(ranks.iter().flatten().copied());
        branches.push((node, ranks));
    }
    Ok(branches)
}

pub(crate) fn sphere_profiles(
    branches: &Branch,
    key: impl Fn(NodeId) -> u8,
) -> FxHashMap<NodeId, Profile> {
    let mut profiles = FxHashMap::default();
    // Resolve sibling order from the leaves, retaining each parent's order
    // when the profiles are joined sphere by sphere.
    for (node, groups) in branches.iter().rev() {
        let mut profile = vec![vec![key(*node)]];
        for group in groups {
            let mut children: Vec<Profile> = group
                .iter()
                .map(|child| profiles.remove(child).expect("child precedes parent"))
                .collect();
            children.sort_by(|a, b| b.cmp(a));
            for levels in children {
                for (depth, keys) in levels.into_iter().enumerate() {
                    if profile.len() <= depth + 1 {
                        profile.push(keys);
                    } else {
                        profile[depth + 1].extend(keys);
                    }
                }
            }
        }
        profiles.insert(*node, profile);
    }
    profiles
}
