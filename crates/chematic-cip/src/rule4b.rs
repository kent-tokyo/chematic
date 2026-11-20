//! Compare auxiliary tetrahedral descriptors in the constitutional digraph.
//!
//! Use Rules 4a, 4b, 4c and then 5 on constitutionally tied ligands. The
//! auxiliaries belong to occurrences in this rooted digraph, so different
//! occurrences of the same ring atom can have different descriptors.

use std::cmp::Ordering;
use std::collections::HashMap;

use chematic_core::{CipCode, Molecule};

use crate::compare::{CipCompareError, CompareContext, rank_children};
use crate::digraph::CipDigraph;
use crate::node::NodeId;
use crate::resolver::resolve_chirality;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum StereoRuleLimit {
    Rule4c,
    Rule5,
}

/// Nodes in breadth-first order and their constitutionally ranked child groups.
struct StereoBranch {
    nodes: Vec<(NodeId, Vec<Vec<NodeId>>)>,
    reference_votes: [usize; 2],
}

fn describe_branch(
    graph: &mut CipDigraph,
    mol: &Molecule,
    start: NodeId,
    auxiliaries: &mut HashMap<NodeId, Option<CipCode>>,
) -> Result<StereoBranch, CipCompareError> {
    let mut level = vec![vec![start]];
    let mut nodes = Vec::new();
    let mut reference_votes = [0; 2];
    let mut ctx = CompareContext::new();

    while !level.is_empty() {
        let mut next = Vec::new();
        for group in level {
            let mut votes = [0; 2];
            let mut child_groups: Vec<Vec<NodeId>> = Vec::new();
            for node in group {
                match resolve_chirality(graph, mol, node, auxiliaries)? {
                    Some(CipCode::S) => votes[0] += 1,
                    Some(CipCode::R) => votes[1] += 1,
                    _ => {}
                }
                let children = graph
                    .expand_child_ids(node)
                    .map_err(CipCompareError::Digraph)?;
                let ranks = rank_children(graph, &children, &mut ctx)?;
                // Keep corresponding constitutional groups together when
                // looking for the first group containing reference descriptors.
                for (index, children) in ranks.iter().enumerate() {
                    if let Some(group) = child_groups.get_mut(index) {
                        group.extend(children);
                    } else {
                        child_groups.push(children.clone());
                    }
                }
                nodes.push((node, ranks));
            }
            if reference_votes == [0; 2] {
                reference_votes = votes;
            }
            next.extend(child_groups);
        }
        level = next;
    }

    Ok(StereoBranch {
        nodes,
        reference_votes,
    })
}

#[derive(Clone, Copy)]
enum DescriptorOrder {
    Class,
    Like(CipCode),
    Pseudo,
    Absolute(CipCode),
}

impl DescriptorOrder {
    fn key(self, code: Option<CipCode>) -> u8 {
        match self {
            Self::Class => match code {
                Some(CipCode::R | CipCode::S) => 2,
                Some(CipCode::LowerR | CipCode::LowerS) => 1,
                _ => 0,
            },
            Self::Like(reference) | Self::Absolute(reference) => u8::from(code == Some(reference)),
            Self::Pseudo => match code {
                Some(CipCode::LowerR) => 2,
                Some(CipCode::LowerS) => 1,
                _ => 0,
            },
        }
    }
}

fn descriptor_sequence(
    start: NodeId,
    children: &HashMap<NodeId, Vec<NodeId>>,
    auxiliaries: &HashMap<NodeId, Option<CipCode>>,
    order: DescriptorOrder,
) -> Vec<u8> {
    let mut queue = vec![start];
    let mut keys = Vec::new();
    while let Some(&node) = queue.get(keys.len()) {
        keys.push(order.key(auxiliaries[&node]));
        queue.extend(&children[&node]);
    }
    keys
}

impl StereoBranch {
    fn profile(
        &self,
        auxiliaries: &HashMap<NodeId, Option<CipCode>>,
        order: DescriptorOrder,
    ) -> Vec<u8> {
        let mut ordered_children = HashMap::with_capacity(self.nodes.len());
        // Resolve sibling ordering bottom-up and preserve each parent's
        // child order during the final breadth-first walk.
        for (node, groups) in self.nodes.iter().rev() {
            let mut children = Vec::new();
            for group in groups {
                if let [child] = group.as_slice() {
                    children.push(*child);
                    continue;
                }
                let mut ranked: Vec<_> = group
                    .iter()
                    .map(|&child| {
                        let sequence = |order| {
                            descriptor_sequence(child, &ordered_children, auxiliaries, order)
                        };
                        let class = sequence(DescriptorOrder::Class);
                        let pseudo = if matches!(order, DescriptorOrder::Absolute(_)) {
                            sequence(DescriptorOrder::Pseudo)
                        } else {
                            Vec::new()
                        };
                        let priority = sequence(order);
                        (child, (class, pseudo, priority))
                    })
                    .collect();
                ranked.sort_by(|(_, a), (_, b)| b.cmp(a));
                children.extend(ranked.into_iter().map(|(child, _)| child));
            }
            ordered_children.insert(*node, children);
        }
        descriptor_sequence(self.nodes[0].0, &ordered_children, auxiliaries, order)
    }

    fn like_profile(&self, auxiliaries: &HashMap<NodeId, Option<CipCode>>) -> Option<Vec<u8>> {
        let [s, r] = self.reference_votes;
        match r.cmp(&s) {
            Ordering::Greater => Some(self.profile(auxiliaries, DescriptorOrder::Like(CipCode::R))),
            Ordering::Less => Some(self.profile(auxiliaries, DescriptorOrder::Like(CipCode::S))),
            Ordering::Equal if r != 0 => Some(
                self.profile(auxiliaries, DescriptorOrder::Like(CipCode::R))
                    .max(self.profile(auxiliaries, DescriptorOrder::Like(CipCode::S))),
            ),
            Ordering::Equal => None,
        }
    }
}

/// Return the order and whether its deciding comparison is pseudoasymmetric.
pub(crate) fn compare_tied_ligands(
    graph: &mut CipDigraph,
    mol: &Molecule,
    left: NodeId,
    right: NodeId,
    auxiliaries: &mut HashMap<NodeId, Option<CipCode>>,
    limit: StereoRuleLimit,
) -> Result<Option<(Ordering, bool)>, CipCompareError> {
    let left = describe_branch(graph, mol, left, auxiliaries)?;
    let right = describe_branch(graph, mol, right, auxiliaries)?;
    let compare = |order| {
        left.profile(auxiliaries, order)
            .cmp(&right.profile(auxiliaries, order))
    };

    let class = compare(DescriptorOrder::Class);
    if class != Ordering::Equal {
        return Ok(Some((class, false)));
    }
    if let (Some(left), Some(right)) = (
        left.like_profile(auxiliaries),
        right.like_profile(auxiliaries),
    ) {
        let like = left.cmp(&right);
        if like != Ordering::Equal {
            return Ok(Some((like, false)));
        }
    }
    let pseudo = compare(DescriptorOrder::Pseudo);
    if pseudo != Ordering::Equal {
        return Ok(Some((pseudo, false)));
    }

    if limit == StereoRuleLimit::Rule5 {
        let r = compare(DescriptorOrder::Absolute(CipCode::R));
        if r != Ordering::Equal {
            let s = compare(DescriptorOrder::Absolute(CipCode::S));
            return Ok(Some((r, r != s)));
        }
    }
    Ok(None)
}
