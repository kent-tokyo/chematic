//! Compare auxiliary tetrahedral descriptors in the constitutional digraph.
//!
//! Use Rules 4a, 4b, 4c and then 5 on constitutionally tied ligands. The
//! auxiliaries belong to occurrences in this rooted digraph, so different
//! occurrences of the same ring atom can have different descriptors.

use std::cell::OnceCell;
use std::cmp::Ordering;
use std::collections::HashMap;

use chematic_core::{CipCode, Molecule};
use smallvec::SmallVec;

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
    has_pseudo_descriptor: bool,
    profiles: [OnceCell<Vec<u8>>; 6],
}

fn describe_branch(
    graph: &mut CipDigraph,
    mol: &Molecule,
    start: NodeId,
    auxiliaries: &mut HashMap<NodeId, Option<CipCode>>,
    ctx: &mut CompareContext,
) -> Result<StereoBranch, CipCompareError> {
    let mut level = vec![vec![start]];
    let mut nodes = Vec::new();
    let mut reference_votes = [0; 2];
    let mut has_pseudo_descriptor = false;

    while !level.is_empty() {
        let mut next = Vec::new();
        for group in level {
            let mut votes = [0; 2];
            let mut child_groups: Vec<Vec<NodeId>> = Vec::new();
            for node in group {
                match resolve_chirality(graph, mol, node, auxiliaries)? {
                    Some(CipCode::S) => votes[0] += 1,
                    Some(CipCode::R) => votes[1] += 1,
                    Some(CipCode::LowerR | CipCode::LowerS) => has_pseudo_descriptor = true,
                    _ => {}
                }
                let children = graph
                    .expand_child_ids(node)
                    .map_err(CipCompareError::Digraph)?;
                let ranks = rank_children(graph, &children, ctx)?;
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
        has_pseudo_descriptor,
        profiles: Default::default(),
    })
}

#[derive(Clone, Copy)]
enum DescriptorOrder {
    Class,
    LikeR,
    LikeS,
    Pseudo,
    AbsoluteR,
    AbsoluteS,
}

impl DescriptorOrder {
    fn key(self, code: Option<CipCode>) -> u8 {
        match self {
            Self::Class => match code {
                Some(CipCode::R | CipCode::S) => 2,
                Some(CipCode::LowerR | CipCode::LowerS) => 1,
                _ => 0,
            },
            Self::LikeR | Self::AbsoluteR => u8::from(code == Some(CipCode::R)),
            Self::LikeS | Self::AbsoluteS => u8::from(code == Some(CipCode::S)),
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
    ) -> &[u8] {
        self.profiles[order as usize].get_or_init(|| {
            let has_absolute = self.reference_votes != [0; 2];
            let has_keys = match order {
                DescriptorOrder::Class => has_absolute || self.has_pseudo_descriptor,
                DescriptorOrder::Pseudo => self.has_pseudo_descriptor,
                _ => has_absolute,
            };
            if has_keys {
                self.build_profile(auxiliaries, order)
            } else {
                // With no matching descriptors every traversal gives the same
                // sequence, including its length.
                vec![0; self.nodes.len()]
            }
        })
    }

    fn build_profile(
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
                        let pseudo = if matches!(
                            order,
                            DescriptorOrder::AbsoluteR | DescriptorOrder::AbsoluteS
                        ) {
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

    fn like_profile(&self, auxiliaries: &HashMap<NodeId, Option<CipCode>>) -> Option<&[u8]> {
        let [s, r] = self.reference_votes;
        match r.cmp(&s) {
            Ordering::Greater => Some(self.profile(auxiliaries, DescriptorOrder::LikeR)),
            Ordering::Less => Some(self.profile(auxiliaries, DescriptorOrder::LikeS)),
            Ordering::Equal if r != 0 => Some(
                self.profile(auxiliaries, DescriptorOrder::LikeR)
                    .max(self.profile(auxiliaries, DescriptorOrder::LikeS)),
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
    compare_tied_ligands_with_reference(graph, mol, left, right, auxiliaries, limit, CipCode::R)
}

fn compare_tied_ligands_with_reference(
    graph: &mut CipDigraph,
    mol: &Molecule,
    left: NodeId,
    right: NodeId,
    auxiliaries: &mut HashMap<NodeId, Option<CipCode>>,
    limit: StereoRuleLimit,
    reference: CipCode,
) -> Result<Option<(Ordering, bool)>, CipCompareError> {
    let mut ctx = CompareContext::new();
    let left = describe_branch(graph, mol, left, auxiliaries, &mut ctx)?;
    let right = describe_branch(graph, mol, right, auxiliaries, &mut ctx)?;
    Ok(compare_branches(
        &left,
        &right,
        auxiliaries,
        limit,
        reference,
    ))
}

fn compare_branches(
    left: &StereoBranch,
    right: &StereoBranch,
    auxiliaries: &HashMap<NodeId, Option<CipCode>>,
    limit: StereoRuleLimit,
    reference: CipCode,
) -> Option<(Ordering, bool)> {
    let compare = |order| {
        left.profile(auxiliaries, order)
            .cmp(right.profile(auxiliaries, order))
    };

    let class = compare(DescriptorOrder::Class);
    if class != Ordering::Equal {
        return Some((class, false));
    }
    if let (Some(left), Some(right)) = (
        left.like_profile(auxiliaries),
        right.like_profile(auxiliaries),
    ) {
        let like = left.cmp(right);
        if like != Ordering::Equal {
            return Some((like, false));
        }
    }
    let pseudo = compare(DescriptorOrder::Pseudo);
    if pseudo != Ordering::Equal {
        return Some((pseudo, false));
    }

    if limit == StereoRuleLimit::Rule5 {
        let (primary, opposite) = if reference == CipCode::R {
            (DescriptorOrder::AbsoluteR, DescriptorOrder::AbsoluteS)
        } else {
            (DescriptorOrder::AbsoluteS, DescriptorOrder::AbsoluteR)
        };
        let order = compare(primary);
        if order != Ordering::Equal {
            return Some((order, order != compare(opposite)));
        }
    }
    None
}

/// Apply the auxiliary rules within each existing tie without changing the
/// order of ligand groups already established by earlier rules.
pub(crate) fn refine_stereo_groups(
    graph: &mut CipDigraph,
    mol: &Molecule,
    groups: Vec<Vec<NodeId>>,
    positions: &[NodeId],
    cache: &mut HashMap<NodeId, Option<CipCode>>,
) -> Result<(Vec<Vec<NodeId>>, usize), CipCompareError> {
    refine_stereo_groups_with_reference(graph, mol, groups, positions, cache, CipCode::R)
}

pub(crate) fn refine_stereo_groups_with_reference(
    graph: &mut CipDigraph,
    mol: &Molecule,
    groups: Vec<Vec<NodeId>>,
    positions: &[NodeId],
    cache: &mut HashMap<NodeId, Option<CipCode>>,
    reference: CipCode,
) -> Result<(Vec<Vec<NodeId>>, usize), CipCompareError> {
    let mut result = Vec::with_capacity(groups.len());
    let mut reflection_flips = 0;
    for group in groups {
        let (ranked, flips) = rank_stereo_group(graph, mol, group, positions, cache, reference)?;
        result.extend(ranked);
        reflection_flips += flips;
    }
    Ok((result, reflection_flips))
}

/// Refine a constitutional tie over all physical ligands in the group.
/// Nonphysical duplicate nodes stay in the last refined group.
fn rank_stereo_group(
    graph: &mut CipDigraph,
    mol: &Molecule,
    mut group: Vec<NodeId>,
    positions: &[NodeId],
    cache: &mut HashMap<NodeId, Option<CipCode>>,
    reference: CipCode,
) -> Result<(Vec<Vec<NodeId>>, usize), CipCompareError> {
    let physical: SmallVec<[NodeId; 4]> = group
        .iter()
        .copied()
        .filter(|n| positions.contains(n))
        .collect();
    if physical.len() < 2 {
        return Ok((vec![group], 0));
    }
    if let [left, right] = physical.as_slice() {
        let ranked = compare_tied_ligands_with_reference(
            graph,
            mol,
            *left,
            *right,
            cache,
            StereoRuleLimit::Rule5,
            reference,
        )?;
        if let Some((order, pseudo)) = ranked {
            let higher = if order == Ordering::Greater {
                *left
            } else {
                *right
            };
            group.retain(|&node| node != higher);
            return Ok((vec![vec![higher], group], usize::from(pseudo)));
        }
        return Ok((vec![group], 0));
    }
    let mut wins = vec![0; physical.len()];
    let mut comparisons = Vec::new();
    let mut reflection_flips = 0;
    let mut ctx = CompareContext::new();
    let branches: Vec<_> = physical
        .iter()
        .map(|&node| describe_branch(graph, mol, node, cache, &mut ctx))
        .collect::<Result<_, _>>()?;
    for (i, left) in branches.iter().enumerate() {
        for (j, right) in branches.iter().enumerate().skip(i + 1) {
            let (order, pseudo) =
                compare_branches(left, right, cache, StereoRuleLimit::Rule5, reference)
                    .unwrap_or((Ordering::Equal, false));
            match order {
                Ordering::Greater => wins[i] += 1,
                Ordering::Less => wins[j] += 1,
                Ordering::Equal => {}
            }
            reflection_flips += usize::from(pseudo);
            comparisons.push((i, j, order));
        }
    }
    // Equal scores are a tie only if all pairwise decisions agree with them.
    if comparisons
        .iter()
        .any(|&(i, j, order)| wins[i].cmp(&wins[j]) != order)
    {
        return Err(CipCompareError::InvalidDigraph(
            "inconsistent auxiliary ligand ordering".to_owned(),
        ));
    }
    if wins.iter().all(|&score| score == 0) {
        return Ok((vec![group], 0));
    }
    let mut ranked: Vec<_> = physical.iter().copied().zip(wins).collect();
    ranked.sort_by_key(|&(_, score)| std::cmp::Reverse(score));
    let mut groups: Vec<Vec<_>> = ranked
        .chunk_by(|&(_, a), &(_, b)| a == b)
        .map(|class| class.iter().map(|&(node, _)| node).collect())
        .collect();
    groups
        .last_mut()
        .expect("a ranked group exists")
        .extend(group.into_iter().filter(|node| !physical.contains(node)));
    Ok((groups, reflection_flips))
}
