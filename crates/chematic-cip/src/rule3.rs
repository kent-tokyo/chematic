//! Sequence Rule 3 for alkene geometry in the rooted CIP digraph.

use std::cmp::Ordering;
use std::collections::HashMap;

use chematic_core::ez_markers::EzFactEnd;
use chematic_core::{AtomIdx, BondOrder};
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

use crate::auxiliary::{BackItem, compare_rule1a_only};
use crate::compare::{BranchComparison, CipCompareError, CompareContext, compare_ligands};
use crate::digraph::CipDigraph;
use crate::ligand_profile::{Profile, collect_branches, sphere_profiles};
use crate::node::{CipNodeKind, NodeId};
use crate::rule4b::{StereoRuleLimit, compare_tied_ligands};

fn ligand_atom(graph: &CipDigraph, node: NodeId) -> Option<AtomIdx> {
    match graph.node(node).kind {
        CipNodeKind::Atom { atom_idx } => Some(atom_idx),
        CipNodeKind::RingDuplicate { closure_atom, .. } => Some(closure_atom),
        _ => None,
    }
}

fn preferred_ligand(
    graph: &mut CipDigraph,
    node: NodeId,
    paired_end: AtomIdx,
    context: &mut CompareContext,
    keys: &FxHashMap<NodeId, u8>,
    include_auxiliary: bool,
) -> Result<Option<(NodeId, bool)>, CipCompareError> {
    let children = graph
        .expand_child_ids(node)
        .map_err(CipCompareError::Digraph)?;
    let mut forward = SmallVec::<[NodeId; 2]>::new();
    forward.extend(children.into_iter().filter(|&child| {
        !matches!(
            graph.node(child).kind,
            CipNodeKind::MultipleBondDuplicate { .. }
        ) && ligand_atom(graph, child) != Some(paired_end)
    }));
    let back = graph
        .node(node)
        .parent
        .filter(|&parent| ligand_atom(graph, parent) != Some(paired_end));
    Ok(match (back, forward.as_slice()) {
        (None, [only]) => Some((*only, false)),
        (Some(only), []) => Some((only, false)),
        (None, [a, b]) => prefer_forward_ligand(graph, [*a, *b], context, keys, include_auxiliary)?,
        (Some(ancestor), [child]) => {
            let back = BackItem::Ascending {
                ancestor,
                came_from: node,
            };
            match compare_rule1a_only(graph, back, *child, context)? {
                Ordering::Greater => Some((ancestor, false)),
                Ordering::Less => Some((*child, false)),
                Ordering::Equal => None,
            }
        }
        _ => None,
    })
}

fn prefer_forward_ligand(
    graph: &mut CipDigraph,
    [a, b]: [NodeId; 2],
    context: &mut CompareContext,
    keys: &FxHashMap<NodeId, u8>,
    include_auxiliary: bool,
) -> Result<Option<(NodeId, bool)>, CipCompareError> {
    let order = match compare_ligands(graph, a, b, context)? {
        BranchComparison::Higher => Ordering::Greater,
        BranchComparison::Lower => Ordering::Less,
        BranchComparison::Equal => compare_known_descriptors(graph, a, b, context, keys)?,
        BranchComparison::Unresolved => return Ok(None),
    };
    let ranked = if order == Ordering::Equal && include_auxiliary {
        compare_tied_ligands(
            graph,
            graph.molecule(),
            a,
            b,
            &mut HashMap::new(),
            StereoRuleLimit::Rule5,
        )?
    } else {
        (order != Ordering::Equal).then_some((order, false))
    };
    Ok(ranked.map(|(order, pseudo)| (if order == Ordering::Greater { a } else { b }, pseudo)))
}

fn compare_known_descriptors(
    graph: &mut CipDigraph,
    left: NodeId,
    right: NodeId,
    context: &mut CompareContext,
    keys: &FxHashMap<NodeId, u8>,
) -> Result<Ordering, CipCompareError> {
    if keys.is_empty() {
        return Ok(Ordering::Equal);
    }
    let branches = collect_branches(graph, [left, right].into_iter(), context)?;
    let profiles = sphere_profiles(&branches, |node| keys.get(&node).copied().unwrap_or(0));
    Ok(profiles[&left].cmp(&profiles[&right]))
}

fn ranked_side(graph: &CipDigraph, end: EzFactEnd, ligand: NodeId) -> bool {
    let reference = graph.molecule().bond(end.reference).other(end.atom);
    end.up ^ (ligand_atom(graph, ligand) != reference)
}

fn sequence_key(
    graph: &mut CipDigraph,
    node: NodeId,
    children: &[Vec<NodeId>],
    context: &mut CompareContext,
    keys: &FxHashMap<NodeId, u8>,
    include_auxiliary: bool,
) -> Result<Option<u8>, CipCompareError> {
    let CipNodeKind::Atom { atom_idx } = graph.node(node).kind else {
        return Ok(Some(0));
    };
    let mut key = 0;
    for &child in children.iter().flatten() {
        let CipNodeKind::Atom { atom_idx: partner } = graph.node(child).kind else {
            continue;
        };
        let Some((bond_idx, bond)) = graph.molecule().bond_between(atom_idx, partner) else {
            continue;
        };
        let mol = graph.molecule();
        if bond.order != BondOrder::Double
            || mol.atom(atom_idx).aromatic && mol.atom(partner).aromatic
            || [atom_idx, partner].into_iter().any(|end| {
                mol.neighbors(end).any(|(_, adjacent)| {
                    adjacent != bond_idx && mol.bond(adjacent).order == BondOrder::Double
                })
            })
        {
            continue;
        }
        let Some((top_near, pseudo_near)) =
            preferred_ligand(graph, node, partner, context, keys, include_auxiliary)?
        else {
            continue;
        };
        let Some((top_far, pseudo_far)) =
            preferred_ligand(graph, child, atom_idx, context, keys, include_auxiliary)?
        else {
            continue;
        };
        if pseudo_near != pseudo_far {
            continue;
        }
        let Some(mut fact) = graph.alkene_geometry().get(&bond_idx).copied() else {
            // A potentially stereogenic bond with unspecified geometry cannot
            // be ranked below a configured bond as if it were nonstereogenic.
            return Ok(None);
        };
        if fact.ends[0].atom != atom_idx {
            fact.ends.swap(0, 1);
        }
        let cis =
            ranked_side(graph, fact.ends[0], top_near) == ranked_side(graph, fact.ends[1], top_far);
        // Put the auxiliary descriptor at the first end reached in this
        // digraph. The same physical bond can have a different label in
        // another branch when the endpoint priorities change.
        key = if cis { 2 } else { 1 };
    }
    Ok(Some(key))
}

fn profiles(
    graph: &mut CipDigraph,
    ligands: &[NodeId],
    context: &mut CompareContext,
    include_auxiliary: bool,
) -> Result<Option<FxHashMap<NodeId, Profile>>, CipCompareError> {
    let branches = collect_branches(graph, ligands.iter().copied(), context)?;
    let mut keys = FxHashMap::default();
    // An end's forward ligands can themselves be E or Z. Resolve those
    // descendants before deciding which one defines this double bond's side.
    for (node, children) in branches.iter().rev() {
        let Some(key) = sequence_key(graph, *node, children, context, &keys, include_auxiliary)?
        else {
            return Ok(None);
        };
        if key != 0 {
            keys.insert(*node, key);
        }
    }
    Ok((!keys.is_empty())
        .then(|| sphere_profiles(&branches, |node| keys.get(&node).copied().unwrap_or(0))))
}

pub(crate) fn compare(
    graph: &mut CipDigraph,
    left: NodeId,
    right: NodeId,
    context: &mut CompareContext,
    include_auxiliary: bool,
) -> Result<Ordering, CipCompareError> {
    if graph.alkene_geometry().is_empty() {
        return Ok(Ordering::Equal);
    }
    Ok(profiles(graph, &[left, right], context, include_auxiliary)?
        .map_or(Ordering::Equal, |ranks| ranks[&left].cmp(&ranks[&right])))
}

pub(crate) fn refine(
    graph: &mut CipDigraph,
    groups: Vec<Vec<NodeId>>,
    context: &mut CompareContext,
) -> Result<Vec<Vec<NodeId>>, CipCompareError> {
    if groups.iter().all(|group| group.len() < 2) || graph.alkene_geometry().is_empty() {
        return Ok(groups);
    }
    let mut ranked = Vec::with_capacity(groups.len());
    for mut group in groups {
        if group.len() > 1
            && let Some(profiles) = profiles(graph, &group, context, false)?
        {
            group.sort_by(|a, b| profiles[b].cmp(&profiles[a]));
            ranked.extend(
                group
                    .chunk_by(|a, b| profiles[a] == profiles[b])
                    .map(<[_]>::to_vec),
            );
        } else {
            ranked.push(group);
        }
    }
    Ok(ranked)
}

#[cfg(test)]
mod tests {
    use chematic_core::CipCode;

    use crate::{CipBudget, SkipReason, assign_cip_accurate_experimental};

    #[test]
    fn authored_rule3_labels_survive_new_traversals() {
        for (smiles, expected) in [
            (
                r"C/C=C\C/C(=C\[C@H](C=C(C/C=C\C)C/C=C\C)O)/C/C=C/C",
                CipCode::S,
            ), // VS191
            (r"C/C=C(/CCC[C@@H](CCC/C(=C\C)/Cl)O)\Cl", CipCode::R), // VS192
            (r"C/C=C(/CCC[C@H](CCC/C(=C\C)/Cl)O)\Cl", CipCode::S),  // VS193
            (r"C/C=C\CCC[C@@H](CCC/C=C/C)O", CipCode::R),           // VS194
            (r"C/C=C\C[C@@H](C/C=C/C)O", CipCode::R),               // VS195
            (r"C/C=C/1\CC[C@H](C)CC1", CipCode::R),                 // VS246
            (r"O/N=C/1\CC[C@@H](CC1)C", CipCode::R),                // VS247
            (r"CCC[C@H]1C/C(/C1)=C/Br", CipCode::S),                // VS248
        ] {
            let molecule = chematic_smiles::parse(smiles).unwrap();
            let variants = std::iter::once(smiles.to_owned())
                .chain((1..=16).map(|seed| chematic_smiles::random_smiles(&molecule, seed)));
            for variant in variants {
                let input = chematic_smiles::parse(&variant).unwrap();
                let result =
                    assign_cip_accurate_experimental(&input, CipBudget::default_budget()).unwrap();
                let labels: Vec<_> = result.assignments.iter().map(|(_, code)| *code).collect();
                assert_eq!(labels, [expected], "{smiles} -> {variant}");
                assert!(result.skipped.is_empty(), "{smiles} -> {variant}");
            }
        }
    }

    #[test]
    fn four_constitutionally_tied_branches_rank_by_configured_geometry() {
        // The proximal double bond breaks ties first. The four branches are
        // ZZ, ZE, EZ, EE in decreasing priority.
        let smiles = r"[C@](C/C=C\C/C=C\C)(C/C=C\C/C=C/C)(C/C=C/C/C=C\C)C/C=C/C/C=C/C";
        let molecule = chematic_smiles::parse(smiles).unwrap();
        let result =
            assign_cip_accurate_experimental(&molecule, CipBudget::default_budget()).unwrap();
        assert_eq!(
            result.assignments,
            [(chematic_core::AtomIdx(0), CipCode::S)]
        );
    }

    #[test]
    fn missing_or_identical_geometry_does_not_create_a_stereocentre() {
        for smiles in [r"C/C=C\C[C@@H](CC=CC)O", r"C/C=C\C[C@@H](C\C=C/C)O"] {
            let molecule = chematic_smiles::parse(smiles).unwrap();
            let result =
                assign_cip_accurate_experimental(&molecule, CipBudget::default_budget()).unwrap();
            assert!(result.assignments.is_empty(), "{smiles}");
            assert_eq!(
                result.skipped,
                [(chematic_core::AtomIdx(4), SkipReason::Tied)]
            );
        }
    }
}
