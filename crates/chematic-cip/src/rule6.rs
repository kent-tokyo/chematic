//! Reference-ligand ranking proposed by Hanson et al. (2018).
//!
//! Apply after the preceding sequence rules leave a qualifying tie in the
//! ranking of the four ligands, including for auxiliary configurations.

use std::collections::HashMap;

use chematic_core::{AtomIdx, Chirality, CipCode, Molecule};
use rustc_hash::FxHashSet;

use crate::assign::resolve_is_r_from_groups;
use crate::compare::{CipCompareError, CompareContext};
use crate::digraph::CipDigraph;
use crate::ligand_profile::{Branch, collect_branches, sphere_profiles};
use crate::node::{CipNodeKind, NodeId};

pub(crate) fn eligible(mol: &Molecule, center: AtomIdx) -> bool {
    let direct: Vec<_> = mol.neighbors(center).map(|(atom, _)| atom).collect();
    if direct.len() < 3 {
        return false;
    }
    let mut visited = FxHashSet::default();
    visited.insert(center);
    let mut returning_ligands = 0;
    for &start in &direct {
        if visited.contains(&start) {
            continue;
        }
        let mut pending = vec![start];
        let mut ligands = 0;
        while let Some(atom) = pending.pop() {
            if !visited.insert(atom) {
                continue;
            }
            ligands += usize::from(direct.contains(&atom));
            for (neighbor, _) in mol.neighbors(atom) {
                if !visited.contains(&neighbor) {
                    pending.push(neighbor);
                }
            }
        }
        if ligands > 1 {
            returning_ligands += ligands;
        }
    }
    returning_ligands >= 3
}

pub(crate) fn has_reference_partition(groups: &[Vec<NodeId>]) -> bool {
    match groups {
        [only] => only.len() == 4,
        [first, second] => first.len() + second.len() == 4,
        _ => false,
    }
}

fn represented_atom(graph: &CipDigraph, node: NodeId) -> Option<AtomIdx> {
    match graph.node(node).kind {
        CipNodeKind::Atom { atom_idx } => Some(atom_idx),
        CipNodeKind::MultipleBondDuplicate {
            duplicated_atom, ..
        } => Some(duplicated_atom),
        CipNodeKind::RingDuplicate { closure_atom, .. } => Some(closure_atom),
        CipNodeKind::ImplicitHydrogen => None,
    }
}

fn refine_branch_order(
    graph: &mut CipDigraph,
    branches: &mut Branch,
    context: &mut CompareContext,
    cache: &mut HashMap<NodeId, Option<CipCode>>,
    reference: Option<CipCode>,
) -> Result<(), CipCompareError> {
    let mol = graph.molecule();
    // Fixed auxiliary descriptors feed earlier rules; promotion never
    // changes a descriptor used when establishing these sibling groups.
    for (_, groups) in branches
        .iter_mut()
        .filter(|(_, groups)| groups.iter().any(|group| group.len() > 1))
    {
        let mut ranked = crate::rule3::refine(graph, std::mem::take(groups), context)?;
        if let Some(reference) = reference {
            let positions: Vec<_> = ranked
                .iter()
                .flatten()
                .copied()
                .filter(|&node| {
                    !matches!(
                        graph.node(node).kind,
                        CipNodeKind::MultipleBondDuplicate { .. }
                    )
                })
                .collect();
            ranked = crate::rule4b::refine_stereo_groups_with_reference(
                graph, mol, ranked, &positions, cache, reference,
            )?
            .0;
        }
        *groups = ranked;
    }
    Ok(())
}

pub(crate) fn assign(
    graph: &mut CipDigraph,
    groups: &[Vec<NodeId>],
    positions: &[NodeId],
    chirality: Chirality,
    context: &mut CompareContext,
    cache: &mut HashMap<NodeId, Option<CipCode>>,
    preceding_order_flips: bool,
) -> Result<Option<CipCode>, CipCompareError> {
    let tied_starts = groups
        .iter()
        .filter(|group| group.len() > 1)
        .flatten()
        .copied();
    let mut branches = collect_branches(graph, tied_starts.clone(), context)?;
    if returning_paths(graph, &branches) < 3 {
        return Ok(None);
    }
    let refs: FxHashSet<_> = tied_starts
        .filter_map(|node| represented_atom(graph, node))
        .collect();
    let has_stereo = branches.iter().any(|&(node, _)| {
        matches!(graph.node(node).kind, CipNodeKind::Atom { atom_idx }
            if graph.molecule().atom(atom_idx).chirality.is_tetrahedral())
    });
    let reflected = has_stereo.then(|| branches.clone());
    if has_stereo || !graph.alkene_geometry().is_empty() {
        refine_branch_order(
            graph,
            &mut branches,
            context,
            cache,
            has_stereo.then_some(CipCode::R),
        )?;
    }
    let Some(is_r) = promoted_parity(graph, &branches, &refs, groups, positions, chirality) else {
        return Ok(None);
    };
    if let Some(mut reflected) = reflected {
        refine_branch_order(graph, &mut reflected, context, cache, Some(CipCode::S))?;
        let mirrored_order =
            promoted_parity(graph, &reflected, &refs, groups, positions, chirality);
        // Reversing the centre's tag must reverse an absolute descriptor.
        // A flip in ligand order would cancel that reversal.
        if mirrored_order != Some(is_r ^ preceding_order_flips) {
            return Ok(None);
        }
    }
    Ok(Some(if is_r { CipCode::R } else { CipCode::S }))
}

fn returning_paths(graph: &CipDigraph, branches: &Branch) -> usize {
    let Some(focus) = branches
        .first()
        .and_then(|&(node, _)| graph.node(node).parent)
    else {
        return 0;
    };
    let CipNodeKind::Atom { atom_idx } = graph.node(focus).kind else {
        return 0;
    };
    let mut incoming = FxHashSet::default();
    for &(mut node, _) in branches {
        if !matches!(graph.node(node).kind, CipNodeKind::RingDuplicate { closure_atom, .. } if closure_atom == atom_idx)
        {
            continue;
        }
        while let Some(parent) = graph.node(node).parent {
            if parent == focus {
                incoming.insert(node);
                break;
            }
            node = parent;
        }
    }
    incoming.len()
}

fn promoted_parity(
    graph: &CipDigraph,
    branches: &Branch,
    refs: &FxHashSet<AtomIdx>,
    groups: &[Vec<NodeId>],
    positions: &[NodeId],
    chirality: Chirality,
) -> Option<bool> {
    let mut parity_balance = 0;
    for &reference in refs {
        let profiles = sphere_profiles(branches, |node| {
            u8::from(represented_atom(graph, node) == Some(reference))
        });
        let mut ranked = Vec::with_capacity(4);
        for group in groups {
            let mut group = group.clone();
            if group.len() > 1 {
                group.sort_by(|a, b| profiles[b].cmp(&profiles[a]));
                if group
                    .windows(2)
                    .any(|pair| profiles[&pair[0]] == profiles[&pair[1]])
                {
                    break;
                }
            }
            ranked.extend(group.into_iter().map(|node| vec![node]));
        }
        if ranked.len() == 4
            && let Some(is_r) = resolve_is_r_from_groups(&ranked, positions, chirality)
        {
            parity_balance += if is_r { 1 } else { -1 };
        }
    }
    (parity_balance != 0).then_some(parity_balance > 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CipBudget, SkipReason, assign_cip_accurate_experimental};

    #[test]
    fn authored_reference_ligand_labels_survive_new_traversals() {
        // Hanson et al., ValidationSuite. The set includes 2+2, 3+1
        // and 4-way ties, with both members of enantiomer pairs.
        for (text, want) in [
            ("C1=C2[C@@]34C(=C1)CCCCC3=CC=C4CCCC2", CipCode::S), // VS281
            ("C1=C2[C@]34C(=C1)CCCCC3=CC=C4CCCC2", CipCode::R),  // VS282
            ("C1=C2[C@@H]3C(CCC=C3CC1)=CCC2", CipCode::S),       // VS283
            ("C1=C2[C@H]3C(CCC=C3CC1)=CCC2", CipCode::R),        // VS284
            ("O=C1[C@]2(NCC1)C(CCN2)=O", CipCode::R),            // VS285
            ("O=C1[C@@]2(CCCC2=O)CCC1", CipCode::S),             // VS288
            ("C1=C2[C@@H]3C(C1)=CCC3=CC2", CipCode::R),          // VS289
            ("C1=C2[C@H]3C(C1)=CCC3=CC2", CipCode::S),           // VS290
            ("C1[C@@]2(OCC1)OCCC2", CipCode::R),                 // VS291
            ("C1=C2[C@@]34C(=C1)CCC3=CC=C4CC2", CipCode::S),     // VS292
        ] {
            let mol = chematic_smiles::parse(text).unwrap();
            let variants = std::iter::once(text.to_owned())
                .chain((1..=16).map(|seed| chematic_smiles::random_smiles(&mol, seed)));
            for variant in variants {
                let input = chematic_smiles::parse(&variant).unwrap();
                let assignments =
                    assign_cip_accurate_experimental(&input, CipBudget::default_budget())
                        .unwrap()
                        .assignments;
                assert_eq!(assignments.len(), 1, "{text} -> {variant}");
                assert_eq!(assignments[0].1, want, "{text} -> {variant}");
            }
        }
    }

    #[test]
    fn symmetry_cancels_reference_assignments() {
        for text in [
            "C1=C2CC=C3[C@]24C(CC=C4C1)=CC3", // VS012, S4 symmetry
            "[C@@H]12C3C4C1C5C4C3C25",        // singly tagged cubane
            "[C@](Cl)(Cl)(Cl)Br",             // repeated acyclic ligands
        ] {
            let mol = chematic_smiles::parse(text).unwrap();
            let result =
                assign_cip_accurate_experimental(&mol, CipBudget::default_budget()).unwrap();
            assert!(result.assignments.is_empty(), "{text}");
            assert_eq!(result.skipped.len(), 1, "{text}");
            assert_eq!(result.skipped[0].1, SkipReason::Tied, "{text}");
        }
    }

    #[test]
    fn rule6_uses_auxiliary_descriptors_after_earlier_rules() {
        let cases: &[(&str, &[(u32, CipCode)])] = &[
            (
                r"C[C@H]([C@H]([C@@]([C@@H]([C@@]([C@@H]([C@H](C)Cl)[C@@H](Cl)C)([C@H]([C@H](C)Cl)[C@@H](C)Cl)[C@@H]([C@H](C)Cl)[C@@H](C)Cl)O)([C@H]([C@H](C)Cl)[C@H](C)Cl)[C@@H]([C@@H](Cl)C)[C@@H](C)Cl)[C@H](C)Cl)Cl",
                &[(3, CipCode::R), (5, CipCode::R)],
            ), // VS273
            (
                r"C[C@H]([C@H]([C@@]([C@@H]([C@@]([C@@H]([C@H](C)Cl)[C@@H](C)Cl)([C@H]([C@@H](C)Cl)[C@@H](C)Cl)[C@@H]([C@@H](C)Cl)[C@H](C)Cl)O)([C@H]([C@@H](C)Cl)[C@H](C)Cl)[C@@H]([C@H](C)Cl)[C@H](C)Cl)[C@@H](Cl)C)Cl",
                &[(3, CipCode::S), (5, CipCode::R)],
            ), // VS274
            (
                r"C1C[C@]23CCC[C@@]24CCC[C@@]45CCC[C@@]56CCC[C@]36C1",
                &[
                    (2, CipCode::R),
                    (6, CipCode::R),
                    (10, CipCode::R),
                    (14, CipCode::R),
                    (18, CipCode::R),
                ],
            ), // VS280
            (
                r"C1[C@@]2(CC[C@]13CCCC3)C[C@@]4(CC2)CCCC4",
                &[(1, CipCode::S)],
            ), // VS286
            (
                r"C\1C/C/2=C/CC/C/3=C/CC/C(=C1)/C[C@@H](C2)C3",
                &[(13, CipCode::S)],
            ), // VS293
            (
                r"O[C@H]1O[C@@]2(O[C@H]1O)O[C@H]([C@H](O2)O)O",
                &[(3, CipCode::R)],
            ), // VS294
            (
                r"O[C@H]1O[C@@]2(O[C@H]1O)O[C@@H]([C@@H](O2)O)O",
                &[(3, CipCode::S)],
            ), // VS295
            (
                r"C1[C@@]2(OCCC1)CC[C@]3(CC2)CC[C@@]4(CC3)OCCCC4",
                &[(8, CipCode::S)],
            ), // VS296
            (
                r"C1[C@@]2(CC[C@]3(C1)OCCC3)CC[C@@]4(CC2)CCCO4",
                &[(1, CipCode::R)],
            ), // VS297
            (r"Cl[C@H]1C[C@]2(C1)C[C@H](C2)Cl", &[(3, CipCode::S)]), // VS298
            (r"C\C=C/1\C[C@]2(C\C(\C2)=C\C)C1", &[(4, CipCode::S)]), // VS299
            (
                r"O[C@H](C[C@@]12C=3CCC=C2CCC=C1CCC3)C[C@]45C=6CCC=C5CCC=C4CCC6",
                &[(1, CipCode::LowerS), (3, CipCode::R), (17, CipCode::S)],
            ), // VS300
        ];
        for &(text, expected) in cases {
            let molecule = chematic_smiles::parse(text).unwrap();
            let n = molecule.atom_count();
            for offset in 0..n.min(8) {
                let order: Vec<_> = (0..n)
                    .map(|i| AtomIdx(((n + offset - i) % n) as u32))
                    .collect();
                let reordered = molecule.renumber_atoms(&order).unwrap();
                let labels =
                    assign_cip_accurate_experimental(&reordered, CipBudget::default_budget())
                        .unwrap();
                for &(idx, code) in expected {
                    assert!(
                        labels
                            .assignments
                            .iter()
                            .any(
                                |&(mapped, seen)| order[mapped.0 as usize].0 == idx && seen == code
                            ),
                        "{text}, offset {offset}, atom {idx}, expected {code:?}: {labels:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn reflection_invariant_reference_order_does_not_assign_absolute_cage_stereo() {
        for text in [
            "[C@]12(C)C[C@H]3C[C@H](C[C@@H](C1)C3)C2",
            "[C@@]12(C)C[C@@H]3C[C@@H](C[C@H](C1)C3)C2",
        ] {
            let molecule = chematic_smiles::parse(text).unwrap();
            let result =
                assign_cip_accurate_experimental(&molecule, CipBudget::default_budget()).unwrap();
            assert!(
                result.skipped.contains(&(AtomIdx(0), SkipReason::Tied)),
                "{text}: {result:?}"
            );
            assert!(
                result
                    .assignments
                    .iter()
                    .all(|(_, code)| *code == CipCode::LowerS),
                "{text}: {result:?}"
            );
        }
    }
}
