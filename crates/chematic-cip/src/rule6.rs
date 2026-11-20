//! Reference-ligand ranking proposed by Hanson et al. (2018).
//!
//! This pass handles tetrahedral configurations whose connected component
//! has no other stereo descriptors to consider. Auxiliary Rule 6 and the
//! interaction with earlier stereochemical rules require a rooted resolver.

use chematic_core::{AtomIdx, BondOrder, Chirality, CipCode, Molecule};
use rustc_hash::FxHashSet;

use crate::assign::resolve_is_r_from_groups;
use crate::compare::{CipCompareError, CompareContext};
use crate::digraph::CipDigraph;
use crate::ligand_profile::{collect_branches, sphere_profiles};
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
            if mol.atom(atom).chirality != Chirality::None {
                return false;
            }
            ligands += usize::from(direct.contains(&atom));
            for (neighbor, bond) in mol.neighbors(atom) {
                if mol.bond_direction(bond).is_some()
                    || matches!(mol.bond(bond).order, BondOrder::Up | BondOrder::Down)
                {
                    // Sequence Rule 3 may be decisive at this centre.
                    return false;
                }
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

pub(crate) fn assign(
    graph: &mut CipDigraph,
    groups: &[Vec<NodeId>],
    positions: &[NodeId],
    chirality: Chirality,
    context: &mut CompareContext,
) -> Result<Option<CipCode>, CipCompareError> {
    let tied_starts = groups
        .iter()
        .filter(|group| group.len() > 1)
        .flatten()
        .copied();
    let branches = collect_branches(graph, tied_starts.clone(), context)?;
    let refs: FxHashSet<_> = tied_starts
        .filter_map(|node| represented_atom(graph, node))
        .collect();
    let mut parity_balance = 0;
    for reference in refs {
        let profiles = sphere_profiles(&branches, |node| {
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
    Ok(match parity_balance {
        value if value > 0 => Some(CipCode::R),
        value if value < 0 => Some(CipCode::S),
        _ => None,
    })
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
}
