//! Validate MCS tetrahedral stereo once the ligand mapping is known.

use chematic_core::{AtomIdx, Chirality, Molecule, STEREO_H_SENTINEL, StereoGroupKind};

use super::PartialMapping;
use crate::match_chirality::permutation_parity;
use crate::query::{AtomPrimitive, AtomQuery, QueryMolecule};

pub(super) fn mapped_chirality_compatible(mols: &[&Molecule], mapping: &PartialMapping) -> bool {
    mapping.query_to_mol.iter().all(|row| {
        let source_atom = row[0];
        let source_chirality = mols[0].atom(source_atom).chirality;
        if !source_chirality.is_tetrahedral() {
            return true;
        }
        mols.iter().enumerate().skip(1).all(|(mi, target)| {
            let target_atom = row[mi];
            let target_chirality = target.atom(target_atom).chirality;
            if !target_chirality.is_tetrahedral() {
                return false;
            }

            // Allene @AL1/@AL2 shares the atom tag variants with tetrahedral
            // stereo, but its six-entry terminal order needs a separate parity
            // rule. Keep the previous raw-tag policy for that class.
            if mols[0].is_allene_stereo(source_atom) || target.is_allene_stereo(target_atom) {
                return mols[0].is_allene_stereo(source_atom)
                    == target.is_allene_stereo(target_atom)
                    && source_chirality == target_chirality;
            }

            let (Some(source_order), Some(target_order)) = (
                mols[0].stereo_neighbor_order(source_atom),
                target.stereo_neighbor_order(target_atom),
            ) else {
                return source_chirality == target_chirality;
            };
            if source_order.len() != 4 || target_order.len() != 4 {
                return source_chirality == target_chirality;
            }

            // The stereocenter may be in a valid smaller MCS without all of
            // its substituents. Its configuration is not constrained until
            // both complete ligand sets participate in the mapping.
            if !target_order.iter().all(|&neighbor| {
                neighbor == STEREO_H_SENTINEL
                    || mapping.query_idx_of(mi, AtomIdx(neighbor)).is_some()
            }) {
                return true;
            }
            let mut projected_order = Vec::with_capacity(4);
            for &neighbor in source_order {
                if neighbor == STEREO_H_SENTINEL {
                    projected_order.push(STEREO_H_SENTINEL);
                } else if let Some(query_idx) = mapping.query_idx_of(0, AtomIdx(neighbor)) {
                    projected_order.push(mapping.query_to_mol[query_idx][mi].0);
                } else {
                    return true;
                }
            }
            permutation_parity(&projected_order, target_order)
                .is_some_and(|even| (source_chirality == target_chirality) == even)
        })
    })
}

pub(super) fn constrain_query_chirality(
    query: &mut QueryMolecule,
    mols: &[&Molecule],
    mapped: impl Fn(usize, usize) -> AtomIdx,
) {
    for (q, neighbors) in query.adj.iter().enumerate() {
        let tag = mols[0].atom(mapped(q, 0)).chirality;
        if tag == Chirality::None {
            continue;
        }

        let source_order = mols[0].stereo_neighbor_order(mapped(q, 0));
        let has_bracket_h = source_order.is_some_and(|order| order.contains(&STEREO_H_SENTINEL));
        let hydrogen_position = usize::from(neighbors.first().is_some_and(|&(_, n)| n < q));
        let mut common_kind = None;
        for (m, mol) in mols.iter().enumerate() {
            let center = mapped(q, m);
            let projections = neighbors.iter().map(|&(_, n)| mapped(n, m).0);
            let kind = projected_tetrahedral_kind(
                mol,
                center,
                projections,
                has_bracket_h,
                hydrogen_position,
            );
            if kind.is_none() || common_kind.is_some_and(|first| Some(first) != kind) {
                common_kind = None;
                break;
            }
            common_kind = kind;
        }

        let atom = &mut query.atoms[q].query;
        let constraint = |predicate| AtomQuery::Primitive(predicate);
        let and = |left, right| AtomQuery::And(Box::new(left), Box::new(right));
        // Kind 3 can check presence in the native matcher. It does not
        // represent partial, non-tetrahedral or enhanced stereo as SMARTS.
        *atom = and(
            atom.clone(),
            constraint(AtomPrimitive::Chirality(common_kind.unwrap_or(3))),
        );
        if common_kind.is_some() && has_bracket_h {
            *atom = and(atom.clone(), constraint(AtomPrimitive::HCount(1)));
        }
    }
}

fn projected_tetrahedral_kind(
    mol: &Molecule,
    center: AtomIdx,
    neighbors: impl Iterator<Item = u32>,
    has_bracket_h: bool,
    hydrogen_position: usize,
) -> Option<u8> {
    let atom = mol.atom(center);
    if !atom.chirality.is_tetrahedral()
        || mol.is_allene_stereo(center)
        || mol.stereo_groups().iter().any(|group| {
            group.kind != StereoGroupKind::Absolute && group.atom_indices.contains(&center)
        })
    {
        return None;
    }
    let order = mol.stereo_neighbor_order(center)?;
    if order.len() != 4
        || has_bracket_h != order.contains(&STEREO_H_SENTINEL)
        || (has_bracket_h && atom.hydrogen_count != Some(1))
    {
        return None;
    }
    let mut projected: Vec<_> = neighbors.collect();
    if has_bracket_h && projected.len() == 3 {
        projected.insert(hydrogen_position, STEREO_H_SENTINEL);
    }
    let even = permutation_parity(order, &projected)?;
    let clockwise = (atom.chirality == Chirality::Clockwise) == even;
    Some(if clockwise { 2 } else { 1 })
}
