//! Validate MCS tetrahedral stereo once the ligand mapping is known.

use chematic_core::{AtomIdx, Molecule, STEREO_H_SENTINEL};

use super::PartialMapping;
use crate::match_chirality::permutation_parity;

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
