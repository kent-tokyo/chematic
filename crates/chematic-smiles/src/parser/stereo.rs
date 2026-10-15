//! Parse-time ligand positions for stereocenters, including vacant sites.

use chematic_core::{Atom, AtomIdx, Chirality, Molecule, STEREO_VACANT_SENTINEL};

use super::{Parser, StereoEntry};

/// Retain the vacant-site marker only when exactly one coordination site is
/// unoccupied. Fully coordinated and more incomplete centers keep their
/// prior ligand-order representation.
pub(super) fn retain_single_vacant_site(molecule: &mut Molecule, tagged_atoms: &[AtomIdx]) {
    for &center in tagged_atoms {
        let Some((class, _)) = molecule.atom(center).chirality.nontetrahedral() else {
            continue;
        };
        let Some(order) = molecule.stereo_neighbor_order(center) else {
            continue;
        };
        if order.len() != class.max_nbors() {
            let order = order
                .iter()
                .copied()
                .filter(|&id| id != STEREO_VACANT_SENTINEL)
                .collect();
            molecule.set_stereo_neighbor_order(center, order);
        }
    }
}

impl Parser<'_> {
    /// Start a stereo record for `atom_idx` if `atom` has chirality.
    pub(super) fn begin_stereo_if_chiral(
        atom: &Atom,
        atom_idx: AtomIdx,
        from_atom: Option<AtomIdx>,
        current_stereo: &mut Option<(AtomIdx, Vec<StereoEntry>)>,
    ) {
        if atom.chirality == Chirality::None {
            return;
        }
        let mut entries: Vec<StereoEntry> = Vec::new();
        if let Some(prev) = from_atom {
            entries.push(StereoEntry::Atom(prev));
        }
        if atom.chirality.nontetrahedral().is_some() {
            entries.push(StereoEntry::Vacant);
        }
        if atom.hydrogen_count.is_some_and(|h| h > 0) {
            entries.push(StereoEntry::ImplicitH);
        }
        *current_stereo = Some((atom_idx, entries));
    }
}
