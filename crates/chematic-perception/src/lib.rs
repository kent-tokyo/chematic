//! `chematic-perception` — molecular perception algorithms.
//!
//! Provides:
//! - [`sssr`]: Smallest Set of Smallest Rings (SSSR) via Balducci-Pearlman algorithm.
//! - [`aromaticity`]: Hückel aromaticity perception for kekulized molecules.

#![forbid(unsafe_code)]

pub mod aromaticity;
pub mod stereo_validation;
pub mod sssr;

pub mod stereo2d;

pub use aromaticity::{AromaticityModel, RingAromaticity, apply_aromaticity, assign_aromaticity};
pub use chematic_core::{ValenceError, validate_valence};
pub use sssr::{RingSet, find_sssr};
pub use stereo_validation::{StereoCompleteness, StereoError, StereoErrorKind, stereo_completeness, validate_stereo};
pub use stereo2d::{
    StereoAssignment2D,
    apply_stereo_from_2d,
    assign_stereo_from_2d,
    assign_ez_from_2d,
    cip_ez_descriptor,
};

use chematic_core::{AtomIdx, Molecule};

// ---------------------------------------------------------------------------
// Ring system helper API
// ---------------------------------------------------------------------------

/// For each atom, return the list of SSSR ring indices that contain it.
///
/// The outer `Vec` is indexed by atom position; each inner `Vec` contains
/// 0-based indices into the SSSR ring list (`find_sssr(mol).rings()`).
/// Atoms that belong to no ring get an empty inner vec.
pub fn ring_membership(mol: &Molecule) -> Vec<Vec<usize>> {
    let ring_set = find_sssr(mol);
    let rings = ring_set.rings();
    let n = mol.atom_count();
    let mut membership: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (ring_idx, ring) in rings.iter().enumerate() {
        for &atom in ring {
            membership[atom.0 as usize].push(ring_idx);
        }
    }
    membership
}

/// Return the sizes of all SSSR rings that contain `atom_idx`.
///
/// Returns an empty vec for acyclic atoms.
pub fn ring_sizes_for_atom(mol: &Molecule, atom_idx: usize) -> Vec<usize> {
    let ring_set = find_sssr(mol);
    let target = AtomIdx(atom_idx as u32);
    ring_set
        .rings()
        .iter()
        .filter(|ring| ring.contains(&target))
        .map(|ring| ring.len())
        .collect()
}

/// Return `true` if the molecule contains a fused ring system.
///
/// Two rings are fused when they share at least one bond (i.e. two adjacent
/// atoms in both rings).  Spiro rings (sharing exactly one atom) return `false`.
pub fn is_fused_ring_system(mol: &Molecule) -> bool {
    let ring_set = find_sssr(mol);
    let rings = ring_set.rings();
    for i in 0..rings.len() {
        for j in (i + 1)..rings.len() {
            // Count shared atoms.
            let shared = rings[i].iter().filter(|a| rings[j].contains(a)).count();
            if shared >= 2 {
                return true; // two rings share an edge → fused
            }
        }
    }
    false
}

/// Apply aromaticity to `mol` in-place (wrapper for [`apply_aromaticity`]).
pub fn aromatize(mol: &mut Molecule) {
    *mol = apply_aromaticity(mol);
}

/// Convert `mol` to Kekulé form in-place (wrapper for `kekulize` + `apply_kekule`).
///
/// Returns `Err` if kekulization fails (e.g. invalid aromatic system).
pub fn kekulize_inplace(mol: &mut Molecule) -> Result<(), chematic_core::KekuleError> {
    use chematic_core::{apply_kekule, kekulize};
    let result = kekulize(mol)?;
    *mol = apply_kekule(mol, &result);
    Ok(())
}
