//! `chematic-perception` — molecular perception algorithms.
//!
//! Provides:
//! - [`sssr`]: Smallest Set of Smallest Rings (SSSR) via Balducci-Pearlman algorithm.
//! - [`aromaticity`]: Hückel aromaticity perception for kekulized molecules.

#![forbid(unsafe_code)]

pub mod aromaticity;
pub mod cip_priority;
pub mod pharmacophore;
pub mod ring_family;
pub mod sssr;
pub mod stereo_validation;

pub mod stereo2d;

pub use aromaticity::{
    AromaticityAlgorithm, AromaticityModel, AtomElectronTrace, ContributionReason, RingAromaticity,
    RingElectronTrace, all_ring_list, apply_aromaticity, apply_aromaticity_ex, aromatic_ring_list,
    assign_aromaticity, assign_aromaticity_ex, augmented_ring_set, count_aromatic_rings,
    ring_bonds_all_aromatic, trace_ring_pi_electrons,
};
pub use chematic_core::{ValenceError, validate_valence};
pub use pharmacophore::{Feature, FeatureType, detect_features, features_to_bitvec};
pub use ring_family::{RingFamily, RingSystemKind, find_ring_families, find_ring_families_over};
pub use sssr::{RingSet, find_sssr};
pub use stereo_validation::{
    StereoCompleteness, StereoError, StereoErrorKind, stereo_completeness, validate_stereo,
};
pub use stereo2d::{
    StereoAssignment2D, apply_stereo_from_2d, assign_ez_from_2d, assign_stereo_from_2d,
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

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(smiles: &str) -> Molecule {
        parse(smiles).expect("valid SMILES")
    }

    #[test]
    fn test_ring_membership_benzene() {
        let m = mol("c1ccccc1");
        let membership = ring_membership(&m);
        assert_eq!(membership.len(), 6);
        for atom_membership in membership.iter().take(6) {
            assert_eq!(
                atom_membership.len(),
                1,
                "each benzene atom in exactly 1 ring"
            );
            assert_eq!(atom_membership[0], 0, "all in ring index 0");
        }
    }

    #[test]
    fn test_ring_membership_naphthalene() {
        let m = mol("c1ccc2ccccc2c1");
        let membership = ring_membership(&m);
        assert_eq!(membership.len(), 10);
        // In naphthalene SSSR, some atoms appear in 2 rings depending on the ring decomposition
        // Just verify all atoms are in at least 1 ring
        for mem in &membership {
            assert!(
                !mem.is_empty(),
                "all naphthalene atoms should be in at least 1 ring"
            );
        }
    }

    #[test]
    fn test_ring_membership_acyclic() {
        let m = mol("CC");
        let membership = ring_membership(&m);
        assert_eq!(membership.len(), 2);
        for mem in &membership {
            assert!(mem.is_empty(), "ethane atoms should not be in rings");
        }
    }

    #[test]
    fn test_ring_sizes_for_atom_benzene() {
        let m = mol("c1ccccc1");
        let sizes = ring_sizes_for_atom(&m, 0);
        assert_eq!(sizes, vec![6]);
    }

    #[test]
    fn test_ring_sizes_for_atom_naphthalene() {
        let m = mol("c1ccc2ccccc2c1");
        // Just verify naphthalene atoms are in rings of size 6
        let sizes = ring_sizes_for_atom(&m, 0);
        assert!(!sizes.is_empty());
        assert!(sizes.contains(&6), "naphthalene has 6-membered rings");
    }

    #[test]
    fn test_ring_sizes_for_atom_acyclic() {
        let m = mol("CC");
        let sizes = ring_sizes_for_atom(&m, 0);
        assert!(sizes.is_empty());
    }

    #[test]
    fn test_is_fused_ring_naphthalene() {
        let m = mol("c1ccc2ccccc2c1");
        assert!(is_fused_ring_system(&m), "naphthalene is fused");
    }

    #[test]
    fn test_is_fused_ring_benzene() {
        let m = mol("c1ccccc1");
        assert!(
            !is_fused_ring_system(&m),
            "single benzene ring is not fused"
        );
    }

    #[test]
    fn test_is_fused_ring_spiro() {
        // Spiro[4.4]nonane has two rings sharing only 1 atom
        let m = mol("C1CCC2(C1)CCCC2");
        assert!(
            !is_fused_ring_system(&m),
            "spiro compound shares only 1 atom, not fused"
        );
    }

    #[test]
    fn test_aromatize_benzene() {
        let mut m = mol("c1ccccc1");
        aromatize(&mut m);
        for (_, atom) in m.atoms() {
            assert!(atom.aromatic, "all benzene atoms should be aromatic");
        }
        for (_, bond) in m.bonds() {
            assert_eq!(
                bond.order,
                chematic_core::BondOrder::Aromatic,
                "all benzene bonds should be aromatic"
            );
        }
    }

    #[test]
    fn test_kekulize_inplace_benzene() {
        let mut m = mol("c1ccccc1");
        kekulize_inplace(&mut m).expect("benzene should kekulize");
        let mut single_count = 0;
        let mut double_count = 0;
        for (_, bond) in m.bonds() {
            match bond.order {
                chematic_core::BondOrder::Single => single_count += 1,
                chematic_core::BondOrder::Double => double_count += 1,
                _ => panic!("unexpected bond order after kekulization"),
            }
        }
        assert_eq!(single_count, 3);
        assert_eq!(double_count, 3);
    }
}
