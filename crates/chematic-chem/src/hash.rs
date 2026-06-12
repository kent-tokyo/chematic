//! Fast molecular structure hashing for duplicate detection.
//!
//! Provides structural fingerprinting and identity checking based on canonical
//! SMILES representation. Useful for detecting duplicate molecules in bulk operations.

use chematic_core::Molecule;
use chematic_smiles::canonical_smiles;

pub(crate) const FNV1A_OFFSET: u64 = 0xcbf29ce484222325;
pub(crate) const FNV1A_PRIME: u64 = 0x100000001b3;

/// FNV-1a 64-bit hash of a byte slice.
pub(crate) fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h = FNV1A_OFFSET;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(FNV1A_PRIME);
    }
    h
}

/// Compute a structural hash for a molecule using FNV-1a on the canonical SMILES.
///
/// The hash is deterministic and invariant to atom ordering. Identical molecules
/// always produce the same hash, but hash collisions are possible (use [`are_identical`]
/// for a definitive check).
///
/// Useful for rapid duplicate detection in large sets; pair with [`are_identical`]
/// for validation.
pub fn mol_hash(mol: &Molecule) -> u64 {
    fnv1a(canonical_smiles(mol).as_bytes())
}

/// Check whether two molecules are structurally identical.
///
/// Returns `true` if both molecules have the same canonical SMILES representation.
/// This is a definitive check (no collisions), unlike [`mol_hash`].
pub fn are_identical(a: &Molecule, b: &Molecule) -> bool {
    canonical_smiles(a) == canonical_smiles(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn mol_hash_same_smiles_same_hash() {
        let a = mol("CC");
        let b = mol("CC");
        assert_eq!(
            mol_hash(&a),
            mol_hash(&b),
            "same SMILES should give same hash"
        );
    }

    #[test]
    fn mol_hash_different_smiles_likely_different_hash() {
        let ethane = mol("CC");
        let propane = mol("CCC");
        // Not guaranteed, but extremely likely to differ
        assert_ne!(
            mol_hash(&ethane),
            mol_hash(&propane),
            "different SMILES should (very likely) differ"
        );
    }

    #[test]
    fn mol_hash_deterministic() {
        // Same molecule should hash to same value on repeated calls
        let m = mol("c1ccccc1");
        let hash1 = mol_hash(&m);
        let hash2 = mol_hash(&m);
        assert_eq!(hash1, hash2, "same molecule should produce same hash");
    }

    #[test]
    fn are_identical_true_for_same_structure() {
        let a = mol("CC");
        let b = mol("CC");
        assert!(
            are_identical(&a, &b),
            "identical structures should return true"
        );
    }

    #[test]
    fn are_identical_false_for_different() {
        let a = mol("CC");
        let b = mol("CCC");
        assert!(
            !are_identical(&a, &b),
            "different structures should return false"
        );
    }

    #[test]
    fn are_identical_aspirin() {
        let aspirin1 = mol("CC(=O)Oc1ccccc1C(=O)O");
        let aspirin2 = mol("O=C(c1ccccc1OC(=O)C)O");
        assert!(
            are_identical(&aspirin1, &aspirin2),
            "different notations of aspirin should be identical"
        );
    }
}
