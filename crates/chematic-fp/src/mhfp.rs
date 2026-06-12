//! C-Series Phase 1: MHFP/SECFP MinHash fingerprints.
//!
//! **NOTE**: Current implementation is a simplified MinHash approximation.
//! Uses ECFP4-derived bit positions with DefaultHasher.
//! Full MHFP (Lowe & Sayle 2013) would extract circular substructure SMILES
//! and compute MinHash on those strings for improved accuracy.
//! See: https://pubs.acs.org/doi/10.1021/ci034236b
//!
//! Useful for:
//! - Fast approximate similarity searches (lower accuracy than true MHFP)
//! - Large molecular library screening
//! - Hash-based clustering
//! - Scalable molecular database queries
//!
//! TODO (v0.1.90+): Upgrade to true MHFP by:
//! 1. Extract circular substructure SMILES (radius 0-4 like ECFP)
//! 2. Hash each SMILES string individually
//! 3. Return K smallest hashes across all radii

use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;

use crate::ecfp::ecfp4;

/// Configuration for MinHash fingerprint generation.
#[derive(Clone, Debug)]
pub struct MhfpConfig {
    /// Number of hash functions (lanes) for MinHash
    pub num_hashes: usize,
    /// Seed offset for hash functions
    pub seed: u64,
}

impl Default for MhfpConfig {
    fn default() -> Self {
        MhfpConfig {
            num_hashes: 128,
            seed: 0,
        }
    }
}

/// MinHash fingerprint for approximate similarity searching.
#[derive(Clone, Debug)]
pub struct MhfpFingerprint {
    /// Hash values (one per hash function / lane)
    pub hashes: Vec<u64>,
    /// Number of hash functions used
    pub num_hashes: usize,
}

impl MhfpFingerprint {
    /// Calculate Tanimoto-like similarity between MinHash fingerprints.
    ///
    /// Uses the Jaccard index approximation from MinHash theory:
    /// Jaccard(A, B) ≈ (agreements / total_hashes)
    pub fn tanimoto(&self, other: &MhfpFingerprint) -> f64 {
        if self.num_hashes == 0 {
            return 1.0;
        }

        let mut agreements = 0;
        for (h1, h2) in self.hashes.iter().zip(other.hashes.iter()) {
            if h1 == h2 {
                agreements += 1;
            }
        }

        agreements as f64 / self.num_hashes as f64
    }
}

/// Generate MinHash fingerprint from a molecule (simplified approximation).
///
/// **Current Implementation**: Uses ECFP4-derived bit positions with DefaultHasher.
/// This is a simplified approximation for fast similarity search.
///
/// **True MHFP** (Lowe & Sayle 2013) would:
/// 1. Extract circular substructure SMILES (radius 0-4)
/// 2. Hash each SMILES string individually
/// 3. Retain K smallest hashes
/// 4. Achieve higher accuracy for similarity searching
///
/// For production systems requiring high accuracy, consider ECFP4 Tanimoto
/// instead, or implement full MHFP per reference paper.
pub fn mhfp(mol: &chematic_core::Molecule) -> MhfpFingerprint {
    mhfp_with_config(mol, &MhfpConfig::default())
}

/// Generate MinHash fingerprint with custom configuration.
pub fn mhfp_with_config(
    mol: &chematic_core::Molecule,
    config: &MhfpConfig,
) -> MhfpFingerprint {
    // Get ECFP4 bitvector
    let ecfp = ecfp4(mol);

    // Extract set of bit positions from ECFP
    let mut bit_set = Vec::new();
    for i in 0..2048 {
        if ecfp.get(i) {
            bit_set.push(i as u64);
        }
    }

    // If no bits set, return all max values
    if bit_set.is_empty() {
        return MhfpFingerprint {
            hashes: vec![u64::MAX; config.num_hashes],
            num_hashes: config.num_hashes,
        };
    }

    // Compute MinHash values
    let mut hashes = vec![u64::MAX; config.num_hashes];

    for h in 0..config.num_hashes {
        let seed = config.seed.wrapping_add(h as u64);

        for &bit_pos in &bit_set {
            let mut hasher = DefaultHasher::new();
            hasher.write_u64(seed);
            hasher.write_u64(bit_pos);
            let hash_val = hasher.finish();

            if hash_val < hashes[h] {
                hashes[h] = hash_val;
            }
        }
    }

    MhfpFingerprint {
        hashes,
        num_hashes: config.num_hashes,
    }
}

/// Convenience function for standard MHFP generation.
pub fn mhfp_128(mol: &chematic_core::Molecule) -> MhfpFingerprint {
    mhfp(mol)
}

/// Calculate Tanimoto-like similarity between two molecules using MHFP.
pub fn tanimoto_mhfp(mol1: &chematic_core::Molecule, mol2: &chematic_core::Molecule) -> f64 {
    let fp1 = mhfp(mol1);
    let fp2 = mhfp(mol2);
    fp1.tanimoto(&fp2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_mhfp_simple() {
        let mol = parse("CC").unwrap();
        let fp = mhfp(&mol);

        assert_eq!(fp.num_hashes, 128);
        assert_eq!(fp.hashes.len(), 128);
    }

    #[test]
    fn test_mhfp_identical() {
        let mol = parse("CC").unwrap();
        let fp1 = mhfp(&mol);
        let fp2 = mhfp(&mol);

        // Identical molecules should have identical hashes
        assert_eq!(fp1.hashes, fp2.hashes);
        assert!((fp1.tanimoto(&fp2) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_mhfp_different_molecules() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("CCC").unwrap();

        let fp1 = mhfp(&mol1);
        let fp2 = mhfp(&mol2);

        // Different molecules should have different hashes (usually)
        let similarity = fp1.tanimoto(&fp2);
        assert!(similarity >= 0.0 && similarity <= 1.0);
    }

    #[test]
    fn test_mhfp_symmetry() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("CCC").unwrap();

        let sim12 = tanimoto_mhfp(&mol1, &mol2);
        let sim21 = tanimoto_mhfp(&mol2, &mol1);

        // Similarity should be symmetric
        assert!((sim12 - sim21).abs() < 1e-10);
    }

    #[test]
    fn test_mhfp_config() {
        let mol = parse("CC").unwrap();
        let config = MhfpConfig {
            num_hashes: 64,
            seed: 42,
        };

        let fp = mhfp_with_config(&mol, &config);
        assert_eq!(fp.num_hashes, 64);
        assert_eq!(fp.hashes.len(), 64);
    }

    #[test]
    fn test_mhfp_similar_molecules() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("CC").unwrap();

        let fp1 = mhfp(&mol1);
        let fp2 = mhfp(&mol2);

        let similarity = fp1.tanimoto(&fp2);
        // Identical SMILES should have high similarity
        assert!(similarity > 0.9);
    }

    #[test]
    fn test_mhfp_empty_molecule() {
        let mol = parse("C").unwrap();
        let fp = mhfp(&mol);

        // Single atom still produces fingerprint
        assert_eq!(fp.num_hashes, 128);
    }

    #[test]
    fn test_mhfp_bounds() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("CCCCCCCC").unwrap();

        let fp1 = mhfp(&mol1);
        let fp2 = mhfp(&mol2);

        let similarity = fp1.tanimoto(&fp2);
        // Jaccard-like similarity should be in [0, 1]
        assert!(similarity >= 0.0 && similarity <= 1.0);
    }
}
