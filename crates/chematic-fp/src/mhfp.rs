//! v0.1.91: True MHFP (MinHash Fingerprint) — Circular Fragment MinHash
//!
//! Implements MinHash algorithm based on Lowe & Sayle 2013:
//! 1. Extract circular fragments (radius 0-4, centered on each atom)
//! 2. Hash each fragment's structural signature
//! 3. Retain K smallest hashes to form the fingerprint
//!
//! Structural signature includes atomic properties and bond connectivity
//! without requiring SMILES generation (to avoid dev-dep on chematic_smiles).
//!
//! Reference: https://pubs.acs.org/doi/10.1021/ci034236b

use std::collections::HashSet;
use chematic_core::{AtomIdx, Molecule, BondOrder};

use crate::ecfp::fnv1a as fnv1a_hash;

/// Configuration for MinHash fingerprint generation.
#[derive(Clone, Debug)]
pub struct MhfpConfig {
    /// Number of hash functions (lanes) for MinHash
    pub num_hashes: usize,
    /// Seed offset for hash functions
    pub seed: u64,
    /// ECFP radius for circular fragment extraction (default 2 = ECFP4)
    pub radius: u32,
}

impl Default for MhfpConfig {
    fn default() -> Self {
        MhfpConfig {
            num_hashes: 128,
            seed: 0,
            radius: 2,
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

/// Collect all atoms within a given radius from a center atom (BFS).
fn atoms_within_radius(mol: &Molecule, center: AtomIdx, radius: u32) -> Vec<AtomIdx> {
    let mut visited = HashSet::new();
    let mut frontier = vec![center];
    visited.insert(center);

    for _ in 0..radius {
        let mut next_frontier = vec![];
        for atom in frontier {
            for (nb, _) in mol.neighbors(atom) {
                if visited.insert(nb) {
                    next_frontier.push(nb);
                }
            }
        }
        frontier = next_frontier;
    }

    visited.into_iter().collect()
}

/// Compute structural hash of a fragment (centered at atom, within radius).
/// Uses atomic properties + bond connectivity as the signature.
fn fragment_structural_hash(mol: &Molecule, atom_set: &[AtomIdx]) -> Vec<u8> {
    let mut sig = Vec::new();

    // Sort atoms by index for consistent ordering
    let mut sorted = atom_set.to_vec();
    sorted.sort_unstable();

    // Add atomic properties
    for &idx in &sorted {
        let atom = mol.atom(idx);
        sig.push(atom.element.atomic_number());
        sig.push(atom.charge.wrapping_add(8) as u8); // shift to [0, 15]
        sig.push(atom.aromatic as u8);
        let degree = mol.neighbors(idx).count() as u8;
        sig.push(degree);
    }

    // Add bond connectivity between atoms in the fragment
    let atom_set_sorted: HashSet<_> = sorted.iter().copied().collect();
    for (_, bond) in mol.bonds() {
        if atom_set_sorted.contains(&bond.atom1) && atom_set_sorted.contains(&bond.atom2) {
            sig.push(bond.atom1.0 as u8 % 255);
            sig.push(bond.atom2.0 as u8 % 255);
            let bond_type = match bond.order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                BondOrder::Quadruple => 5,
                _ => 0,
            };
            sig.push(bond_type);
        }
    }

    sig
}

/// Extract circular fragment signatures from a molecule at all radii.
fn extract_fragment_signatures(mol: &Molecule, radius: u32) -> Vec<Vec<u8>> {
    let mut signatures = Vec::new();

    for i in 0..mol.atom_count() {
        let center = AtomIdx(i as u32);
        for r in 0..=radius {
            let atom_set = atoms_within_radius(mol, center, r);
            let sig = fragment_structural_hash(mol, &atom_set);
            signatures.push(sig);
        }
    }

    signatures
}

/// Generate MinHash fingerprint from a molecule (true MHFP, v0.1.91+).
pub fn mhfp(mol: &Molecule) -> MhfpFingerprint {
    mhfp_with_config(mol, &MhfpConfig::default())
}

/// Generate MinHash fingerprint with custom configuration.
pub fn mhfp_with_config(mol: &Molecule, config: &MhfpConfig) -> MhfpFingerprint {
    let signatures = extract_fragment_signatures(mol, config.radius);

    // If no fragments, return all max values
    if signatures.is_empty() {
        return MhfpFingerprint {
            hashes: vec![u64::MAX; config.num_hashes],
            num_hashes: config.num_hashes,
        };
    }

    // Compute MinHash: for each hash function, find the minimum hash value
    let mut hashes = vec![u64::MAX; config.num_hashes];

    for (h, hash_slot) in hashes.iter_mut().enumerate() {
        let seed = config.seed.wrapping_add(h as u64);

        for sig in &signatures {
            let mut hash_data = Vec::new();
            hash_data.extend_from_slice(&seed.to_le_bytes());
            hash_data.extend_from_slice(sig);

            let hash_val = fnv1a_hash(&hash_data);

            if hash_val < *hash_slot {
                *hash_slot = hash_val;
            }
        }
    }

    MhfpFingerprint {
        hashes,
        num_hashes: config.num_hashes,
    }
}

/// Convenience function for standard MHFP generation.
pub fn mhfp_128(mol: &Molecule) -> MhfpFingerprint {
    mhfp(mol)
}

/// Calculate Tanimoto-like similarity between two molecules using MHFP.
pub fn tanimoto_mhfp(mol1: &Molecule, mol2: &Molecule) -> f64 {
    let fp1 = mhfp(mol1);
    let fp2 = mhfp(mol2);
    fp1.tanimoto(&fp2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_mhfp_consistency() {
        let mol = parse("CC").unwrap();
        let fp1 = mhfp(&mol);
        let fp2 = mhfp(&mol);

        assert_eq!(fp1.hashes, fp2.hashes);
        assert!((fp1.tanimoto(&fp2) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_mhfp_different_molecules() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("CCC").unwrap();

        let fp1 = mhfp(&mol1);
        let fp2 = mhfp(&mol2);

        let similarity = fp1.tanimoto(&fp2);
        assert!((0.0..=1.0).contains(&similarity));
    }

    #[test]
    fn test_mhfp_symmetry() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("CCC").unwrap();

        let sim12 = tanimoto_mhfp(&mol1, &mol2);
        let sim21 = tanimoto_mhfp(&mol2, &mol1);

        assert!((sim12 - sim21).abs() < 1e-10);
    }

    #[test]
    fn test_mhfp_config() {
        let mol = parse("CC").unwrap();
        let config = MhfpConfig {
            num_hashes: 64,
            seed: 42,
            radius: 2,
        };

        let fp = mhfp_with_config(&mol, &config);
        assert_eq!(fp.num_hashes, 64);
    }

    #[test]
    fn test_atoms_within_radius() {
        let mol = parse("CCCC").unwrap();
        let center = AtomIdx(1);

        let r0 = atoms_within_radius(&mol, center, 0);
        assert_eq!(r0.len(), 1); // just the center

        let r1 = atoms_within_radius(&mol, center, 1);
        assert!(r1.len() >= 2); // center + neighbors
    }

    #[test]
    fn test_extract_fragment_signatures() {
        let mol = parse("CC").unwrap();
        let sigs = extract_fragment_signatures(&mol, 2);

        // Should have signatures for each atom at each radius
        assert!(!sigs.is_empty());
    }
}
