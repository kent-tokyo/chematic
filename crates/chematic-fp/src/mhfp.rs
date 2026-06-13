//! MHFP (MinHash Fingerprint) — Canonical Circular Fragment MinHash
//!
//! Implements the MinHash algorithm (Lowe & Sayle 2013) using Morgan-style
//! circular fragment hashes as shingles:
//! 1. For each (center_atom, radius=0..r), compute a canonical circular hash
//!    via iterative Morgan expansion (atomic invariants only, no atom indices).
//! 2. MinHash over the set of u64 shingle hashes.
//! 3. Tanimoto similarity ≈ Jaccard(A, B) via MinHash theory.
//!
//! Key improvement over v0.1.91: fragment hashes no longer include raw atom
//! indices, so the fingerprint is canonical across different SMILES orderings.

use std::collections::{HashMap, HashSet};
use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};

use crate::ecfp::fnv1a as fnv1a_hash;

/// Configuration for MinHash fingerprint generation.
#[derive(Clone, Debug)]
pub struct MhfpConfig {
    /// Number of hash functions (lanes) for MinHash
    pub num_hashes: usize,
    /// Seed offset for hash functions
    pub seed: u64,
    /// Morgan radius for circular fragment extraction (default 2 = ECFP4 equivalent)
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

        let agreements = self.hashes.iter().zip(other.hashes.iter())
            .filter(|(h1, h2)| h1 == h2)
            .count();

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

/// Compute a canonical circular hash for the subgraph at (center, radius).
///
/// Uses Morgan-style iterative expansion over atomic invariants (element, charge,
/// degree, aromaticity) — no atom indices, so the result is canonical regardless
/// of the order atoms were added to the molecule.
fn circular_fragment_hash(mol: &Molecule, center: AtomIdx, radius: u32) -> u64 {
    let atom_set: Vec<AtomIdx> = atoms_within_radius(mol, center, radius);
    let atom_set_hs: HashSet<AtomIdx> = atom_set.iter().copied().collect();

    // Initial invariant: element, charge, in-fragment degree, aromaticity, implicit H count
    // H count distinguishes NH2/NH/N and OH/O, reducing false-positive collisions.
    let mut ids: HashMap<AtomIdx, u64> = atom_set.iter().map(|&idx| {
        let atom = mol.atom(idx);
        let frag_degree = mol.neighbors(idx)
            .filter(|(nb, _)| atom_set_hs.contains(nb))
            .count() as u8;
        let h_count = implicit_hcount(mol, idx).min(255) as u8;
        let h = fnv1a_hash(&[
            atom.element.atomic_number(),
            (atom.charge.wrapping_add(8)) as u8,
            frag_degree,
            atom.aromatic as u8,
            h_count,
        ]);
        (idx, h)
    }).collect();

    // Morgan expansion (radius iterations)
    for _ in 0..radius {
        let prev = ids.clone();
        for &idx in &atom_set {
            let mut nb_info: Vec<(u8, u64)> = mol.neighbors(idx)
                .filter(|(nb, _)| atom_set_hs.contains(nb))
                .map(|(nb_idx, bond_idx)| {
                    let bond_t = match mol.bond(bond_idx).order {
                        BondOrder::Single => 1u8,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        BondOrder::Aromatic => 4,
                        _ => 0,
                    };
                    (bond_t, prev[&nb_idx])
                })
                .collect();
            nb_info.sort_unstable(); // canonical order

            let mut buf = prev[&idx].to_le_bytes().to_vec();
            for (bt, nh) in &nb_info {
                buf.push(*bt);
                buf.extend_from_slice(&nh.to_le_bytes());
            }
            ids.insert(idx, fnv1a_hash(&buf));
        }
    }

    ids[&center]
}

/// Compute all circular fragment hashes for a molecule (one per center × radius).
fn extract_fragment_hashes(mol: &Molecule, radius: u32) -> Vec<u64> {
    (0..mol.atom_count())
        .flat_map(|i| {
            let center = AtomIdx(i as u32);
            (0..=radius).map(move |r| circular_fragment_hash(mol, center, r))
        })
        .collect()
}

/// Generate MinHash fingerprint from a molecule.
pub fn mhfp(mol: &Molecule) -> MhfpFingerprint {
    mhfp_with_config(mol, &MhfpConfig::default())
}

/// Generate MinHash fingerprint with custom configuration.
pub fn mhfp_with_config(mol: &Molecule, config: &MhfpConfig) -> MhfpFingerprint {
    let shingles = extract_fragment_hashes(mol, config.radius);

    if shingles.is_empty() {
        return MhfpFingerprint {
            hashes: vec![u64::MAX; config.num_hashes],
            num_hashes: config.num_hashes,
        };
    }

    // For each hash lane h, find min over all shingles of fnv1a(seed_h || shingle)
    let mut minhashes = vec![u64::MAX; config.num_hashes];
    for (h, slot) in minhashes.iter_mut().enumerate() {
        let seed = config.seed.wrapping_add(h as u64);
        let seed_bytes = seed.to_le_bytes();
        for &shingle in &shingles {
            let mut buf = seed_bytes.to_vec();
            buf.extend_from_slice(&shingle.to_le_bytes());
            let v = fnv1a_hash(&buf);
            if v < *slot {
                *slot = v;
            }
        }
    }

    MhfpFingerprint { hashes: minhashes, num_hashes: config.num_hashes }
}

/// Convenience function for standard MHFP generation.
pub fn mhfp_128(mol: &Molecule) -> MhfpFingerprint {
    mhfp(mol)
}

/// Calculate Tanimoto-like similarity between two molecules using MHFP.
pub fn tanimoto_mhfp(mol1: &Molecule, mol2: &Molecule) -> f64 {
    mhfp(mol1).tanimoto(&mhfp(mol2))
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
        let config = MhfpConfig { num_hashes: 64, seed: 42, radius: 2 };
        let fp = mhfp_with_config(&mol, &config);
        assert_eq!(fp.num_hashes, 64);
    }

    #[test]
    fn test_atoms_within_radius() {
        let mol = parse("CCCC").unwrap();
        let center = AtomIdx(1);
        let r0 = atoms_within_radius(&mol, center, 0);
        assert_eq!(r0.len(), 1);
        let r1 = atoms_within_radius(&mol, center, 1);
        assert!(r1.len() >= 2);
    }

    #[test]
    fn test_extract_fragment_hashes() {
        let mol = parse("CC").unwrap();
        let hashes = extract_fragment_hashes(&mol, 2);
        assert!(!hashes.is_empty());
    }

    /// Canonical test: the same molecule parsed from different SMILES orderings
    /// should produce the same fingerprint.
    #[test]
    fn test_mhfp_canonical_same_mol_different_parse() {
        // Ethanol: same molecule, two valid SMILES
        let mol1 = parse("CCO").unwrap();
        let mol2 = parse("OCC").unwrap();
        let fp1 = mhfp(&mol1);
        let fp2 = mhfp(&mol2);
        // With canonical Morgan hashes, identical topology → identical fingerprint
        assert!((fp1.tanimoto(&fp2) - 1.0).abs() < 1e-6,
            "Same molecule with reversed SMILES should have Tanimoto=1.0, got {}", fp1.tanimoto(&fp2));
    }

    /// Similar molecules should score higher than dissimilar ones.
    #[test]
    fn test_mhfp_similar_mols_higher_than_dissimilar() {
        let ethane = parse("CC").unwrap();
        let propane = parse("CCC").unwrap();
        let benzene = parse("c1ccccc1").unwrap();

        let sim_ec = tanimoto_mhfp(&ethane, &propane);    // similar (aliphatic C chain)
        let sim_eb = tanimoto_mhfp(&ethane, &benzene);   // dissimilar (aromatic vs aliphatic)

        assert!(sim_ec > sim_eb,
            "ethane~propane ({:.3}) should be > ethane~benzene ({:.3})", sim_ec, sim_eb);
    }

    /// H count in invariants: pyridine (no NH) vs pyrrole (has NH) should differ.
    #[test]
    fn test_mhfp_pyridine_vs_pyrrole() {
        let pyridine = parse("c1ccncc1").unwrap();
        let pyrrole  = parse("c1cc[nH]c1").unwrap();
        let sim = tanimoto_mhfp(&pyridine, &pyrrole);
        // They share aromatic ring structure, but pyridine N has no H, pyrrole does
        // Should be similar but NOT identical
        assert!(sim < 1.0,
            "pyridine and pyrrole should NOT be identical (got sim={sim:.3})");
    }

    /// H count distinguishes methylamine (NH2) from dimethylamine (NH)
    #[test]
    fn test_mhfp_amine_h_count() {
        let methylamine    = parse("CN").unwrap();   // NH2
        let dimethylamine  = parse("CNC").unwrap();  // NH
        let trimethylamine = parse("CN(C)C").unwrap(); // N (no H)
        let fp_m  = mhfp(&methylamine);
        let fp_dm = mhfp(&dimethylamine);
        let fp_tm = mhfp(&trimethylamine);
        // All three should be distinct fingerprints
        assert!(fp_m.tanimoto(&fp_dm) < 1.0);
        assert!(fp_dm.tanimoto(&fp_tm) < 1.0);
    }

    /// Larger radius should produce finer-grained differentiation for larger molecules.
    #[test]
    fn test_mhfp_radius_effect() {
        let mol = parse("c1ccccc1CC").unwrap(); // toluene-like
        let fp_r1 = mhfp_with_config(&mol, &MhfpConfig { radius: 1, ..Default::default() });
        let fp_r3 = mhfp_with_config(&mol, &MhfpConfig { radius: 3, ..Default::default() });
        // Different radii → different fingerprints
        assert_ne!(fp_r1.hashes, fp_r3.hashes);
    }
}
