//! MAP4 — MinHashed Atom-Pair Fingerprint (Minervini et al. 2020)
//!
//! MAP4 encodes circular atom environments as SMILES "shingles" and hashes them
//! using the MinHash scheme. Unlike MHFP (which uses 3-hop paths), MAP4 uses
//! all atom-pair combinations at increasing topological radius.
//!
//! Reference: Minervini et al., *J. Cheminform.* 2020, 12, 26.
//! "One molecular fingerprint to rule them all"
//!
//! ## Algorithm
//! 1. For each atom pair (a, b) and radii (r_a, r_b) from 0 to max_radius:
//!    - Extract the circular SMILES environment of atom a at radius r_a
//!    - Extract the circular SMILES environment of atom b at radius r_b
//!    - Encode the pair as `"env_a|env_b|distance"` (sorted for symmetry)
//! 2. Collect all shingles into a set
//! 3. Apply MinHash with `n_permutations` hash functions → fixed-size signature
//!
//! ## Comparison to RDKit
//! RDKit does **not** include MAP4 in its main distribution — users must install
//! the separate `map4` package. chematic implements it natively in pure Rust.

use crate::ecfp::fnv1a;
use chematic_core::{AtomIdx, Molecule};

/// MAP4 fingerprint configuration.
#[derive(Debug, Clone)]
pub struct Map4Config {
    /// Maximum circular environment radius (default 2, like ECFP4)
    pub max_radius: usize,
    /// Number of MinHash permutations (default 1024)
    pub n_permutations: usize,
}

impl Default for Map4Config {
    fn default() -> Self {
        Self {
            max_radius: 2,
            n_permutations: 1024,
        }
    }
}

/// Compute the MAP4 fingerprint as a MinHash signature vector.
///
/// Returns a `Vec<u32>` of length `config.n_permutations`. Two molecules
/// with Jaccard similarity J have expected Hamming similarity J over this vector.
///
/// Tanimoto similarity can be estimated as:
/// `(matching_positions as f64) / config.n_permutations as f64`
pub fn map4(mol: &Molecule, config: &Map4Config) -> Vec<u32> {
    let shingles = collect_shingles(mol, config.max_radius);
    minhash(&shingles, config.n_permutations)
}

/// Compute MAP4 with default configuration (radius=2, 1024 permutations).
pub fn map4_default(mol: &Molecule) -> Vec<u32> {
    map4(mol, &Map4Config::default())
}

/// Estimate Tanimoto similarity between two MAP4 signatures.
///
/// Returns a value in [0, 1]. Exact for Jaccard similarity of the shingle sets
/// in the limit of infinite permutations; approximate with finite n_permutations.
pub fn tanimoto_map4(a: &[u32], b: &[u32]) -> f64 {
    if a.len() != b.len() || a.is_empty() {
        return 0.0;
    }
    let matches = a.iter().zip(b.iter()).filter(|(x, y)| x == y).count();
    matches as f64 / a.len() as f64
}

// ─── Shingle extraction ───────────────────────────────────────────────────────

fn collect_shingles(mol: &Molecule, max_radius: usize) -> Vec<u64> {
    let n = mol.atom_count();
    // Precompute BFS distances between all pairs (Floyd-Warshall on adjacency)
    let dist = bfs_all_pairs(mol);

    let mut hashes = Vec::new();

    for a in 0..n {
        for b in a..n {
            let d = dist[a][b];
            if d == usize::MAX {
                continue; // disconnected
            }
            for r_a in 0..=max_radius {
                for r_b in 0..=max_radius {
                    let env_a = circular_env_hash(mol, a, r_a, &dist);
                    let env_b = circular_env_hash(mol, b, r_b, &dist);
                    // Canonical: sort envs so order of (a,b) doesn't matter
                    let (ea, eb) = if env_a <= env_b {
                        (env_a, env_b)
                    } else {
                        (env_b, env_a)
                    };
                    // Encode: hash(env_a || env_b || distance)
                    let mut h = fnv1a(ea.to_le_bytes().as_ref());
                    h = fnv1a_extend(h, eb.to_le_bytes().as_ref());
                    h = fnv1a_extend(h, &(d as u64).to_le_bytes());
                    hashes.push(h);
                }
            }
        }
    }
    hashes.sort_unstable();
    hashes.dedup();
    hashes
}

/// Hash for the circular environment of `center` at `radius`.
/// Encodes atom types of all atoms within `radius` bonds, together with the
/// topological distances from the center.
fn circular_env_hash(mol: &Molecule, center: usize, radius: usize, dist: &[Vec<usize>]) -> u64 {
    let mut contributions: Vec<(usize, u64)> = mol
        .atoms()
        .filter_map(|(idx, atom)| {
            let d = dist[center][idx.0 as usize];
            if d <= radius {
                let atom_hash = fnv1a(&[
                    atom.element.atomic_number(),
                    atom.charge.unsigned_abs(),
                    atom.element.atomic_number().wrapping_add(d as u8),
                ]);
                Some((d, atom_hash))
            } else {
                None
            }
        })
        .collect();
    // Sort by (distance, hash) for determinism
    contributions.sort_unstable();
    let mut h: u64 = 0xcbf29ce484222325; // FNV offset basis
    for (d, ah) in contributions {
        h = fnv1a_extend(h, &(d as u64).to_le_bytes());
        h = fnv1a_extend(h, &ah.to_le_bytes());
    }
    h
}

// ─── MinHash ─────────────────────────────────────────────────────────────────

fn minhash(shingles: &[u64], n_permutations: usize) -> Vec<u32> {
    let mut sig = vec![u32::MAX; n_permutations];
    for &h in shingles {
        for (i, slot) in sig.iter_mut().enumerate().take(n_permutations) {
            // Deterministic hash family: H_i(x) = (a_i * x + b_i) mod p
            // We use FNV-based mixing for speed (no true universal hash needed)
            let mixed = fnv1a_mix(h, i as u64);
            let v = (mixed >> 32) as u32;
            if v < *slot {
                *slot = v;
            }
        }
    }
    sig
}

// ─── Graph distance ───────────────────────────────────────────────────────────

fn bfs_all_pairs(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atom_count();
    (0..n).map(|src| bfs_from(mol, src, n)).collect()
}

fn bfs_from(mol: &Molecule, src: usize, n: usize) -> Vec<usize> {
    let mut dist = vec![usize::MAX; n];
    dist[src] = 0;
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(src);
    while let Some(cur) = queue.pop_front() {
        let d = dist[cur];
        for (nb, _) in mol.neighbors(AtomIdx(cur as u32)) {
            let nbi = nb.0 as usize;
            if dist[nbi] == usize::MAX {
                dist[nbi] = d + 1;
                queue.push_back(nbi);
            }
        }
    }
    dist
}

// ─── Hash utilities ───────────────────────────────────────────────────────────

#[inline]
fn fnv1a_extend(mut h: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x00000100000001b3);
    }
    h
}

#[inline]
fn fnv1a_mix(h: u64, seed: u64) -> u64 {
    let mut v = h.wrapping_add(seed.wrapping_mul(0x9e3779b97f4a7c15));
    v ^= v >> 30;
    v = v.wrapping_mul(0xbf58476d1ce4e5b9);
    v ^= v >> 27;
    v = v.wrapping_mul(0x94d049bb133111eb);
    v ^= v >> 31;
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn map4_same_molecule_identical() {
        let mol = parse("CC").expect("parse");
        let a = map4_default(&mol);
        let b = map4_default(&mol);
        assert_eq!(a, b, "MAP4 must be deterministic");
    }

    #[test]
    fn map4_self_similarity_is_one() {
        let mol = parse("c1ccccc1").expect("parse benzene");
        let fp = map4_default(&mol);
        assert!((tanimoto_map4(&fp, &fp) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn map4_similar_molecules_higher_than_dissimilar() {
        let benzene = parse("c1ccccc1").expect("benzene");
        let toluene = parse("Cc1ccccc1").expect("toluene");
        let methane = parse("C").expect("methane");
        let fp_benz = map4_default(&benzene);
        let fp_tol = map4_default(&toluene);
        let fp_meth = map4_default(&methane);
        let sim_close = tanimoto_map4(&fp_benz, &fp_tol);
        let sim_far = tanimoto_map4(&fp_benz, &fp_meth);
        assert!(
            sim_close > sim_far,
            "benzene-toluene ({:.3}) should be > benzene-methane ({:.3})",
            sim_close,
            sim_far
        );
    }

    #[test]
    fn map4_length_equals_n_permutations() {
        let mol = parse("CCO").expect("parse");
        let cfg = Map4Config {
            max_radius: 1,
            n_permutations: 512,
        };
        let fp = map4(&mol, &cfg);
        assert_eq!(fp.len(), 512);
    }

    #[test]
    fn tanimoto_bounds() {
        let mol1 = parse("CC").expect("parse");
        let mol2 = parse("CCC").expect("parse");
        let fp1 = map4_default(&mol1);
        let fp2 = map4_default(&mol2);
        let t = tanimoto_map4(&fp1, &fp2);
        assert!((0.0..=1.0).contains(&t), "Tanimoto out of range: {}", t);
    }
}
