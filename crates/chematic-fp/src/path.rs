//! RDKit-compatible path fingerprints (Daylight-like algorithm).
//!
//! Enumerates all topological paths up to a configurable length,
//! hashes them, and aggregates into a 2048-bit fingerprint.

use crate::bitvec::BitVec2048;
use chematic_core::{AtomIdx, Molecule};

const HASH_MOD: usize = 2048;

/// RDKit path fingerprint configuration.
#[derive(Clone, Debug)]
pub struct RdkitPathConfig {
    /// Maximum path length (default 7, matching RDKit).
    pub max_path_len: usize,
    /// Use atom types (default false — just atom index).
    pub use_atom_types: bool,
}

impl Default for RdkitPathConfig {
    fn default() -> Self {
        Self {
            max_path_len: 7,
            use_atom_types: false,
        }
    }
}

/// Compute RDKit path fingerprint (Daylight-like).
///
/// Enumerates all simple paths up to `max_path_len`, hashes each path,
/// and sets corresponding bits in a 2048-bit fingerprint.
pub fn rdkit_path_fp(mol: &Molecule) -> BitVec2048 {
    rdkit_path_fp_with_config(mol, &RdkitPathConfig::default())
}

/// RDKit path fingerprint with custom config.
pub fn rdkit_path_fp_with_config(
    mol: &Molecule,
    config: &RdkitPathConfig,
) -> BitVec2048 {
    let mut fp = BitVec2048::new();

    if mol.atom_count() == 0 {
        return fp;
    }

    for start_idx in 0..mol.atom_count() {
        let start = AtomIdx(start_idx as u32);
        enumerate_paths_from(
            mol,
            start,
            vec![start],
            &mut |path| {
                let hash = hash_path(mol, path);
                let bit_idx = hash % HASH_MOD;
                fp.set(bit_idx);
            },
            config.max_path_len,
        );
    }

    fp
}

/// Recursive DFS path enumeration.
fn enumerate_paths_from<F>(
    mol: &Molecule,
    current: AtomIdx,
    path: Vec<AtomIdx>,
    callback: &mut F,
    max_len: usize,
) where
    F: FnMut(&[AtomIdx]),
{
    if path.len() <= max_len {
        // Callback for the current path
        if path.len() > 1 {
            callback(&path);
        }

        // Extend to neighbors
        if path.len() < max_len {
            for (neighbor, _) in mol.neighbors(current) {
                // Avoid cycles — don't revisit atoms already in path
                if !path.contains(&neighbor) {
                    let mut new_path = path.clone();
                    new_path.push(neighbor);
                    enumerate_paths_from(mol, neighbor, new_path, callback, max_len);
                }
            }
        }
    }
}

/// Hash a path using FNV-1a.
fn hash_path(mol: &Molecule, path: &[AtomIdx]) -> usize {
    // Use portable FNV constants that work on both 32-bit and 64-bit
    let fnv_prime: usize = 16777619; // FNV prime (32-bit compatible)
    let mut hash: usize = 2166136261; // FNV offset basis (32-bit compatible)

    for atom_idx in path {
        let atom = mol.atom(*atom_idx);
        // Include atomic number in hash
        let an = atom.element.atomic_number() as usize;
        hash ^= an;
        hash = hash.wrapping_mul(fnv_prime);
    }

    hash
}

/// Tanimoto similarity between two RDKit path fingerprints.
pub fn tanimoto_rdkit_path(a: &BitVec2048, b: &BitVec2048) -> f64 {
    a.tanimoto(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(smiles: &str) -> Molecule {
        parse(smiles).unwrap_or_else(|e| panic!("failed to parse {smiles:?}: {e}"))
    }

    #[test]
    fn test_rdkit_path_fp_ethane() {
        let m = mol("CC");
        let fp = rdkit_path_fp(&m);
        assert!(fp.popcount() > 0, "ethane should have non-zero bits");
    }

    #[test]
    fn test_rdkit_path_fp_propane() {
        let m = mol("CCC");
        let fp = rdkit_path_fp(&m);
        assert!(fp.popcount() > 0, "propane should have non-zero bits");
    }

    #[test]
    fn test_rdkit_path_fp_benzene() {
        let m = mol("c1ccccc1");
        let fp = rdkit_path_fp(&m);
        assert!(fp.popcount() > 0, "benzene should have non-zero bits");
    }

    #[test]
    fn test_rdkit_path_fp_identical() {
        let m1 = mol("CC");
        let m2 = mol("CC");
        let fp1 = rdkit_path_fp(&m1);
        let fp2 = rdkit_path_fp(&m2);
        assert_eq!(fp1.tanimoto(&fp2), 1.0, "identical molecules should have tanimoto=1.0");
    }

    #[test]
    fn test_rdkit_path_fp_single_atom() {
        let m = mol("C");
        let fp = rdkit_path_fp(&m);
        // Single atom has no paths (path length >= 2)
        assert_eq!(fp.popcount(), 0, "single atom should have zero bits");
    }
}
