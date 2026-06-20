//! v0.1.92: RDKit-compatible path fingerprints with bond type inclusion.
//!
//! Enumerates all topological paths up to a configurable length,
//! hashes each path (including bond types), and aggregates into a 2048-bit fingerprint.
//! Bond types (single/double/triple/aromatic) are interleaved with atomic numbers in hash.

use crate::bitvec::BitVec2048;
use chematic_core::{AtomIdx, BondOrder, Molecule};

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

/// Compute RDKit path fingerprint (Daylight-like, with bond types).
///
/// Enumerates all simple paths up to `max_path_len`, hashes each path
/// (including bond order types), and sets corresponding bits in a 2048-bit fingerprint.
pub fn rdkit_path_fp(mol: &Molecule) -> BitVec2048 {
    rdkit_path_fp_with_config(mol, &RdkitPathConfig::default())
}

/// RDKit path fingerprint with custom config.
pub fn rdkit_path_fp_with_config(mol: &Molecule, config: &RdkitPathConfig) -> BitVec2048 {
    let mut fp = BitVec2048::new();

    if mol.atom_count() == 0 {
        return fp;
    }

    for start_idx in 0..mol.atom_count() {
        let start = AtomIdx(start_idx as u32);
        let initial_path = PathWithBonds {
            atoms: vec![start],
            bonds: vec![],
        };
        enumerate_paths_from(
            mol,
            start,
            initial_path,
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

/// Path with atoms and their connecting bond orders.
struct PathWithBonds {
    atoms: Vec<AtomIdx>,
    bonds: Vec<BondOrder>, // len = atoms.len() - 1
}

/// Recursive DFS path enumeration with bond tracking.
fn enumerate_paths_from<F>(
    mol: &Molecule,
    current: AtomIdx,
    path: PathWithBonds,
    callback: &mut F,
    max_len: usize,
) where
    F: FnMut(&PathWithBonds),
{
    if path.atoms.len() <= max_len {
        // Callback for the current path
        if path.atoms.len() > 1 {
            callback(&path);
        }

        // Extend to neighbors
        if path.atoms.len() < max_len {
            for (neighbor, bond_idx) in mol.neighbors(current) {
                // Avoid cycles — don't revisit atoms already in path
                if !path.atoms.contains(&neighbor) {
                    let bond_order = mol.bond(bond_idx).order;
                    let mut new_path = path.clone();
                    new_path.atoms.push(neighbor);
                    new_path.bonds.push(bond_order);
                    enumerate_paths_from(mol, neighbor, new_path, callback, max_len);
                }
            }
        }
    }
}

impl Clone for PathWithBonds {
    fn clone(&self) -> Self {
        PathWithBonds {
            atoms: self.atoms.clone(),
            bonds: self.bonds.clone(),
        }
    }
}

/// Hash a path using FNV-1a, interleaving atoms and bond types.
fn hash_path(mol: &Molecule, path: &PathWithBonds) -> usize {
    let mut bytes: Vec<u8> = Vec::with_capacity(path.atoms.len() * 2);
    for (i, &atom_idx) in path.atoms.iter().enumerate() {
        if i > 0 {
            bytes.push(crate::ecfp::bond_type_int(path.bonds[i - 1]));
        }
        bytes.push(mol.atom(atom_idx).element.atomic_number());
    }
    crate::ecfp::fnv1a(&bytes) as usize
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
        assert_eq!(
            fp1.tanimoto(&fp2),
            1.0,
            "identical molecules should have tanimoto=1.0"
        );
    }

    #[test]
    fn test_rdkit_path_fp_single_atom() {
        let m = mol("C");
        let fp = rdkit_path_fp(&m);
        // Single atom has no paths (path length >= 2)
        assert_eq!(fp.popcount(), 0, "single atom should have zero bits");
    }

    #[test]
    fn test_rdkit_path_fp_bond_type_distinction() {
        // Single vs double bond: same atoms, different bond types → different FP
        let m_single = mol("CC");
        let m_double = mol("C=C");
        let fp_single = rdkit_path_fp(&m_single);
        let fp_double = rdkit_path_fp(&m_double);
        // FPs should differ due to bond type
        assert!(
            tanimoto_rdkit_path(&fp_single, &fp_double) < 1.0,
            "single and double bonds should produce different fingerprints"
        );
    }

    #[test]
    fn test_rdkit_path_fp_bond_type_consistency() {
        // Identical molecules with same bond types should have identical FPs
        let m1 = mol("C=C");
        let m2 = mol("C=C");
        let fp1 = rdkit_path_fp(&m1);
        let fp2 = rdkit_path_fp(&m2);
        assert_eq!(
            tanimoto_rdkit_path(&fp1, &fp2),
            1.0,
            "identical bond types should produce identical fingerprints"
        );
    }

    #[test]
    fn test_rdkit_path_fp_triple_bond() {
        let m = mol("C#C");
        let fp = rdkit_path_fp(&m);
        assert!(fp.popcount() > 0, "triple bond should have non-zero bits");
    }

    #[test]
    fn test_rdkit_path_fp_aromatic_distinction() {
        // Benzene (aromatic) vs cyclohexatriene-like (not aromatic)
        let m_aromatic = mol("c1ccccc1");
        let m_aliphatic = mol("C1CCCC=CC=C1"); // approximation, just for structure
        let fp_aromatic = rdkit_path_fp(&m_aromatic);
        let fp_aliphatic = rdkit_path_fp(&m_aliphatic);
        // Should differ due to aromatic vs single/double bonds
        assert!(
            tanimoto_rdkit_path(&fp_aromatic, &fp_aliphatic) < 1.0,
            "aromatic and aliphatic should produce different fingerprints"
        );
    }
}
