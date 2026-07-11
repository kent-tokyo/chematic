//! Pattern fingerprints — substructure feature hashing.
//!
//! Enumerates atom-centric patterns (neighborhood topology) and hashes them
//! into a 2048-bit fingerprint for fast substructure-based screening.

use crate::bitvec::BitVec2048;
use chematic_core::{AtomIdx, BondOrder, Molecule};

const HASH_MOD: usize = 2048;

/// Compute pattern fingerprint (2048-bit).
///
/// Hashes atom-centric patterns: atomic number + neighbor count + bond types.
/// Suitable for fast substructure screening and molecular similarity.
pub fn pattern_fp(mol: &Molecule) -> BitVec2048 {
    let mut fp = BitVec2048::new();

    if mol.atom_count() == 0 {
        return fp;
    }

    for (idx, _atom) in mol.atoms() {
        // Compute pattern hash for this atom's local neighborhood
        let pattern_hash = compute_pattern_hash(mol, idx);
        let bit_idx = pattern_hash % HASH_MOD;
        fp.set(bit_idx);
    }

    fp
}

/// Compute hash for atom's pattern: atomic number + degree + neighbor types.
fn compute_pattern_hash(mol: &Molecule, idx: AtomIdx) -> usize {
    let fnv_prime: usize = 16777619;
    let mut hash: usize = 2166136261; // FNV offset basis

    let atom = mol.atom(idx);
    let an = atom.element.atomic_number() as usize;

    // Hash atomic number
    hash ^= an;
    hash = hash.wrapping_mul(fnv_prime);

    // Hash neighbor information
    let mut neighbor_info: Vec<(usize, usize)> = mol
        .neighbors(idx)
        .map(|(neighbor_idx, bond_idx)| {
            let neighbor_an = mol.atom(neighbor_idx).element.atomic_number() as usize;
            let bond_order = match mol.bond(bond_idx).order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                _ => 1,
            };
            (neighbor_an, bond_order)
        })
        .collect();
    // Sort so the hash depends only on the multiset of neighbor
    // (atomic_number, bond_order) pairs, not on `mol.neighbors`' iteration
    // order (which reflects parse/insertion order, not anything canonical).
    neighbor_info.sort_unstable();
    let degree = neighbor_info.len();

    hash ^= degree;
    hash = hash.wrapping_mul(fnv_prime);

    for (neighbor_an, bond_order) in neighbor_info {
        hash ^= neighbor_an;
        hash = hash.wrapping_mul(fnv_prime);
        hash ^= bond_order;
        hash = hash.wrapping_mul(fnv_prime);
    }

    // Hash aromaticity
    if atom.aromatic {
        hash ^= 1;
        hash = hash.wrapping_mul(fnv_prime);
    }

    hash
}

/// Tanimoto similarity between two pattern fingerprints.
pub fn tanimoto_pattern(a: &BitVec2048, b: &BitVec2048) -> f64 {
    a.tanimoto(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(smiles: &str) -> Molecule {
        parse(smiles).unwrap_or_else(|e| panic!("parse '{smiles}': {e}"))
    }

    #[test]
    fn test_pattern_fp_ethane() {
        let m = mol("CC");
        let fp = pattern_fp(&m);
        assert!(fp.popcount() > 0, "ethane should have non-zero bits");
    }

    #[test]
    fn test_pattern_fp_order_independent() {
        // Acetamide's carbonyl carbon has 3 chemically distinct neighbors
        // (=O, -N, -C) -- their traversal order differs between these two
        // spellings (which neighbor gets bonded to the carbonyl carbon
        // first during parsing differs), so this is real tie material for
        // compute_pattern_hash's neighbor fold, not a harmless symmetric
        // case. Same molecule must give the same fingerprint regardless of
        // input spelling.
        let a = pattern_fp(&mol("CC(=O)N"));
        let b = pattern_fp(&mol("NC(=O)C"));
        assert_eq!(a, b, "pattern_fp must not depend on SMILES atom order");
    }

    #[test]
    fn test_pattern_fp_benzene() {
        let m = mol("c1ccccc1");
        let fp = pattern_fp(&m);
        assert!(fp.popcount() > 0, "benzene should have non-zero bits");
    }

    #[test]
    fn test_pattern_fp_similarity() {
        let m1 = mol("CC");
        let m2 = mol("CC");
        let fp1 = pattern_fp(&m1);
        let fp2 = pattern_fp(&m2);
        assert_eq!(
            fp1.tanimoto(&fp2),
            1.0,
            "identical molecules should have tanimoto=1.0"
        );
    }

    #[test]
    fn test_pattern_fp_single_atom() {
        let m = mol("C");
        let fp = pattern_fp(&m);
        assert!(fp.popcount() > 0, "single atom should have bits");
    }
}
