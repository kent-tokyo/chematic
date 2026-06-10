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
    let neighbors: Vec<_> = mol.neighbors(idx).collect();
    let degree = neighbors.len();

    hash ^= degree;
    hash = hash.wrapping_mul(fnv_prime);

    // Hash neighbor atomic numbers and bond types
    for (neighbor_idx, bond_idx) in neighbors {
        let neighbor = mol.atom(neighbor_idx);
        let neighbor_an = neighbor.element.atomic_number() as usize;
        let bond = mol.bond(bond_idx);
        let bond_order = match bond.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            _ => 1,
        };

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
        assert_eq!(fp1.tanimoto(&fp2), 1.0, "identical molecules should have tanimoto=1.0");
    }

    #[test]
    fn test_pattern_fp_single_atom() {
        let m = mol("C");
        let fp = pattern_fp(&m);
        assert!(fp.popcount() > 0, "single atom should have bits");
    }
}
