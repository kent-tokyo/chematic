//! Avalon-style structural fingerprint.
//!
//! Loosely modelled on the Avalon Cheminformatics Toolkit fingerprint
//! (`rdkit.Avalon.pyAvalonTools.GetAvalonFP`): a broad mix of atom, bond,
//! ring, and path features hashed into a fixed-width bitset. This is not a
//! bit-exact reimplementation of RDKit's Avalon fingerprint — like Morgan
//! (FNV-1a here vs RDKit's MurmurHash), bit positions differ; the bar is
//! deterministic, isomorphism-invariant, similarity-preserving hashing, not
//! cross-library bit parity.

use chematic_core::Molecule;
use chematic_perception::find_sssr;

use crate::bitvec::BitVec2048;
use crate::ecfp::{EcfpInvariantMode, bond_type_int, fnv1a, initial_atom_id};

/// Configuration for Avalon fingerprint computation.
#[derive(Debug, Clone)]
pub struct AvalonConfig {
    /// Maximum number of atoms in a path feature (default: 7).
    pub max_path_len: usize,
    /// Bitvector size in bits (default: 2048).
    pub nbits: usize,
}

impl Default for AvalonConfig {
    fn default() -> Self {
        Self {
            max_path_len: 7,
            nbits: 2048,
        }
    }
}

// Category tags keep identical byte sequences from different feature kinds
// from colliding into the same hash input.
const TAG_ATOM: u8 = 0x01;
const TAG_BOND: u8 = 0x02;
const TAG_RING: u8 = 0x03;
const TAG_PATH: u8 = 0x04;

fn set_feature_bytes(tag: u8, bytes: &[u8], fp: &mut BitVec2048, nbits: usize) {
    let mut tagged = Vec::with_capacity(bytes.len() + 1);
    tagged.push(tag);
    tagged.extend_from_slice(bytes);
    let hash = fnv1a(&tagged);
    fp.set((hash % nbits as u64) as usize);
}

/// Compute an Avalon-style fingerprint for `mol` using the default configuration.
pub fn avalon_fp(mol: &Molecule) -> BitVec2048 {
    avalon_fp_with_config(mol, &AvalonConfig::default())
}

/// Compute an Avalon-style fingerprint for `mol` with explicit configuration.
pub fn avalon_fp_with_config(mol: &Molecule, config: &AvalonConfig) -> BitVec2048 {
    let mut fp = BitVec2048::new();
    if mol.atom_count() == 0 {
        return fp;
    }
    let ring_set = find_sssr(mol);
    let nbits = config.nbits;

    // 1. Atom features.
    for (idx, _) in mol.atoms() {
        let id = initial_atom_id(mol, idx, &ring_set, false, EcfpInvariantMode::Chematic);
        set_feature_bytes(TAG_ATOM, &id.to_le_bytes(), &mut fp, nbits);
    }

    // 2. Bond features: (bond type, aromatic, in-ring, sorted endpoint degrees).
    for (_, bond) in mol.bonds() {
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        let deg1 = mol.neighbors(a1).count().min(255) as u8;
        let deg2 = mol.neighbors(a2).count().min(255) as u8;
        let (deg_lo, deg_hi) = if deg1 <= deg2 {
            (deg1, deg2)
        } else {
            (deg2, deg1)
        };
        let in_ring = (ring_set.contains_atom(a1) && ring_set.contains_atom(a2)) as u8;
        let bytes = [
            bond_type_int(bond.order),
            (bond.order == chematic_core::BondOrder::Aromatic) as u8,
            in_ring,
            deg_lo,
            deg_hi,
        ];
        set_feature_bytes(TAG_BOND, &bytes, &mut fp, nbits);
    }

    // 3. Ring features: (size, all-aromatic flag, sorted atomic-number multiset).
    for ring in ring_set.rings() {
        let size = ring.len().min(255) as u8;
        let n = ring.len();
        let all_aromatic = (0..n).all(|i| {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            mol.bond_between(a, b)
                .is_some_and(|(_, bd)| bd.order == chematic_core::BondOrder::Aromatic)
        });
        let mut elems: Vec<u8> = ring
            .iter()
            .map(|&a| mol.atom(a).element.atomic_number())
            .collect();
        elems.sort_unstable();
        let mut bytes = Vec::with_capacity(2 + elems.len());
        bytes.push(size);
        bytes.push(all_aromatic as u8);
        bytes.extend_from_slice(&elems);
        set_feature_bytes(TAG_RING, &bytes, &mut fp, nbits);
    }

    // 4. Path features: simple paths of 2..=max_path_len atoms, canonicalised
    // forward-vs-reverse (same scheme as topo_path::hash_path).
    let n = mol.atom_count();
    let mut path_atoms: Vec<u8> = Vec::with_capacity(config.max_path_len);
    let mut path_bonds: Vec<u8> = Vec::with_capacity(config.max_path_len.saturating_sub(1));
    let mut visited: Vec<bool> = vec![false; n];
    for start in 0..n {
        let start_idx = chematic_core::AtomIdx(start as u32);
        path_atoms.push(mol.atom(start_idx).element.atomic_number());
        visited[start] = true;
        path_dfs(
            mol,
            start_idx,
            None,
            &mut path_atoms,
            &mut path_bonds,
            &mut visited,
            &mut fp,
            config,
        );
        path_atoms.pop();
        visited[start] = false;
    }

    fp
}

#[allow(clippy::too_many_arguments)]
fn path_dfs(
    mol: &Molecule,
    atom: chematic_core::AtomIdx,
    from_bond: Option<chematic_core::BondIdx>,
    path_atoms: &mut Vec<u8>,
    path_bonds: &mut Vec<u8>,
    visited: &mut Vec<bool>,
    fp: &mut BitVec2048,
    config: &AvalonConfig,
) {
    if path_atoms.len() >= 2 {
        hash_path(path_atoms, path_bonds, fp, config.nbits);
    }
    if path_atoms.len() >= config.max_path_len {
        return;
    }
    for (nbr, bid) in mol.neighbors(atom) {
        if Some(bid) == from_bond {
            continue;
        }
        if visited[nbr.0 as usize] {
            continue;
        }
        let order = mol.bond(bid).order;
        path_bonds.push(bond_type_int(order));
        path_atoms.push(mol.atom(nbr).element.atomic_number());
        visited[nbr.0 as usize] = true;
        path_dfs(
            mol,
            nbr,
            Some(bid),
            path_atoms,
            path_bonds,
            visited,
            fp,
            config,
        );
        visited[nbr.0 as usize] = false;
        path_atoms.pop();
        path_bonds.pop();
    }
}

fn hash_path(atoms: &[u8], bonds: &[u8], fp: &mut BitVec2048, nbits: usize) {
    let len = atoms.len() * 2 - 1;
    let mut fwd = Vec::with_capacity(len);
    let mut rev = Vec::with_capacity(len);
    for i in 0..atoms.len() {
        fwd.push(atoms[i]);
        if i + 1 < atoms.len() {
            fwd.push(bonds[i]);
        }
    }
    for i in (0..atoms.len()).rev() {
        rev.push(atoms[i]);
        if i > 0 {
            rev.push(bonds[i - 1]);
        }
    }
    let canonical = if fwd <= rev { fwd } else { rev };
    set_feature_bytes(TAG_PATH, &canonical, fp, nbits);
}

/// Tanimoto similarity between two molecules using Avalon-style fingerprints.
pub fn tanimoto_avalon(a: &Molecule, b: &Molecule) -> f64 {
    avalon_fp(a).tanimoto(&avalon_fp(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn deterministic() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        assert_eq!(avalon_fp(&mol), avalon_fp(&mol));
    }

    #[test]
    fn isomorphism_invariant() {
        // Same molecule (aspirin), two different SMILES atom orderings.
        let a = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let b = parse("O=C(C)Oc1ccccc1C(=O)O").unwrap();
        assert_eq!(avalon_fp(&a), avalon_fp(&b));
    }

    #[test]
    fn non_degenerate() {
        for smi in ["c1ccccc1", "Cc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"] {
            let mol = parse(smi).unwrap();
            let fp = avalon_fp(&mol);
            assert!(fp.popcount() > 0, "{smi}: fingerprint must be non-empty");
            assert!(fp.popcount() < 2048, "{smi}: fingerprint must not saturate");
        }
    }

    #[test]
    fn similarity_ordering() {
        let benzene = parse("c1ccccc1").unwrap();
        let toluene = parse("Cc1ccccc1").unwrap();
        let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let t_toluene = avalon_fp(&benzene).tanimoto(&avalon_fp(&toluene));
        let t_aspirin = avalon_fp(&benzene).tanimoto(&avalon_fp(&aspirin));
        assert!(
            t_toluene > t_aspirin,
            "toluene ({t_toluene}) should be closer to benzene than aspirin ({t_aspirin})"
        );
    }

    #[test]
    fn bond_order_sensitivity() {
        let single = parse("CC").unwrap();
        let double = parse("C=C").unwrap();
        assert_ne!(avalon_fp(&single), avalon_fp(&double));
    }

    #[test]
    fn aromaticity_sensitivity() {
        let aromatic = parse("c1ccccc1").unwrap();
        let aliphatic = parse("C1CCCCC1").unwrap();
        assert_ne!(avalon_fp(&aromatic), avalon_fp(&aliphatic));
    }
}
