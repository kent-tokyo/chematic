//! Topological path fingerprints.
//!
//! Enumerates all simple paths of length 1..=`max_len` bonds and hashes
//! them with FNV-1a. Paths are canonicalised (forward vs. reverse) so that
//! each unique path is hashed only once.

use chematic_core::{AtomIdx, BondIdx, Molecule};

use crate::bitvec::BitVec2048;
use crate::ecfp::{bond_type_int as bond_byte, fnv1a};

/// Configuration for topological path fingerprints.
#[derive(Debug, Clone)]
pub struct TopoPathConfig {
    /// Maximum number of atoms in a path (default: 7).
    pub max_len: usize,
    /// Bitvector size in bits (default: 2048).
    pub nbits: usize,
}

impl Default for TopoPathConfig {
    fn default() -> Self {
        Self {
            max_len: 7,
            nbits: 2048,
        }
    }
}

/// Compute a topological path fingerprint for `mol`.
///
/// All simple paths of 2..=`config.max_len` atoms are enumerated,
/// encoded as interleaved `[atomic_num, bond_order, atomic_num, ...]` bytes,
/// canonicalised (min of forward/reverse encoding), hashed with FNV-1a, and
/// folded into a `BitVec2048` via `hash % nbits`.
pub fn topo_path(mol: &Molecule, config: &TopoPathConfig) -> BitVec2048 {
    let mut fp = BitVec2048::new();
    if mol.atom_count() == 0 {
        return fp;
    }
    let n = mol.atom_count();
    let mut path_atoms: Vec<u8> = Vec::with_capacity(config.max_len);
    let mut path_bonds: Vec<u8> = Vec::with_capacity(config.max_len.saturating_sub(1));
    let mut visited: Vec<bool> = vec![false; n];

    for start in 0..n {
        let start_idx = AtomIdx(start as u32);
        path_atoms.push(mol.atom(start_idx).element.atomic_number());
        visited[start] = true;
        dfs(
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
fn dfs(
    mol: &Molecule,
    atom: AtomIdx,
    from_bond: Option<BondIdx>,
    path_atoms: &mut Vec<u8>,
    path_bonds: &mut Vec<u8>,
    visited: &mut Vec<bool>,
    fp: &mut BitVec2048,
    config: &TopoPathConfig,
) {
    // Record path when it has at least 2 atoms (1 bond)
    if path_atoms.len() >= 2 {
        hash_path(path_atoms, path_bonds, fp, config.nbits);
    }
    if path_atoms.len() >= config.max_len {
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
        path_bonds.push(bond_byte(order));
        path_atoms.push(mol.atom(nbr).element.atomic_number());
        visited[nbr.0 as usize] = true;
        dfs(
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
    // Interleaved encoding: [a0, b0, a1, b1, ..., a_n-1]
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
    let hash = fnv1a(&canonical);
    fp.set((hash % nbits as u64) as usize);
}

/// Tanimoto similarity between two molecules using topological path fingerprints.
pub fn tanimoto_topo_path(a: &Molecule, b: &Molecule) -> f64 {
    let cfg = TopoPathConfig::default();
    topo_path(a, &cfg).tanimoto(&topo_path(b, &cfg))
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn cfg() -> TopoPathConfig {
        TopoPathConfig::default()
    }

    #[test]
    fn ethane_nonzero() {
        let mol = parse("CC").unwrap();
        assert!(topo_path(&mol, &cfg()).popcount() > 0);
    }

    #[test]
    fn benzene_nonzero() {
        let mol = parse("c1ccccc1").unwrap();
        assert!(topo_path(&mol, &cfg()).popcount() > 0);
    }

    #[test]
    fn aspirin_many_bits() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let fp = topo_path(&mol, &cfg());
        assert!(
            fp.popcount() > 10,
            "aspirin topo_path got {}",
            fp.popcount()
        );
    }

    #[test]
    fn deterministic() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(topo_path(&mol, &cfg()), topo_path(&mol, &cfg()));
    }

    #[test]
    fn longer_max_len_more_bits() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let fp2 = topo_path(
            &mol,
            &TopoPathConfig {
                max_len: 2,
                nbits: 2048,
            },
        );
        let fp7 = topo_path(&mol, &cfg());
        assert!(
            fp2.popcount() <= fp7.popcount(),
            "max_len=2 ({}) should have <= bits than max_len=7 ({})",
            fp2.popcount(),
            fp7.popcount()
        );
    }

    #[test]
    fn ethane_vs_benzene_differ() {
        let e = parse("CC").unwrap();
        let b = parse("c1ccccc1").unwrap();
        let t = topo_path(&e, &cfg()).tanimoto(&topo_path(&b, &cfg()));
        assert!(t < 1.0, "ethane and benzene must differ (tanimoto={t})");
    }
}
