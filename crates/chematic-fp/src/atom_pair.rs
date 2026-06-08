//! Atom-Pair and Topological Torsion fingerprints.
//!
//! **AtomPair** (Carhart et al. 1985): encodes pairs of heavy atoms with their
//! shortest-path topological distance and compact atom-type codes.
//!
//! **TopoTorsion** (Nilakantan et al. 1987): encodes four-atom paths that
//! correspond to rotatable torsion angles in the molecular graph.
//!
//! Both fingerprints use FNV-1a 64-bit hashing into a [`BitVec2048`].

use std::collections::{HashSet, VecDeque};

use chematic_core::{AtomIdx, Molecule};
use chematic_perception::find_sssr;

use crate::bitvec::BitVec2048;
use crate::ecfp::fnv1a;

/// Compact atom-type code used for pair/torsion encoding.
///
/// Bit layout: `[atomic_number:8][in_ring:1][aromatic:1]` — 10 bits, fits u16.
fn atom_type(mol: &Molecule, idx: AtomIdx, ring_set: &HashSet<u32>) -> u16 {
    let atom = mol.atom(idx);
    let an = atom.element.atomic_number() as u16;
    let in_ring = ring_set.contains(&idx.0) as u16;
    let aromatic = atom.aromatic as u16;
    (an << 2) | (in_ring << 1) | aromatic
}

/// Collect all atom indices that belong to at least one SSSR ring.
fn ring_atom_set(mol: &Molecule) -> HashSet<u32> {
    let mut set = HashSet::new();
    for ring in find_sssr(mol).rings() {
        for &idx in ring {
            set.insert(idx.0);
        }
    }
    set
}

/// Maximum topological distance tracked for AtomPair.
const MAX_DIST: u8 = 7;

/// Returns `dist[i][j]` = BFS shortest-path distance between atom i and j.
/// Entries are `u8::MAX` for disconnected or distances > `MAX_DIST`.
fn all_pairs_dist(mol: &Molecule) -> Vec<Vec<u8>> {
    let n = mol.atom_count();
    let mut dist = vec![vec![u8::MAX; n]; n];

    for (start, row) in dist.iter_mut().enumerate() {
        let s = AtomIdx(start as u32);
        row[start] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(s);
        while let Some(cur) = queue.pop_front() {
            let d = row[cur.0 as usize];
            if d >= MAX_DIST {
                continue;
            }
            for (nb, _) in mol.neighbors(cur) {
                let ni = nb.0 as usize;
                if row[ni] == u8::MAX {
                    row[ni] = d + 1;
                    queue.push_back(nb);
                }
            }
        }
    }
    dist
}

/// Compute the AtomPair fingerprint (2048-bit).
///
/// For each pair of heavy atoms `(i, j)` with `i < j` and shortest-path
/// distance `d ∈ [1, MAX_DIST]`, hashes `(type_i, d, type_j)` (types ordered
/// canonically) and sets the corresponding bit.
pub fn atom_pair_fp(mol: &Molecule) -> BitVec2048 {
    let ring_set = ring_atom_set(mol);
    let dist = all_pairs_dist(mol);
    let mut fp = BitVec2048::new();

    for (i, row) in dist.iter().enumerate() {
        for (j, &d) in row.iter().enumerate().skip(i + 1) {
            if d == u8::MAX {
                continue;
            }
            let ai = AtomIdx(i as u32);
            let aj = AtomIdx(j as u32);
            let ti = atom_type(mol, ai, &ring_set);
            let tj = atom_type(mol, aj, &ring_set);
            // Canonical: smaller type first.
            let (ta, tb) = if ti <= tj { (ti, tj) } else { (tj, ti) };
            let bytes: [u8; 5] = [(ta >> 8) as u8, ta as u8, d, (tb >> 8) as u8, tb as u8];
            let h = fnv1a(&bytes);
            fp.set((h % 2048) as usize);
        }
    }
    fp
}

/// Compute the Topological Torsion fingerprint (2048-bit).
///
/// Enumerates all four-atom paths `(a–b–c–d)` where `b` and `c` both have
/// degree ≥ 2.  The path is canonicalised as `min(fwd, rev)` and hashed.
pub fn torsion_fp(mol: &Molecule) -> BitVec2048 {
    let n = mol.atom_count();
    let ring_set = ring_atom_set(mol);
    let mut fp = BitVec2048::new();

    // Build adjacency as Vec<Vec<AtomIdx>> for fast access.
    let adj: Vec<Vec<AtomIdx>> = (0..n)
        .map(|i| mol.neighbors(AtomIdx(i as u32)).map(|(nb, _)| nb).collect())
        .collect();

    for b_i in 0..n {
        let b = AtomIdx(b_i as u32);
        if adj[b_i].len() < 2 {
            continue; // b must have degree ≥ 2
        }
        for &c in &adj[b_i] {
            let c_i = c.0 as usize;
            // Avoid double-counting b–c vs c–b.
            if c_i <= b_i {
                continue;
            }
            if adj[c_i].len() < 2 {
                continue; // c must have degree ≥ 2
            }
            // Enumerate (a, d): a ≠ c, d ≠ b, d ≠ a.
            for &a in &adj[b_i] {
                if a == c {
                    continue;
                }
                for &d in &adj[c_i] {
                    if d == b || d == a {
                        continue;
                    }
                    let ta = atom_type(mol, a, &ring_set);
                    let tb = atom_type(mol, b, &ring_set);
                    let tc = atom_type(mol, c, &ring_set);
                    let td = atom_type(mol, d, &ring_set);

                    // Canonicalise direction: pick lexicographically smaller of
                    // (ta,tb,tc,td) and (td,tc,tb,ta).
                    let fwd = [ta, tb, tc, td];
                    let rev = [td, tc, tb, ta];
                    let canonical = if fwd <= rev { fwd } else { rev };

                    let mut bytes = [0u8; 8];
                    for (k, &v) in canonical.iter().enumerate() {
                        bytes[2 * k] = (v >> 8) as u8;
                        bytes[2 * k + 1] = v as u8;
                    }
                    let h = fnv1a(&bytes);
                    fp.set((h % 2048) as usize);
                }
            }
        }
    }
    fp
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> chematic_core::Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    // --- AtomPair ---

    #[test]
    fn test_atom_pair_benzene_nonzero() {
        let fp = atom_pair_fp(&mol("c1ccccc1"));
        assert!(fp.popcount() > 0, "benzene AtomPair FP should be non-zero");
    }

    #[test]
    fn test_atom_pair_identical() {
        let m = mol("CCO");
        assert_eq!(
            atom_pair_fp(&m).tanimoto(&atom_pair_fp(&m)),
            1.0,
            "same molecule → same FP"
        );
    }

    #[test]
    fn test_atom_pair_benzene_vs_cyclohexane() {
        let b = atom_pair_fp(&mol("c1ccccc1"));
        let c = atom_pair_fp(&mol("C1CCCCC1"));
        assert!(
            b.tanimoto(&c) < 1.0,
            "benzene and cyclohexane AtomPair FPs should differ"
        );
    }

    #[test]
    fn test_atom_pair_single_atom_zero() {
        // Single atom: no pairs.
        let fp = atom_pair_fp(&mol("C"));
        assert_eq!(fp.popcount(), 0, "single atom → no pairs");
    }

    // --- Torsion ---

    #[test]
    fn test_torsion_benzene_nonzero() {
        let fp = torsion_fp(&mol("c1ccccc1"));
        assert!(fp.popcount() > 0, "benzene torsion FP should be non-zero");
    }

    #[test]
    fn test_torsion_identical() {
        let m = mol("CCCC");
        assert_eq!(
            torsion_fp(&m).tanimoto(&torsion_fp(&m)),
            1.0,
            "same molecule → same torsion FP"
        );
    }

    #[test]
    fn test_torsion_ethane_zero() {
        // CC: only 2 atoms, no 4-atom paths.
        let fp = torsion_fp(&mol("CC"));
        assert_eq!(fp.popcount(), 0, "ethane has no 4-atom torsion paths");
    }

    #[test]
    fn test_torsion_propane_zero() {
        // CCC: 3 atoms, no 4-atom paths.
        let fp = torsion_fp(&mol("CCC"));
        assert_eq!(fp.popcount(), 0, "propane has no 4-atom torsion paths");
    }

    #[test]
    fn test_torsion_butane_nonzero() {
        // CCCC: exactly one 4-atom path.
        let fp = torsion_fp(&mol("CCCC"));
        assert!(
            fp.popcount() > 0,
            "butane should have at least one torsion path"
        );
    }
}
