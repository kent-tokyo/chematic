//! Shared feature hashing, retaining the existing six-layer bit contract.
use super::*;
use smallvec::SmallVec;

pub(super) struct Prepared<'a> {
    mol: &'a Molecule,
    aromaticity: Option<AromaticityModel>,
    anums: Vec<u32>,
    atom_aromatic: Vec<bool>,
    ring_info: RingInfo,
}
pub(super) struct Hashes {
    pub seeds: [Option<u32>; NUM_LAYERS],
    pub atoms: SmallVec<[usize; 14]>,
}
impl<'a> Prepared<'a> {
    pub fn new(mol: &'a Molecule, trust_view: bool) -> Self {
        let has_literal_aromatic_bond =
            trust_view || mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic);
        let aromaticity = if has_literal_aromatic_bond {
            None
        } else {
            Some(assign_aromaticity_ex(mol, AromaticityAlgorithm::RdkitLike))
        };
        let anums = (0..mol.atom_count())
            .map(|i| crate::rdkit_atomic_num(mol, AtomIdx(i as u32)))
            .collect();
        let atom_aromatic = (0..mol.atom_count())
            .map(|i| is_atom_aromatic(mol, aromaticity.as_ref(), AtomIdx(i as u32)))
            .collect();
        Self {
            mol,
            aromaticity,
            anums,
            atom_aromatic,
            ring_info: compute_ring_info(mol),
        }
    }
    pub fn hashes(&self, path: &[u32]) -> Hashes {
        let (degrees, atoms, bond_nbrs) = path_info(self.mol, path);
        // Default paths have at most seven bonds; configured longer paths spill.
        let mut layers: [SmallVec<[u32; 9]>; NUM_LAYERS] = Default::default();
        for (i, &bidx) in path.iter().enumerate() {
            let bond = self.mol.bond(BondIdx(bidx));
            let a = bond.atom1.0;
            let b = bond.atom2.0;
            let degree = |atom| degrees.iter().find(|&&(a, _)| a == atom).unwrap().1;
            let d1 = degree(a);
            let d2 = degree(b);
            let n = bond_nbrs[i];
            layers[0].push(topology(d1, d2, n));
            let code = if is_bond_aromatic(self.mol, self.aromaticity.as_ref(), BondIdx(bidx)) {
                1
            } else {
                layer2_bond_type_code(bond.order)
            };
            layers[1].push(bond_type(d1, d2, n, code));
            layers[2].push(atom_type(
                self.anums[a as usize],
                self.anums[b as usize],
                d1,
                d2,
                n,
            ));
            if self.ring_info.bond_in_ring[bidx as usize] {
                layers[3].push(1);
            }
            layers[4].push(self.ring_info.bond_min_ring_size[bidx as usize] % 8);
            layers[5].push(aromaticity(
                self.atom_aromatic[a as usize],
                self.atom_aromatic[b as usize],
                n,
            ));
        }
        let mut seeds = [None; NUM_LAYERS];
        for (i, layer) in layers.iter_mut().enumerate() {
            if layer.is_empty() {
                continue;
            }
            layer.sort_unstable();
            layer.push(atoms.len() as u32);
            layer.push((i + 1) as u32);
            seeds[i] = Some(layer.iter().fold(0, |h, &v| hash_combine(h, v)));
        }
        Hashes { seeds, atoms }
    }
}
type PathInfo = (
    SmallVec<[(u32, u32); 14]>,
    SmallVec<[usize; 14]>,
    SmallVec<[u32; 7]>,
);
fn path_info(mol: &Molecule, path: &[u32]) -> PathInfo {
    let mut degrees: SmallVec<[(u32, u32); 14]> = SmallVec::new();
    for &idx in path {
        let b = mol.bond(BondIdx(idx));
        for atom in [b.atom1.0, b.atom2.0] {
            if let Some((_, count)) = degrees.iter_mut().find(|(a, _)| *a == atom) {
                *count += 1;
            } else {
                degrees.push((atom, 1));
            }
        }
    }
    let mut nbrs: SmallVec<[u32; 7]> = std::iter::repeat_n(0, path.len()).collect();
    for i in 0..path.len() {
        let a = mol.bond(BondIdx(path[i]));
        for j in i + 1..path.len() {
            let b = mol.bond(BondIdx(path[j]));
            if a.atom1 == b.atom1 || a.atom1 == b.atom2 || a.atom2 == b.atom1 || a.atom2 == b.atom2
            {
                nbrs[i] += 1;
                nbrs[j] += 1;
            }
        }
    }
    let mut atoms: SmallVec<[usize; 14]> = degrees.iter().map(|&(a, _)| a as usize).collect();
    atoms.sort_unstable();
    (degrees, atoms, nbrs)
}
fn topology(mut a: u32, mut b: u32, n: u32) -> u32 {
    if a < b {
        std::mem::swap(&mut a, &mut b);
    }
    (n % 8) | ((a % 8) << 3) | ((b % 8) << 6)
}
fn bond_type(mut a: u32, mut b: u32, n: u32, code: u32) -> u32 {
    if a < b {
        std::mem::swap(&mut a, &mut b);
    }
    (code % 8) | ((n % 8) << 3) | ((a % 8) << 6) | ((b % 8) << 9)
}
fn atom_type(a: u32, b: u32, mut da: u32, mut db: u32, n: u32) -> u32 {
    let mut a = a % 128;
    let mut b = b % 128;
    if a < b || (a == b && da < db) {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut da, &mut db);
    }
    a | (b << 7) | ((da % 8) << 14) | ((db % 8) << 17) | ((n % 8) << 20)
}
fn aromaticity(mut a: bool, mut b: bool, n: u32) -> u32 {
    if !a && b {
        std::mem::swap(&mut a, &mut b);
    }
    (a as u32) | ((b as u32) << 1) | ((n % 8) << 5)
}

#[cfg(test)]
mod oracle_tests;
