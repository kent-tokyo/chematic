//! Previous implementation retained only as a seed/atom-membership oracle.
use super::*;
impl Prepared<'_> {
    fn frozen_hashes(&self, path: &[u32]) -> ([Option<u32>; NUM_LAYERS], Vec<usize>) {
        let (degrees, atoms, bond_nbrs) = frozen_path_info(self.mol, path);
        let mut layers: [Vec<u32>; NUM_LAYERS] = Default::default();
        for (i, &bidx) in path.iter().enumerate() {
            let bond = self.mol.bond(BondIdx(bidx));
            let a = bond.atom1.0;
            let b = bond.atom2.0;
            let d1 = degrees[&a];
            let d2 = degrees[&b];
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
        (seeds, atoms)
    }
}
type FrozenPathInfo = (FxHashMap<u32, u32>, Vec<usize>, Vec<u32>);
fn frozen_path_info(mol: &Molecule, path: &[u32]) -> FrozenPathInfo {
    let mut degrees = FxHashMap::default();
    let mut atoms = std::collections::BTreeSet::new();
    for &idx in path {
        let b = mol.bond(BondIdx(idx));
        *degrees.entry(b.atom1.0).or_insert(0) += 1;
        *degrees.entry(b.atom2.0).or_insert(0) += 1;
        atoms.insert(b.atom1.0 as usize);
        atoms.insert(b.atom2.0 as usize);
    }
    let mut nbrs = vec![0; path.len()];
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
    (degrees, atoms.into_iter().collect(), nbrs)
}

#[test]
fn inline_layer_storage_keeps_frozen_seeds_and_atom_membership() {
    let corpus = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
    ));
    let mut paths = 0;
    for text in corpus
        .lines()
        .take(256)
        .map(|s| s.split_whitespace().next().unwrap())
        .chain([
            "[2H]CC([3H])O",
            "C12C3C4C1C5C2C3C45",
            "c1ccccc#1",
            "C1=CC=CC=C1",
            "CC(C)(C)CCCCCCCCCCC",
        ])
    {
        let mol = chematic_smiles::parse_template(text).unwrap();
        let prepared = Prepared::new(&mol, false);
        for plist in find_all_subgraphs(&mol, 1, 9).values() {
            for path in plist {
                let candidate = prepared.hashes(path);
                let expected = prepared.frozen_hashes(path);
                assert_eq!(candidate.seeds, expected.0, "{text} {path:?}");
                assert_eq!(candidate.atoms.as_slice(), expected.1, "{text} {path:?}");
                paths += 1;
            }
        }
    }
    assert!(paths > 10000);
    // Configurable paths are not limited to the default seven bonds. Check
    // both inline buffers spilling without imposing a new size limit.
    let mol = chematic_smiles::parse_template("CCCCCCCCCCCCCCCCCCC").unwrap();
    let prepared = Prepared::new(&mol, false);
    let path: Vec<_> = (0..mol.bond_count() as u32).collect();
    let candidate = prepared.hashes(&path);
    let expected = prepared.frozen_hashes(&path);
    assert_eq!(candidate.seeds, expected.0);
    assert_eq!(candidate.atoms.as_slice(), expected.1);
}
