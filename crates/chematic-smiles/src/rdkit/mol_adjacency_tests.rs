//! Frozen c84c68f2 Vec adjacency mutations, with unchanged bodies.
use super::{Atom, Bond, BondType, Mol, RingInfo};

#[derive(Clone, Debug, Default)]
struct ReferenceMol {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    atom_bonds: Vec<Vec<usize>>,
    rings: Option<RingInfo>,
}

impl ReferenceMol {
    pub(crate) fn add_atom(&mut self, atom: Atom) -> usize {
        self.atoms.push(atom);
        self.atom_bonds.push(Vec::new());
        self.atoms.len() - 1
    }

    pub(crate) fn add_bond(&mut self, bond: Bond) -> usize {
        let idx = self.bonds.len();
        self.atom_bonds[bond.begin].push(idx);
        self.atom_bonds[bond.end].push(idx);
        self.bonds.push(bond);
        idx
    }

    pub(crate) fn remove_atom(&mut self, idx: usize) {
        let removed: Vec<usize> = self.atom_bonds[idx].clone();
        let mut keep = vec![true; self.bonds.len()];
        for &b in &removed {
            keep[b] = false;
        }
        let mut new_bond_idx = vec![usize::MAX; self.bonds.len()];
        let mut bonds = Vec::with_capacity(self.bonds.len());
        for (b, bond) in self.bonds.drain(..).enumerate() {
            if keep[b] {
                new_bond_idx[b] = bonds.len();
                bonds.push(bond);
            }
        }
        let fix_atom = |x: usize| if x > idx { x - 1 } else { x };
        for bond in &mut bonds {
            bond.begin = fix_atom(bond.begin);
            bond.end = fix_atom(bond.end);
            bond.stereo_atoms = bond.stereo_atoms.iter().map(|&x| fix_atom(x)).collect();
        }
        self.bonds = bonds;
        self.atoms.remove(idx);
        self.atom_bonds.remove(idx);
        for list in &mut self.atom_bonds {
            *list = list
                .iter()
                .filter(|&&b| keep[b])
                .map(|&b| new_bond_idx[b])
                .collect();
        }
        self.rings = None;
    }

    pub(crate) fn remove_bond(&mut self, b: usize) {
        self.bonds.remove(b);
        for list in &mut self.atom_bonds {
            list.retain(|&x| x != b);
            for x in list.iter_mut() {
                if *x > b {
                    *x -= 1;
                }
            }
        }
        self.rings = None;
    }

    pub(crate) fn remove_atoms(&mut self, removed: &[usize]) {
        if removed.is_empty() {
            return;
        }
        let n = self.atoms.len();
        let mut gone = vec![false; n];
        for &a in removed {
            gone[a] = true;
        }
        // Each removal shifts the indices above it down by one; a stereo
        // atom equal to a removed index ends up shifted the same way.
        let mut shift = vec![0usize; n + 1];
        for i in 0..n {
            shift[i + 1] = shift[i] + usize::from(gone[i]);
        }
        let fix_atom = |x: usize| x - shift[x.min(n)];
        let mut new_bond_idx = vec![usize::MAX; self.bonds.len()];
        let mut bonds = Vec::with_capacity(self.bonds.len());
        for (b, mut bond) in self.bonds.drain(..).enumerate() {
            if gone[bond.begin] || gone[bond.end] {
                continue;
            }
            new_bond_idx[b] = bonds.len();
            bond.begin = fix_atom(bond.begin);
            bond.end = fix_atom(bond.end);
            for x in &mut bond.stereo_atoms {
                *x = fix_atom(*x);
            }
            bonds.push(bond);
        }
        self.bonds = bonds;
        let mut i = 0;
        self.atoms.retain(|_| {
            i += 1;
            !gone[i - 1]
        });
        let mut i = 0;
        self.atom_bonds.retain(|_| {
            i += 1;
            !gone[i - 1]
        });
        for list in &mut self.atom_bonds {
            list.retain(|&b| new_bond_idx[b] != usize::MAX);
            for b in list.iter_mut() {
                *b = new_bond_idx[*b];
            }
        }
        self.rings = None;
    }
}

fn assert_same(actual: &Mol, expected: &ReferenceMol) {
    assert_eq!(
        format!("{:?}", actual.atoms),
        format!("{:?}", expected.atoms)
    );
    assert_eq!(
        format!("{:?}", actual.bonds),
        format!("{:?}", expected.bonds)
    );
    let adjacency: Vec<Vec<usize>> = actual.atom_bonds.iter().map(|b| b.to_vec()).collect();
    assert_eq!(adjacency, expected.atom_bonds);
    assert_eq!(
        format!("{:?}", actual.rings),
        format!("{:?}", expected.rings)
    );
    for (a, bonds) in actual.atom_bonds.iter().enumerate() {
        assert!(bonds.windows(2).all(|w| w[0] <= w[1]));
        for &b in bonds {
            assert!(actual.bonds[b].begin == a || actual.bonds[b].end == a);
        }
    }
}

fn graph(n: usize, seed: usize) -> (Mol, ReferenceMol) {
    let mut actual = Mol::default();
    let mut expected = ReferenceMol::default();
    for a in 0..n {
        let mut atom = Atom::new(if a % 3 == 0 { 7 } else { 6 });
        atom.isotope = (a % 4) as u32;
        actual.add_atom(atom.clone());
        expected.add_atom(atom);
    }
    for a in 0..n {
        for b in a..n {
            if !(a * 17 + b * 7 + seed).is_multiple_of(5) && a != 0 {
                continue;
            }
            let mut bond = Bond::new(a, b, BondType::Single);
            bond.stereo_atoms = vec![a, b];
            actual.add_bond(bond.clone());
            expected.add_bond(bond);
        }
    }
    actual.rings = Some(RingInfo::default());
    expected.rings = Some(RingInfo::default());
    (actual, expected)
}

#[test]
fn inline_spill_clone_and_frozen_mutations_preserve_every_index() {
    for n in 1..25 {
        for seed in 0..7 {
            let (actual, expected) = graph(n, seed);
            assert_same(&actual, &expected);
            for (a, bonds) in actual.atom_bonds.iter().enumerate() {
                assert_eq!(bonds.spilled(), expected.atom_bonds[a].len() > 4);
            }
            for a in 0..n {
                let (mut x, mut y) = (actual.clone(), expected.clone());
                x.remove_atom(a);
                y.remove_atom(a);
                assert_same(&x, &y);
            }
            for b in 0..actual.bonds.len() {
                let (mut x, mut y) = (actual.clone(), expected.clone());
                x.remove_bond(b);
                y.remove_bond(b);
                assert_same(&x, &y);
            }
            for step in 1..5 {
                let (mut x, mut y) = (actual.clone(), expected.clone());
                let removed: Vec<_> = (0..n).filter(|a| a % step == 0).collect();
                x.remove_atoms(&removed);
                y.remove_atoms(&removed);
                assert_same(&x, &y);
                let atom = Atom::new(8);
                x.add_atom(atom.clone());
                y.add_atom(atom);
                if x.atoms.len() > 1 {
                    let bond = Bond::new(0, x.atoms.len() - 1, BondType::Single);
                    x.add_bond(bond.clone());
                    y.add_bond(bond);
                }
                assert_same(&x.clone(), &y.clone());
            }
            let (mut x, mut y) = (actual.clone(), expected.clone());
            x.remove_atoms(&[]);
            y.remove_atoms(&[]);
            assert_same(&x, &y);
            while !x.atoms.is_empty() {
                let a = (seed + x.atoms.len()) % x.atoms.len();
                x.remove_atom(a);
                y.remove_atom(a);
                assert_same(&x, &y);
            }
        }
    }
}
