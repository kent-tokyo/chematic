//! Detect a charge-explicit, two-ring haptic encoding before metal disconnection.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondOrder, Molecule};

use super::is_metal;

pub(super) fn specified_haptic_iron_centers(mol: &Molecule) -> HashSet<usize> {
    mol.atoms()
        .filter(|&(idx, atom)| {
            atom.element.atomic_number() == 26
                && atom.charge == 2
                && mol.degree(idx) == 10
                && attaches_two_cyclopentadienides(mol, idx)
        })
        .map(|(idx, _)| idx.0 as usize)
        .collect()
}

fn attaches_two_cyclopentadienides(mol: &Molecule, iron: AtomIdx) -> bool {
    let mut bound = HashSet::with_capacity(10);
    for (ligand, bond) in mol.neighbors(iron) {
        let atom = mol.atom(ligand);
        if atom.wildcard
            || atom.element.atomic_number() != 6
            || !matches!(mol.bond(bond).order, BondOrder::Single | BondOrder::Dative)
        {
            return false;
        }
        bound.insert(ligand);
    }
    let mut remaining = bound;
    for _ in 0..2 {
        let Some(&start) = remaining.iter().next() else {
            return false;
        };
        let mut pending = vec![start];
        let mut ring = Vec::with_capacity(5);
        while let Some(idx) = pending.pop() {
            if remaining.remove(&idx) {
                ring.push(idx);
                pending.extend(
                    mol.neighbors(idx).filter_map(|(neighbor, _)| {
                        remaining.contains(&neighbor).then_some(neighbor)
                    }),
                );
            }
        }
        if ring.len() != 5 || !is_cyclopentadienide(mol, iron, &ring) {
            return false;
        }
    }
    remaining.is_empty()
}

fn is_cyclopentadienide(mol: &Molecule, iron: AtomIdx, ring: &[AtomIdx]) -> bool {
    let mut anions = 0;
    for &idx in ring {
        let atom = mol.atom(idx);
        if !matches!(atom.charge, -1 | 0) {
            return false;
        }
        let anionic = atom.charge == -1;
        anions += usize::from(anionic);
        let mut ring_neighbors = 0;
        let mut doubles = 0;
        let mut substituents = 0;
        for (neighbor, bond_idx) in mol.neighbors(idx) {
            if neighbor == iron {
                continue;
            }
            let order = mol.bond(bond_idx).order;
            if ring.contains(&neighbor) {
                ring_neighbors += 1;
                match order {
                    BondOrder::Single => {}
                    BondOrder::Double => doubles += 1,
                    _ => return false,
                }
            } else if order == BondOrder::Single && !is_metal(mol.atom(neighbor).element) {
                substituents += 1;
            } else {
                return false;
            }
        }
        if ring_neighbors != 2 || doubles != usize::from(!anionic) || substituents > 1 {
            return false;
        }
        let expected_h = u8::from(substituents == 0);
        if atom
            .hydrogen_count
            .is_some_and(|count| count != expected_h && !(anionic && count == 0))
        {
            return false;
        }
    }
    anions == 1
}
