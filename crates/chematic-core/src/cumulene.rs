//! OpenSMILES directional geometry on a linear odd cumulene.

use crate::ez_markers::{EzFactEnd, marker_side};
use crate::{AtomIdx, BondIdx, BondOrder, Molecule};

/// One configured chain of at least three consecutive double bonds.
/// `atoms` lists the whole axis from one terminal atom to the other.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlanarCumulene {
    pub atoms: Vec<AtomIdx>,
    pub middle_bond: BondIdx,
    pub ends: [EzFactEnd; 2],
}

fn double_neighbors(mol: &Molecule, atom: AtomIdx) -> impl Iterator<Item = AtomIdx> + '_ {
    mol.neighbors(atom)
        .filter(|&(_, bond)| mol.bond(bond).order == BondOrder::Double)
        .map(|(neighbor, _)| neighbor)
}

fn marked_end(mol: &Molecule, atom: AtomIdx) -> Option<EzFactEnd> {
    let mut substituents = mol
        .neighbors(atom)
        .filter(|&(_, bond)| mol.bond(bond).order != BondOrder::Double);
    let side = |(other, reference)| -> Option<(EzFactEnd, bool)> {
        let fact = EzFactEnd {
            atom,
            reference,
            up: marker_side(mol, reference, atom)?,
        };
        Some((fact, double_neighbors(mol, other).next().is_some()))
    };
    let left = substituents.next().and_then(side);
    let right = substituents.next().and_then(side);
    if substituents.next().is_some() {
        return None;
    }
    match (left, right) {
        (Some((a, _)), None) | (None, Some((a, _))) => Some(a),
        (Some((a, false)), Some((_, true))) => Some(a),
        (Some((_, true)), Some((b, false))) => Some(b),
        (Some((a, _)), Some((b, _))) if a.up != b.up => Some(a),
        _ => None,
    }
}

/// Read a planar cumulene starting at a terminal atom. A terminal has one
/// double-bond neighbor, and each internal atom has exactly two axis bonds.
pub fn planar_cumulene_at(mol: &Molecule, terminal: AtomIdx) -> Option<PlanarCumulene> {
    let mut neighbors = double_neighbors(mol, terminal);
    let adjacent = neighbors.next()?;
    if neighbors.next().is_some()
        || mol.degree(adjacent) != 2
        || double_neighbors(mol, adjacent).count() != 2
    {
        return None;
    }
    let first = marked_end(mol, terminal)?;
    let mut atoms = vec![terminal, adjacent];
    loop {
        let current = atoms[atoms.len() - 1];
        let previous = atoms[atoms.len() - 2];
        let mut onward = double_neighbors(mol, current).filter(|&atom| atom != previous);
        let Some(next) = onward.next() else { break };
        if mol.degree(current) != 2 || onward.next().is_some() {
            return None;
        }
        atoms.push(next);
    }
    if atoms.len() < 4 || !atoms.len().is_multiple_of(2) {
        return None;
    }
    let last = marked_end(mol, *atoms.last()?)?;
    let middle = atoms.len() / 2;
    let middle_bond = mol.bond_between(atoms[middle - 1], atoms[middle])?.0;
    Some(PlanarCumulene {
        atoms,
        middle_bond,
        ends: [first, last],
    })
}

/// Configured odd cumulenes in bond order. Each axis is emitted once.
pub fn planar_cumulenes(mol: &Molecule) -> Vec<PlanarCumulene> {
    let mut result = Vec::new();
    let stashed = mol.has_bond_directions();
    for (idx, bond) in mol.bonds() {
        let marked = matches!(bond.order, BondOrder::Up | BondOrder::Down)
            || stashed && mol.bond_direction(idx).is_some();
        if !marked || bond.order == BondOrder::Double {
            continue;
        }
        for terminal in [bond.atom1, bond.atom2] {
            if let Some(axis) = planar_cumulene_at(mol, terminal)
                && terminal < axis.ends[1].atom
                && axis.ends[0].reference == idx
            {
                result.push(axis);
            }
        }
    }
    result.sort_by_key(|axis| axis.middle_bond);
    result
}
