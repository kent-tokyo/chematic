//! Molecule graph: atoms, bonds, and adjacency list.

use crate::atom::Atom;
use crate::bond::{BondEntry, BondOrder};

/// Newtype index for an atom in a Molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AtomIdx(pub u32);

/// Newtype index for a bond in a Molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct BondIdx(pub u32);

/// Error types for molecule construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MolError {
    /// Atom index out of range.
    InvalidAtomIdx(AtomIdx),
    /// Duplicate bond between the same pair of atoms.
    DuplicateBond(AtomIdx, AtomIdx),
}

impl core::fmt::Display for MolError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidAtomIdx(idx) => write!(f, "invalid atom index: {}", idx.0),
            Self::DuplicateBond(a, b) => write!(f, "duplicate bond between atoms {} and {}", a.0, b.0),
        }
    }
}

/// An immutable molecular graph built via [`MoleculeBuilder`].
///
/// Representation: atom list + bond list + per-atom adjacency list.
/// No external graph library is used; all graph traversal is domain-aware.
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<BondEntry>,
    /// adjacency[atom_idx] = list of (neighbor_atom_idx, bond_idx)
    adjacency: Vec<Vec<(AtomIdx, BondIdx)>>,
}

impl Molecule {
    /// Number of heavy atoms (does not count implicit H).
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Number of bonds (edges).
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Borrow atom by index.
    ///
    /// # Panics
    /// Panics if `idx` is out of range (should not happen with indices from this molecule).
    pub fn atom(&self, idx: AtomIdx) -> &Atom {
        &self.atoms[idx.0 as usize]
    }

    /// Borrow bond by index.
    pub fn bond(&self, idx: BondIdx) -> &BondEntry {
        &self.bonds[idx.0 as usize]
    }

    /// Iterate over all atoms as `(AtomIdx, &Atom)`.
    pub fn atoms(&self) -> impl Iterator<Item = (AtomIdx, &Atom)> {
        self.atoms
            .iter()
            .enumerate()
            .map(|(i, a)| (AtomIdx(i as u32), a))
    }

    /// Iterate over all bonds as `(BondIdx, &BondEntry)`.
    pub fn bonds(&self) -> impl Iterator<Item = (BondIdx, &BondEntry)> {
        self.bonds
            .iter()
            .enumerate()
            .map(|(i, b)| (BondIdx(i as u32), b))
    }

    /// Iterate over neighbors of `idx` as `(neighbor_atom_idx, bond_idx)`.
    pub fn neighbors(&self, idx: AtomIdx) -> impl Iterator<Item = (AtomIdx, BondIdx)> + '_ {
        self.adjacency[idx.0 as usize].iter().copied()
    }

    /// Degree (number of connected bonds) of atom `idx`.
    pub fn degree(&self, idx: AtomIdx) -> usize {
        self.adjacency[idx.0 as usize].len()
    }

    /// Return the bond between `a` and `b`, or `None` if not connected.
    pub fn bond_between(&self, a: AtomIdx, b: AtomIdx) -> Option<(BondIdx, &BondEntry)> {
        self.adjacency[a.0 as usize]
            .iter()
            .find(|&&(nb, _)| nb == b)
            .map(|&(_, bidx)| (bidx, &self.bonds[bidx.0 as usize]))
    }

    /// Molecular formula as a Hill-order string (C first, H second, then alphabetical).
    pub fn formula(&self) -> String {
        use std::collections::BTreeMap;
        let mut counts: BTreeMap<&str, u32> = BTreeMap::new();
        for (_, atom) in self.atoms() {
            *counts.entry(atom.element.symbol()).or_insert(0) += 1;
        }

        let mut result = String::new();
        let push_count = |sym: &str, n: u32, out: &mut String| {
            out.push_str(sym);
            if n > 1 {
                out.push_str(&n.to_string());
            }
        };

        // Hill order: C, then H, then the remaining symbols alphabetically.
        if let Some(c) = counts.remove("C") {
            push_count("C", c, &mut result);
        }
        if let Some(h) = counts.remove("H") {
            push_count("H", h, &mut result);
        }
        for (sym, count) in &counts {
            push_count(sym, *count, &mut result);
        }
        result
    }
}

/// Builder for constructing a [`Molecule`] incrementally.
///
/// Usage: add atoms, add bonds, then call `build()`.
#[derive(Default)]
pub struct MoleculeBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<BondEntry>,
    adjacency: Vec<Vec<(AtomIdx, BondIdx)>>,
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    /// Read-only reference to an atom already added to the builder.
    ///
    /// Used by the SMILES parser to infer implicit bond types without
    /// consuming the builder (e.g. aromatic-aromatic → Aromatic bond).
    ///
    /// # Panics
    /// Panics if `idx` is out of range.
    pub fn atom_at(&self, idx: AtomIdx) -> &Atom {
        &self.atoms[idx.0 as usize]
    }

    /// Number of atoms added so far.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Iterate over already-added neighbors of `idx` as `(bond_idx, neighbor_atom_idx)`.
    /// Used by kekulization tests to check whether a bond already exists in the builder.
    pub fn atom_neighbors(&self, idx: AtomIdx) -> impl Iterator<Item = (BondIdx, AtomIdx)> + '_ {
        self.adjacency[idx.0 as usize].iter().map(|&(nb, bidx)| (bidx, nb))
    }

    /// Add an atom and return its index.
    pub fn add_atom(&mut self, atom: Atom) -> AtomIdx {
        let idx = AtomIdx(self.atoms.len() as u32);
        self.atoms.push(atom);
        self.adjacency.push(Vec::new());
        idx
    }

    /// Add a bond between two existing atoms.
    ///
    /// Returns an error if either atom index is invalid or if the bond already exists.
    pub fn add_bond(&mut self, a: AtomIdx, b: AtomIdx, order: BondOrder) -> Result<BondIdx, MolError> {
        let n = self.atoms.len() as u32;
        if a.0 >= n { return Err(MolError::InvalidAtomIdx(a)); }
        if b.0 >= n { return Err(MolError::InvalidAtomIdx(b)); }

        // Check for duplicate
        for &(nb, _) in &self.adjacency[a.0 as usize] {
            if nb == b {
                return Err(MolError::DuplicateBond(a, b));
            }
        }

        let bidx = BondIdx(self.bonds.len() as u32);
        self.bonds.push(BondEntry { atom1: a, atom2: b, order });
        self.adjacency[a.0 as usize].push((b, bidx));
        self.adjacency[b.0 as usize].push((a, bidx));
        Ok(bidx)
    }

    /// Consume the builder and return an immutable [`Molecule`].
    pub fn build(self) -> Molecule {
        Molecule {
            atoms: self.atoms,
            bonds: self.bonds,
            adjacency: self.adjacency,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::element::Element;

    fn ethane() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn test_basic_molecule() {
        let mol = ethane();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn test_adjacency() {
        let mol = ethane();
        let neighbors: Vec<_> = mol.neighbors(AtomIdx(0)).collect();
        assert_eq!(neighbors.len(), 1);
        assert_eq!(neighbors[0].0, AtomIdx(1));
    }

    #[test]
    fn test_bond_between() {
        let mol = ethane();
        assert!(mol.bond_between(AtomIdx(0), AtomIdx(1)).is_some());
        assert!(mol.bond_between(AtomIdx(1), AtomIdx(0)).is_some());
    }

    #[test]
    fn test_duplicate_bond_error() {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let err = b.add_bond(c1, c2, BondOrder::Double);
        assert!(matches!(err, Err(MolError::DuplicateBond(_, _))));
    }

    #[test]
    fn test_formula() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let n = b.add_atom(Atom::new(Element::N));
        b.add_bond(c, n, BondOrder::Single).unwrap();
        let mol = b.build();
        assert_eq!(mol.formula(), "CN");
    }
}
