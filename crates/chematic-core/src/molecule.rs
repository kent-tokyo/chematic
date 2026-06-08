//! Molecule graph: atoms, bonds, and adjacency list.

use crate::atom::Atom;
use crate::bond::{BondEntry, BondOrder};
use crate::element::Element;
use crate::stereo_group::StereoGroup;

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
            Self::DuplicateBond(a, b) => {
                write!(f, "duplicate bond between atoms {} and {}", a.0, b.0)
            }
        }
    }
}

impl std::error::Error for MolError {}

/// An immutable molecular graph built via [`MoleculeBuilder`].
///
/// Representation: atom list + bond list + per-atom adjacency list.
/// No external graph library is used; all graph traversal is domain-aware.
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<BondEntry>,
    /// adjacency[atom_idx] = list of (neighbor_atom_idx, bond_idx)
    adjacency: Vec<Vec<(AtomIdx, BondIdx)>>,
    /// Enhanced stereo groups (ChemDraw V3000 Absolute / Or / And).
    stereo_groups: Vec<StereoGroup>,
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

// ---------------------------------------------------------------------------
// Immutable update methods (functional-style editing)
// ---------------------------------------------------------------------------

impl Molecule {
    /// Return a new `Molecule` with one extra atom appended, along with the
    /// index that the new atom will have in the returned molecule.
    pub fn with_atom_added(&self, atom: Atom) -> (Molecule, AtomIdx) {
        let mut builder = MoleculeBuilder::new();
        for (_, a) in self.atoms() {
            builder.add_atom(a.clone());
        }
        for (_, b) in self.bonds() {
            let _ = builder.add_bond(b.atom1, b.atom2, b.order);
        }
        let new_idx = builder.add_atom(atom);
        (builder.build(), new_idx)
    }

    /// Return a new `Molecule` with one extra bond added, along with the index
    /// of the newly added bond in the returned molecule.
    ///
    /// Returns `Err` if `a == b` or the bond already exists (same semantics as
    /// [`MoleculeBuilder::add_bond`]).
    pub fn with_bond_added(
        &self,
        a: AtomIdx,
        b: AtomIdx,
        order: BondOrder,
    ) -> Result<(Molecule, BondIdx), MolError> {
        let mut builder = MoleculeBuilder::new();
        for (_, atom) in self.atoms() {
            builder.add_atom(atom.clone());
        }
        for (_, bond) in self.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        let bond_idx = builder.add_bond(a, b, order)?;
        Ok((builder.build(), bond_idx))
    }

    /// Return a new `Molecule` with the formal charge of atom `idx` changed.
    pub fn with_atom_charge(&self, idx: AtomIdx, charge: i8) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        for (aidx, atom) in self.atoms() {
            let mut a = atom.clone();
            if aidx == idx {
                a.charge = charge;
            }
            builder.add_atom(a);
        }
        for (_, bond) in self.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        builder.build()
    }

    /// Return a new `Molecule` with the element of atom `idx` changed.
    ///
    /// Chirality and hydrogen count are reset to `None` when the element
    /// changes, since those properties are element-specific.
    pub fn with_atom_element(&self, idx: AtomIdx, el: Element) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        for (aidx, atom) in self.atoms() {
            let mut a = atom.clone();
            if aidx == idx {
                a.element = el;
                // Reset element-specific fields so valence stays consistent.
                a.chirality = crate::atom::Chirality::None;
                a.hydrogen_count = None;
                a.aromatic = false;
            }
            builder.add_atom(a);
        }
        for (_, bond) in self.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        builder.build()
    }

    /// Return a new `Molecule` with atom `idx` and all bonds involving it
    /// removed.  Atom indices of survivors shift down past the removed slot.
    ///
    /// The returned tuple also includes a mapping from **old** `AtomIdx` to
    /// **new** `AtomIdx` (indices that fall below `idx` are unchanged; indices
    /// above `idx` decrease by 1).
    pub fn with_atom_removed(&self, idx: AtomIdx) -> (Molecule, Vec<Option<AtomIdx>>) {
        let n = self.atom_count();
        let removed = idx.0 as usize;

        // Build old→new index table.
        let mut remap: Vec<Option<AtomIdx>> = vec![None; n];
        let mut new_pos = 0u32;
        for (old, slot) in remap.iter_mut().enumerate() {
            if old == removed {
                continue;
            }
            *slot = Some(AtomIdx(new_pos));
            new_pos += 1;
        }

        let mut builder = MoleculeBuilder::new();
        for (aidx, atom) in self.atoms() {
            if aidx == idx {
                continue;
            }
            builder.add_atom(atom.clone());
        }
        for (_, bond) in self.bonds() {
            if bond.atom1 == idx || bond.atom2 == idx {
                continue;
            }
            if let (Some(a1), Some(a2)) =
                (remap[bond.atom1.0 as usize], remap[bond.atom2.0 as usize])
            {
                let _ = builder.add_bond(a1, a2, bond.order);
            }
        }
        (builder.build(), remap)
    }

    /// Implicit hydrogen count for atom `idx` based on valence rules.
    ///
    /// Delegates to [`crate::valence::implicit_hcount`].
    pub fn implicit_hydrogen_count(&self, idx: AtomIdx) -> u8 {
        crate::valence::implicit_hcount(self, idx)
    }

    /// Hill-order molecular formula including implicit hydrogens.
    ///
    /// Unlike [`Self::formula`] (which counts only explicit heavy atoms),
    /// this method adds the implicit H count for every atom so the result
    /// reflects the true molecular composition (e.g. methane → "CH4").
    pub fn total_formula(&self) -> String {
        use std::collections::BTreeMap;
        let mut counts: BTreeMap<&str, u32> = BTreeMap::new();
        let mut implicit_h: u32 = 0;
        for (aidx, atom) in self.atoms() {
            *counts.entry(atom.element.symbol()).or_insert(0) += 1;
            implicit_h += crate::valence::implicit_hcount(self, aidx) as u32;
        }
        *counts.entry("H").or_insert(0) += implicit_h;

        let mut result = String::new();
        let push_count = |sym: &str, n: u32, out: &mut String| {
            out.push_str(sym);
            if n > 1 {
                out.push_str(&n.to_string());
            }
        };

        if let Some(c) = counts.remove("C") {
            push_count("C", c, &mut result);
        }
        if let Some(h) = counts.remove("H")
            && h > 0
        {
            push_count("H", h, &mut result);
        }
        for (sym, count) in &counts {
            push_count(sym, *count, &mut result);
        }
        result
    }

    /// Hill-order molecular formula with isotope labels.
    ///
    /// Like [`Self::formula`] but prefixes each element symbol with its
    /// isotope number when `atom.isotope` is `Some(n)`.
    /// Example: a molecule with one `¹³C` and one `O` → `"¹³CO"`.
    pub fn formula_with_isotopes(&self) -> String {
        use std::collections::BTreeMap;
        // Collect (isotope_prefix + symbol) counts, heavy atoms only.
        let mut counts: BTreeMap<String, u32> = BTreeMap::new();
        let mut has_carbon = false;
        let mut has_explicit_h = false;
        for (_, atom) in self.atoms() {
            let sym = atom.element.symbol();
            let key = match atom.isotope {
                Some(n) => format!("{n}{sym}"),
                None => sym.to_string(),
            };
            if sym == "C" && atom.isotope.is_none() {
                has_carbon = true;
            }
            if sym == "H" {
                has_explicit_h = true;
            }
            *counts.entry(key).or_insert(0) += 1;
        }

        let push_count = |key: &str, n: u32, out: &mut String| {
            out.push_str(key);
            if n > 1 {
                out.push_str(&n.to_string());
            }
        };

        let mut result = String::new();
        // Hill order: C first (if unlabelled C present), then H, then rest alphabetically.
        if has_carbon && let Some(c) = counts.remove("C") {
            push_count("C", c, &mut result);
        }
        if has_explicit_h && let Some(h) = counts.remove("H") {
            push_count("H", h, &mut result);
        }
        for (key, count) in &counts {
            push_count(key, *count, &mut result);
        }
        result
    }

    /// Return a new `Molecule` with atom `idx`'s aromatic flag changed.
    pub fn with_atom_aromatic(&self, idx: AtomIdx, aromatic: bool) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        for (aidx, atom) in self.atoms() {
            let mut a = atom.clone();
            if aidx == idx {
                a.aromatic = aromatic;
            }
            builder.add_atom(a);
        }
        for (_, bond) in self.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        builder.build()
    }

    /// Return a new `Molecule` with bond `idx`'s order changed.
    pub fn with_bond_order(&self, idx: BondIdx, order: BondOrder) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        for (_, atom) in self.atoms() {
            builder.add_atom(atom.clone());
        }
        for (bidx, bond) in self.bonds() {
            let o = if bidx == idx { order } else { bond.order };
            let _ = builder.add_bond(bond.atom1, bond.atom2, o);
        }
        builder.build()
    }

    /// Return a new `Molecule` with bond `idx` removed.
    ///
    /// Atom indices are unchanged.  Bond indices of survivors shift down.
    pub fn with_bond_removed(&self, idx: BondIdx) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        for (_, atom) in self.atoms() {
            builder.add_atom(atom.clone());
        }
        for (bidx, bond) in self.bonds() {
            if bidx == idx {
                continue;
            }
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        builder.build()
    }
}

// ---------------------------------------------------------------------------
// In-place mutation methods
// ---------------------------------------------------------------------------

impl Molecule {
    /// Append a new atom and return its index.
    pub fn add_atom(&mut self, atom: Atom) -> AtomIdx {
        let idx = AtomIdx(self.atoms.len() as u32);
        self.atoms.push(atom);
        self.adjacency.push(vec![]);
        idx
    }

    /// Remove atom `idx` and all bonds involving it.
    ///
    /// Returns a remapping table: `remap[old_idx]` gives the new `AtomIdx`
    /// for surviving atoms, or `None` for the removed atom.  Atom indices
    /// of atoms after the removed slot shift down by 1.
    pub fn remove_atom(&mut self, idx: AtomIdx) -> Vec<Option<AtomIdx>> {
        let n = self.atoms.len();
        let removed = idx.0 as usize;

        let mut remap: Vec<Option<AtomIdx>> = vec![None; n];
        let mut new_pos = 0u32;
        for (old, slot) in remap.iter_mut().enumerate() {
            if old == removed {
                continue;
            }
            *slot = Some(AtomIdx(new_pos));
            new_pos += 1;
        }

        self.atoms.remove(removed);

        // Keep only bonds not involving the removed atom; remap endpoints.
        let mut new_bonds: Vec<BondEntry> = Vec::new();
        for bond in &self.bonds {
            if bond.atom1 == idx || bond.atom2 == idx {
                continue;
            }
            if let (Some(a1), Some(a2)) =
                (remap[bond.atom1.0 as usize], remap[bond.atom2.0 as usize])
            {
                new_bonds.push(BondEntry {
                    atom1: a1,
                    atom2: a2,
                    order: bond.order,
                });
            }
        }
        self.bonds = new_bonds;

        // Rebuild adjacency from scratch.
        let new_n = self.atoms.len();
        self.adjacency = vec![vec![]; new_n];
        for (bidx, bond) in self.bonds.iter().enumerate() {
            let bi = BondIdx(bidx as u32);
            self.adjacency[bond.atom1.0 as usize].push((bond.atom2, bi));
            self.adjacency[bond.atom2.0 as usize].push((bond.atom1, bi));
        }

        remap
    }

    /// Add a bond between `a` and `b` with the given `order`.
    ///
    /// Returns `Err` if `a == b` or the bond already exists.
    pub fn add_bond(
        &mut self,
        a: AtomIdx,
        b: AtomIdx,
        order: BondOrder,
    ) -> Result<BondIdx, MolError> {
        let n = self.atoms.len() as u32;
        if a.0 >= n {
            return Err(MolError::InvalidAtomIdx(a));
        }
        if b.0 >= n {
            return Err(MolError::InvalidAtomIdx(b));
        }
        if self.adjacency[a.0 as usize].iter().any(|&(nb, _)| nb == b) {
            return Err(MolError::DuplicateBond(a, b));
        }
        let bidx = BondIdx(self.bonds.len() as u32);
        self.bonds.push(BondEntry {
            atom1: a,
            atom2: b,
            order,
        });
        self.adjacency[a.0 as usize].push((b, bidx));
        self.adjacency[b.0 as usize].push((a, bidx));
        Ok(bidx)
    }

    /// Remove bond `idx`.  Atom indices are unchanged; bond indices of
    /// surviving bonds shift down past the removed slot.
    pub fn remove_bond(&mut self, idx: BondIdx) {
        let removed = idx.0 as usize;
        if removed >= self.bonds.len() {
            return;
        }
        self.bonds.remove(removed);
        // Rebuild adjacency with renumbered bond indices.
        let n = self.atoms.len();
        self.adjacency = vec![vec![]; n];
        for (bidx, bond) in self.bonds.iter().enumerate() {
            let bi = BondIdx(bidx as u32);
            self.adjacency[bond.atom1.0 as usize].push((bond.atom2, bi));
            self.adjacency[bond.atom2.0 as usize].push((bond.atom1, bi));
        }
    }

    /// Set the formal charge of atom `idx` in-place.
    pub fn set_charge(&mut self, idx: AtomIdx, charge: i8) {
        self.atoms[idx.0 as usize].charge = charge;
    }

    /// Set the element of atom `idx` in-place.
    ///
    /// Chirality and hydrogen count are reset (element-specific properties).
    pub fn set_element(&mut self, idx: AtomIdx, el: Element) {
        let a = &mut self.atoms[idx.0 as usize];
        a.element = el;
        a.chirality = crate::atom::Chirality::None;
        a.hydrogen_count = None;
        a.aromatic = false;
    }

    /// Set the CIP stereo code of atom `idx` in-place.
    pub fn set_cip_code(&mut self, idx: AtomIdx, code: Option<crate::atom::CipCode>) {
        self.atoms[idx.0 as usize].cip_code = code;
    }

    /// Return the enhanced stereo groups attached to this molecule.
    pub fn stereo_groups(&self) -> &[StereoGroup] {
        &self.stereo_groups
    }

    /// Replace the stereo group list in-place.
    pub fn set_stereo_groups(&mut self, groups: Vec<StereoGroup>) {
        self.stereo_groups = groups;
    }

    /// Add a single stereo group in-place.
    pub fn add_stereo_group(&mut self, group: StereoGroup) {
        self.stereo_groups.push(group);
    }
}

// ---------------------------------------------------------------------------
// Connectivity utilities
// ---------------------------------------------------------------------------

impl Molecule {
    /// Return `true` if the molecule has exactly one connected component
    /// (i.e. every atom can be reached from every other atom).
    pub fn is_connected(&self) -> bool {
        let n = self.atoms.len();
        if n == 0 {
            return true;
        }
        let mut visited = vec![false; n];
        let mut stack = vec![AtomIdx(0)];
        visited[0] = true;
        let mut count = 1;
        while let Some(cur) = stack.pop() {
            for (nb, _) in self.neighbors(cur) {
                if !visited[nb.0 as usize] {
                    visited[nb.0 as usize] = true;
                    count += 1;
                    stack.push(nb);
                }
            }
        }
        count == n
    }

    /// Split the molecule into its connected components.
    ///
    /// Returns a `Vec` of sub-molecules, one per component.  Atoms are
    /// renumbered within each sub-molecule starting at index 0.
    pub fn fragments(&self) -> Vec<Molecule> {
        let n = self.atoms.len();
        if n == 0 {
            return vec![];
        }

        let mut component: Vec<usize> = vec![usize::MAX; n];
        let mut comp_id = 0;

        for start in 0..n {
            if component[start] != usize::MAX {
                continue;
            }
            let mut stack = vec![start];
            component[start] = comp_id;
            while let Some(cur) = stack.pop() {
                for (nb, _) in self.neighbors(AtomIdx(cur as u32)) {
                    let ni = nb.0 as usize;
                    if component[ni] == usize::MAX {
                        component[ni] = comp_id;
                        stack.push(ni);
                    }
                }
            }
            comp_id += 1;
        }

        (0..comp_id)
            .map(|cid| {
                let mut builder = MoleculeBuilder::new();
                let mut old_to_new: std::collections::HashMap<AtomIdx, AtomIdx> =
                    std::collections::HashMap::new();
                for (aidx, atom) in self.atoms() {
                    if component[aidx.0 as usize] == cid {
                        let new_idx = builder.add_atom(atom.clone());
                        old_to_new.insert(aidx, new_idx);
                    }
                }
                for (_, bond) in self.bonds() {
                    if let (Some(&a1), Some(&a2)) =
                        (old_to_new.get(&bond.atom1), old_to_new.get(&bond.atom2))
                    {
                        let _ = builder.add_bond(a1, a2, bond.order);
                    }
                }
                builder.build()
            })
            .collect()
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
    stereo_groups: Vec<StereoGroup>,
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a builder pre-populated with all atoms and bonds from `mol`.
    ///
    /// Use this to make incremental edits to an existing molecule instead of
    /// reconstructing it from scratch.
    pub fn from_molecule(mol: &Molecule) -> Self {
        let mut b = Self::new();
        for (_, atom) in mol.atoms() {
            b.add_atom(atom.clone());
        }
        for (_, bond) in mol.bonds() {
            let _ = b.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        b.stereo_groups = mol.stereo_groups.clone();
        b
    }

    /// Append a stereo group to this builder.
    pub fn add_stereo_group(&mut self, group: StereoGroup) {
        self.stereo_groups.push(group);
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
        self.adjacency[idx.0 as usize]
            .iter()
            .map(|&(nb, bidx)| (bidx, nb))
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
    pub fn add_bond(
        &mut self,
        a: AtomIdx,
        b: AtomIdx,
        order: BondOrder,
    ) -> Result<BondIdx, MolError> {
        let n = self.atoms.len() as u32;
        if a.0 >= n {
            return Err(MolError::InvalidAtomIdx(a));
        }
        if b.0 >= n {
            return Err(MolError::InvalidAtomIdx(b));
        }

        // Check for duplicate
        for &(nb, _) in &self.adjacency[a.0 as usize] {
            if nb == b {
                return Err(MolError::DuplicateBond(a, b));
            }
        }

        let bidx = BondIdx(self.bonds.len() as u32);
        self.bonds.push(BondEntry {
            atom1: a,
            atom2: b,
            order,
        });
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
            stereo_groups: self.stereo_groups,
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

    #[test]
    fn test_implicit_hydrogen_count() {
        // Isolated C atom (sp3, 4 bonds available): 4 implicit H
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::organic(Element::C));
        let mol = b.build();
        assert_eq!(mol.implicit_hydrogen_count(AtomIdx(0)), 4);
    }

    #[test]
    fn test_total_formula_methane() {
        // Organic C atom with 0 explicit bonds → 4 implicit H → CH4
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::organic(Element::C));
        let mol = b.build();
        assert_eq!(mol.total_formula(), "CH4");
    }

    #[test]
    fn test_total_formula_no_hydrogen() {
        // NaCl — neither Na nor Cl is in the organic subset, no implicit H
        let mut b = MoleculeBuilder::new();
        let na = b.add_atom(Atom::new(Element::NA));
        let cl = b.add_atom(Atom::new(Element::CL));
        b.add_bond(na, cl, BondOrder::Single).unwrap();
        let mol = b.build();
        assert_eq!(mol.total_formula(), "ClNa");
    }

    #[test]
    fn test_with_atom_aromatic() {
        let mol = ethane();
        let updated = mol.with_atom_aromatic(AtomIdx(0), true);
        assert!(updated.atom(AtomIdx(0)).aromatic);
        assert!(!updated.atom(AtomIdx(1)).aromatic);
    }

    #[test]
    fn test_with_bond_order() {
        let mol = ethane();
        let updated = mol.with_bond_order(BondIdx(0), BondOrder::Double);
        assert_eq!(updated.bond(BondIdx(0)).order, BondOrder::Double);
    }

    // --- mutable API ---

    #[test]
    fn test_add_remove_atom() {
        let mut mol = ethane();
        let n_idx = mol.add_atom(Atom::new(Element::N));
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.atom(n_idx).element.atomic_number(), 7);

        let remap = mol.remove_atom(n_idx);
        assert_eq!(mol.atom_count(), 2);
        assert!(remap[n_idx.0 as usize].is_none());
    }

    #[test]
    fn test_add_remove_bond() {
        let mut mol = ethane();
        let n_idx = mol.add_atom(Atom::new(Element::N));
        let bidx = mol.add_bond(AtomIdx(0), n_idx, BondOrder::Single).unwrap();
        assert_eq!(mol.bond_count(), 2);
        mol.remove_bond(bidx);
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn test_set_charge_element() {
        let mut mol = ethane();
        mol.set_charge(AtomIdx(0), 1);
        assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
        mol.set_element(AtomIdx(0), Element::N);
        assert_eq!(mol.atom(AtomIdx(0)).element.atomic_number(), 7);
    }

    #[test]
    fn test_is_connected() {
        let mol = ethane();
        assert!(mol.is_connected());

        // Two separate atoms — disconnected
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::new(Element::C));
        b.add_atom(Atom::new(Element::N));
        let disconnected = b.build();
        assert!(!disconnected.is_connected());
    }

    #[test]
    fn test_fragments() {
        // "CC.N" — two components
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::organic(Element::C));
        let c2 = b.add_atom(Atom::organic(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_atom(Atom::new(Element::N)); // disconnected N
        let mol = b.build();
        let frags = mol.fragments();
        assert_eq!(frags.len(), 2);
        let sizes: std::collections::HashSet<usize> =
            frags.iter().map(|f| f.atom_count()).collect();
        assert!(sizes.contains(&2));
        assert!(sizes.contains(&1));
    }

    #[test]
    fn test_builder_from_molecule() {
        let mol = ethane();
        let mut b = MoleculeBuilder::from_molecule(&mol);
        b.add_atom(Atom::new(Element::O));
        let mol2 = b.build();
        assert_eq!(mol2.atom_count(), 3);
        assert_eq!(mol2.bond_count(), 1); // original bond preserved
    }
}
