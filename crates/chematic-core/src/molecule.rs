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
/// Sentinel used in `stereo_neighbor_order` to represent the implicit H in a bracket atom.
pub const STEREO_H_SENTINEL: u32 = u32::MAX;

#[derive(Clone)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<BondEntry>,
    /// adjacency[atom_idx] = list of (neighbor_atom_idx, bond_idx)
    adjacency: Vec<Vec<(AtomIdx, BondIdx)>>,
    /// Enhanced stereo groups (ChemDraw V3000 Absolute / Or / And).
    stereo_groups: Vec<StereoGroup>,
    /// SMILES-text-order neighbor sequence for chiral atoms.
    ///
    /// Keyed by atom index.  Each value lists the atom indices of neighbors in
    /// the order they appeared in the SMILES string (including ring-closure
    /// partners), with [`STEREO_H_SENTINEL`] (`u32::MAX`) standing in for the
    /// implicit bracket H.  Populated by the SMILES parser; absent for atoms
    /// not parsed from SMILES or without recorded stereo.
    stereo_neighbor_order: std::collections::HashMap<u32, Vec<u32>>,
    /// Directional (`/`, `\`) bond marker, stashed for bonds whose `order` was
    /// overwritten to `Aromatic` (e.g. an exocyclic C=N adjacent to an
    /// aromatic ring atom — `order` must stay `Aromatic` for SMARTS `:a`
    /// matching, but the E/Z direction would otherwise be lost). Keyed by
    /// bond index; value is `BondOrder::Up` or `BondOrder::Down`. Absent for
    /// bonds whose direction is already carried directly by `order`.
    bond_directions: std::collections::HashMap<u32, BondOrder>,
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
    ///
    /// For a non-panicking variant, use [`Self::atom_opt`].
    pub fn atom(&self, idx: AtomIdx) -> &Atom {
        let i = idx.0 as usize;
        if i >= self.atoms.len() {
            panic!("atom index out of range");
        }
        &self.atoms[i]
    }

    /// Borrow atom by index, returning `None` if out of range.
    pub fn atom_opt(&self, idx: AtomIdx) -> Option<&Atom> {
        let i = idx.0 as usize;
        if i < self.atoms.len() {
            Some(&self.atoms[i])
        } else {
            None
        }
    }

    /// Borrow bond by index.
    ///
    /// # Panics
    /// Panics if `idx` is out of range (should not happen with indices from this molecule).
    ///
    /// For a non-panicking variant, use [`Self::bond_opt`].
    pub fn bond(&self, idx: BondIdx) -> &BondEntry {
        let i = idx.0 as usize;
        if i >= self.bonds.len() {
            panic!("bond index out of range");
        }
        &self.bonds[i]
    }

    /// Borrow bond by index, returning `None` if out of range.
    pub fn bond_opt(&self, idx: BondIdx) -> Option<&BondEntry> {
        let i = idx.0 as usize;
        if i < self.bonds.len() {
            Some(&self.bonds[i])
        } else {
            None
        }
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
    ///
    /// # Panics
    /// Panics if `idx` is out of range (should not happen with indices from this molecule).
    ///
    /// For a non-panicking variant, use [`Self::neighbors_opt`].
    pub fn neighbors(&self, idx: AtomIdx) -> impl Iterator<Item = (AtomIdx, BondIdx)> + '_ {
        let i = idx.0 as usize;
        if i >= self.adjacency.len() {
            // Keep the panic diagnostic free of input-derived values. Callers
            // handling untrusted atom indices should use `degree_opt`.
            panic!("atom index out of range");
        }
        self.adjacency[i].iter().copied()
    }

    /// Iterate over neighbors of `idx` as `(neighbor_atom_idx, bond_idx)`, returning `None` if out of range.
    pub fn neighbors_opt(&self, idx: AtomIdx) -> Option<Vec<(AtomIdx, BondIdx)>> {
        let i = idx.0 as usize;
        if i < self.adjacency.len() {
            Some(self.adjacency[i].to_vec())
        } else {
            None
        }
    }

    /// Degree (number of connected bonds) of atom `idx`.
    ///
    /// # Panics
    /// Panics if `idx` is out of range (should not happen with indices from this molecule).
    ///
    /// For a non-panicking variant, use [`Self::degree_opt`].
    pub fn degree(&self, idx: AtomIdx) -> usize {
        let i = idx.0 as usize;
        if i >= self.adjacency.len() {
            panic!("atom index out of range");
        }
        self.adjacency[i].len()
    }

    /// Degree (number of connected bonds) of atom `idx`, returning `None` if out of range.
    pub fn degree_opt(&self, idx: AtomIdx) -> Option<usize> {
        let i = idx.0 as usize;
        if i < self.adjacency.len() {
            Some(self.adjacency[i].len())
        } else {
            None
        }
    }

    /// Return the bond between `a` and `b`, or `None` if not connected or indices are out of bounds.
    pub fn bond_between(&self, a: AtomIdx, b: AtomIdx) -> Option<(BondIdx, &BondEntry)> {
        let a_idx = a.0 as usize;
        let b_idx = b.0 as usize;
        if a_idx >= self.adjacency.len() || b_idx >= self.atoms.len() {
            return None;
        }
        self.adjacency[a_idx]
            .iter()
            .find(|&&(nb, _)| nb == b)
            .and_then(|&(_, bidx)| {
                let bond_idx = bidx.0 as usize;
                if bond_idx < self.bonds.len() {
                    Some((bidx, &self.bonds[bond_idx]))
                } else {
                    None
                }
            })
    }

    /// Molecular formula as a Hill-order string (C first, H second, then alphabetical).
    pub fn formula(&self) -> String {
        use std::collections::BTreeMap;
        let mut counts: BTreeMap<&str, u32> = BTreeMap::new();
        for (_, atom) in self.atoms() {
            *counts.entry(atom.element.symbol()).or_insert(0) += 1;
        }
        let mut result = Self::format_hill_order_formula(&counts);
        let total_charge: i32 = self.atoms().map(|(_, a)| a.charge as i32).sum();
        match total_charge {
            0 => {}
            1 => result.push('+'),
            -1 => result.push('-'),
            n if n > 0 => result.push_str(&format!("+{n}")),
            n => result.push_str(&n.to_string()),
        }
        result
    }
}

// ---------------------------------------------------------------------------
// Immutable update methods (functional-style editing)
// ---------------------------------------------------------------------------

impl Molecule {
    /// Format element counts in Hill order: C, H, then alphabetically.
    fn format_hill_order_formula(counts: &std::collections::BTreeMap<&str, u32>) -> String {
        let mut counts = counts.clone();
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

    /// Return a new `Molecule` with one extra atom appended, along with the
    /// index that the new atom will have in the returned molecule.
    pub fn with_atom_added(&self, atom: Atom) -> (Molecule, AtomIdx) {
        let mut builder = MoleculeBuilder::from_molecule(self);
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
        let mut builder = MoleculeBuilder::from_molecule(self);
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
        builder.copy_stereo_from(self);
        builder.copy_bond_directions_from(self);
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
        builder.copy_stereo_from(self);
        builder.copy_bond_directions_from(self);
        // Chirality was cleared for the changed atom; remove its stereo order too.
        builder.clear_stereo_neighbor_order(idx);
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
        // Track old→new BOND index so `bond_directions` (keyed by bond index,
        // not atom index) can be remapped the same way `stereo_neighbor_order`
        // is remapped below — a prior version of this method dropped
        // `bond_directions` entirely on atom removal (silent loss, not
        // misattribution, but still a real gap the E/Z-direction side
        // channel needs closed: see `Molecule::bond_direction`).
        let mut bond_remap: Vec<Option<BondIdx>> = vec![None; self.bonds.len()];
        for (old_bidx, bond) in self.bonds() {
            if bond.atom1 == idx || bond.atom2 == idx {
                continue;
            }
            if let (Some(a1), Some(a2)) =
                (remap[bond.atom1.0 as usize], remap[bond.atom2.0 as usize])
                && let Ok(new_bidx) = builder.add_bond(a1, a2, bond.order)
            {
                bond_remap[old_bidx.0 as usize] = Some(new_bidx);
            }
        }
        // Remap stereo neighbor order: drop removed atom's entry, remap neighbor indices.
        for (old_key, order) in &self.stereo_neighbor_order {
            let old_atom = *old_key as usize;
            if old_atom == removed {
                continue; // removed atom's stereo is gone
            }
            if let Some(Some(new_key)) = remap.get(old_atom) {
                let new_order: Vec<u32> = order
                    .iter()
                    .filter_map(|&v| {
                        if v == STEREO_H_SENTINEL {
                            Some(STEREO_H_SENTINEL)
                        } else if v as usize == removed {
                            None // neighbor was the removed atom — stereo is now invalid
                        } else {
                            remap.get(v as usize).and_then(|r| r.map(|a| a.0))
                        }
                    })
                    .collect();
                builder.set_stereo_neighbor_order(*new_key, new_order);
            }
        }
        for (old_bidx, direction) in &self.bond_directions {
            if let Some(Some(new_bidx)) = bond_remap.get(*old_bidx as usize) {
                builder.set_bond_direction(*new_bidx, *direction);
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
        Self::format_hill_order_formula(&counts)
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
        builder.copy_stereo_from(self);
        builder.copy_bond_directions_from(self);
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
        builder.copy_stereo_from(self);
        builder.copy_bond_directions_from(self);
        builder.build()
    }

    /// Return a new `Molecule` with bond `idx` removed.
    ///
    /// Atom indices are unchanged.  Bond indices of survivors shift down, so
    /// `bond_directions` (keyed by bond index) is remapped bond-by-bond
    /// rather than copied wholesale — `copy_bond_directions_from` would
    /// misattribute directions to the wrong bond for every survivor after
    /// the removed one.
    pub fn with_bond_removed(&self, idx: BondIdx) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        for (_, atom) in self.atoms() {
            builder.add_atom(atom.clone());
        }
        for (bidx, bond) in self.bonds() {
            if bidx == idx {
                continue;
            }
            if let Ok(new_bidx) = builder.add_bond(bond.atom1, bond.atom2, bond.order)
                && let Some(direction) = self.bond_direction(bidx)
            {
                builder.set_bond_direction(new_bidx, direction);
            }
        }
        builder.copy_stereo_from(self);
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

        // Keep only bonds not involving the removed atom; remap endpoints and
        // track each surviving bond's new index so `bond_directions` (keyed
        // by bond index) can be remapped the same way as `stereo_neighbor_order`
        // is remapped by atom index below.
        let mut new_bonds: Vec<BondEntry> = Vec::new();
        let mut bond_remap: Vec<Option<u32>> = vec![None; self.bonds.len()];
        for (old_bidx, bond) in self.bonds.iter().enumerate() {
            if bond.atom1 == idx || bond.atom2 == idx {
                continue;
            }
            if let (Some(a1), Some(a2)) =
                (remap[bond.atom1.0 as usize], remap[bond.atom2.0 as usize])
            {
                bond_remap[old_bidx] = Some(new_bonds.len() as u32);
                new_bonds.push(BondEntry {
                    atom1: a1,
                    atom2: a2,
                    order: bond.order,
                });
            }
        }
        self.bonds = new_bonds;

        // Remap bond directions in-place, same shift as bond_remap above.
        let old_bond_directions = std::mem::take(&mut self.bond_directions);
        for (old_key, direction) in old_bond_directions {
            if let Some(Some(new_key)) = bond_remap.get(old_key as usize) {
                self.bond_directions.insert(*new_key, direction);
            }
        }

        // Rebuild adjacency from scratch.
        let new_n = self.atoms.len();
        self.adjacency = vec![vec![]; new_n];
        for (bidx, bond) in self.bonds.iter().enumerate() {
            let bi = BondIdx(bidx as u32);
            self.adjacency[bond.atom1.0 as usize].push((bond.atom2, bi));
            self.adjacency[bond.atom2.0 as usize].push((bond.atom1, bi));
        }

        // Remap stereo neighbor order in-place.
        let old_stereo = std::mem::take(&mut self.stereo_neighbor_order);
        for (old_key, order) in old_stereo {
            let old_atom = old_key as usize;
            if old_atom == removed {
                continue;
            }
            if let Some(Some(new_key)) = remap.get(old_atom) {
                let new_order: Vec<u32> = order
                    .iter()
                    .filter_map(|&v| {
                        if v == STEREO_H_SENTINEL {
                            Some(STEREO_H_SENTINEL)
                        } else if v as usize == removed {
                            None
                        } else {
                            remap.get(v as usize).and_then(|r| r.map(|a| a.0))
                        }
                    })
                    .collect();
                self.stereo_neighbor_order.insert(new_key.0, new_order);
            }
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
        // Remap `bond_directions` (keyed by bond index) for the same index
        // shift the bond list itself just underwent: entries below `removed`
        // are unchanged, the removed bond's own entry (if any) is dropped,
        // and entries above `removed` shift down by 1. A prior version of
        // this method left `bond_directions` untouched, which silently
        // MISATTRIBUTED a direction to whichever bond happened to shift into
        // the vacated slot -- worse than losing it, since it looks like valid
        // data pointing at the wrong physical bond.
        let old_bond_directions = std::mem::take(&mut self.bond_directions);
        for (old_key, direction) in old_bond_directions {
            let old = old_key as usize;
            match old.cmp(&removed) {
                std::cmp::Ordering::Less => {
                    self.bond_directions.insert(old_key, direction);
                }
                std::cmp::Ordering::Equal => {} // this bond itself was removed
                std::cmp::Ordering::Greater => {
                    self.bond_directions.insert(old_key - 1, direction);
                }
            }
        }
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

    /// Set the isotope label of atom `idx` in-place. `None` = natural
    /// isotope abundance (no label).
    pub fn set_isotope(&mut self, idx: AtomIdx, isotope: Option<u16>) {
        self.atoms[idx.0 as usize].isotope = isotope;
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

    /// Set the tetrahedral chirality (`@`/`@@`) of atom `idx` in-place.
    pub fn set_chirality(&mut self, idx: AtomIdx, chirality: crate::atom::Chirality) {
        self.atoms[idx.0 as usize].chirality = chirality;
    }

    /// Set the bond order of bond `idx` in-place. Endpoints (`atom1`/
    /// `atom2`) and adjacency are untouched -- order alone doesn't affect
    /// connectivity, so unlike [`Self::remove_bond`] + [`Self::add_bond`],
    /// this never perturbs any atom's `neighbors()` iteration order (some
    /// callers, e.g. 2D-wedge tetrahedral-parity perception, rely on that
    /// order to pick an "apex" neighbor -- a remove+re-add would silently
    /// change which neighbor that is).
    pub fn set_bond_order(&mut self, idx: BondIdx, order: BondOrder) {
        self.bonds[idx.0 as usize].order = order;
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

    /// SMILES-text-order neighbor sequence for a chiral atom.
    ///
    /// Returns `None` for atoms not parsed from SMILES or without stereo.
    /// The slice contains neighbor atom indices in SMILES text order;
    /// [`STEREO_H_SENTINEL`] (`u32::MAX`) marks the implicit bracket-H slot.
    ///
    /// # Invariant
    ///
    /// For any atom with `chirality != Chirality::None`, once this table is
    /// populated it must stay populated and correct relative to that atom's
    /// *current* neighbor set for as long as the chirality flag is set. A
    /// function that rebuilds a `Molecule` while keeping the same surviving
    /// atom/bond set (even if nothing actually changes) must carry this
    /// table forward explicitly (`MoleculeBuilder::copy_stereo_from`, plus
    /// `copy_bond_directions_from`/`copy_stereo_groups_from` for the other
    /// two stereo side tables) rather than leaving a downstream consumer to
    /// reconstruct it from raw adjacency — that reconstruction is only a
    /// best-effort fallback and is provably wrong for ring-opening
    /// stereocenters (see `chematic-chem`'s `hydrogen::declared_neighbor_order`
    /// and issue #399). A function that genuinely removes atoms or bonds
    /// must use the index-remap-with-sentinel-substitution pattern in
    /// [`Self::with_atom_removed`]/[`Self::with_bond_removed`], not a bare
    /// rebuild.
    pub fn stereo_neighbor_order(&self, idx: AtomIdx) -> Option<&[u32]> {
        self.stereo_neighbor_order.get(&idx.0).map(|v| v.as_slice())
    }

    /// Set the SMILES stereo neighbor order for atom `idx`.
    pub fn set_stereo_neighbor_order(&mut self, idx: AtomIdx, order: Vec<u32>) {
        self.stereo_neighbor_order.insert(idx.0, order);
    }

    /// Directional (`/`, `\`) marker stashed for bond `idx`, if its `order`
    /// was overwritten to `Aromatic` while it still carried E/Z direction.
    /// Returns `BondOrder::Up` or `BondOrder::Down` when present.
    pub fn bond_direction(&self, idx: BondIdx) -> Option<BondOrder> {
        self.bond_directions.get(&idx.0).copied()
    }

    /// Stash a directional marker for bond `idx` (see [`Self::bond_direction`]).
    pub fn set_bond_direction(&mut self, idx: BondIdx, direction: BondOrder) {
        self.bond_directions.insert(idx.0, direction);
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
    stereo_neighbor_order: std::collections::HashMap<u32, Vec<u32>>,
    bond_directions: std::collections::HashMap<u32, BondOrder>,
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
        b.stereo_neighbor_order = mol.stereo_neighbor_order.clone();
        b.bond_directions = mol.bond_directions.clone();
        b
    }

    /// Set the SMILES stereo neighbor order for atom `idx`.
    pub fn set_stereo_neighbor_order(&mut self, idx: AtomIdx, order: Vec<u32>) {
        self.stereo_neighbor_order.insert(idx.0, order);
    }

    /// Remove the stereo neighbor order entry for atom `idx`.
    pub fn clear_stereo_neighbor_order(&mut self, idx: AtomIdx) {
        self.stereo_neighbor_order.remove(&idx.0);
    }

    /// Append a stereo group to this builder.
    pub fn add_stereo_group(&mut self, group: StereoGroup) {
        self.stereo_groups.push(group);
    }

    /// Copy all enhanced stereo groups from `mol` into this builder verbatim.
    ///
    /// Only valid when atom indices are unchanged from `mol` (atoms re-added
    /// in the same order, none removed) — same caveat as
    /// [`Self::copy_bond_directions_from`].
    pub fn copy_stereo_groups_from(&mut self, mol: &Molecule) {
        self.stereo_groups = mol.stereo_groups.clone();
    }

    /// Copy all stereo neighbor order entries from `mol` into this builder.
    pub fn copy_stereo_from(&mut self, mol: &Molecule) {
        self.stereo_neighbor_order = mol.stereo_neighbor_order.clone();
    }

    /// Stash a directional marker for bond `idx` (see [`Molecule::bond_direction`]).
    pub fn set_bond_direction(&mut self, idx: BondIdx, direction: BondOrder) {
        self.bond_directions.insert(idx.0, direction);
    }

    /// Copy all bond-direction entries from `mol` into this builder verbatim.
    ///
    /// Only valid when bond indices are unchanged from `mol` (atoms/bonds
    /// re-added in the same order, none skipped) — e.g. a rebuild that only
    /// touches atom fields or promotes bond order to `Aromatic`. A rebuild
    /// that removes or reorders bonds must remap directions bond-by-bond
    /// instead (see `Molecule::with_bond_removed`).
    pub fn copy_bond_directions_from(&mut self, mol: &Molecule) {
        self.bond_directions = mol.bond_directions.clone();
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
            stereo_neighbor_order: self.stereo_neighbor_order,
            bond_directions: self.bond_directions,
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

    // --- bond_directions remap correctness (not just presence) ---

    /// 4-atom chain A-B-C-D with a `bond_direction` stash on the LAST bond
    /// (C-D, index 2). Removing the FIRST bond (A-B, index 0) shifts every
    /// surviving bond's index down by one; the stash must follow the C-D
    /// bond to its new index (1), not stay pinned to numeric index 2 (which
    /// would now point at a different physical bond) and not vanish.
    fn chain_with_direction_on_last_bond() -> (Molecule, BondIdx) {
        let mut b = MoleculeBuilder::new();
        let a = b.add_atom(Atom::new(Element::C));
        let bb = b.add_atom(Atom::new(Element::C));
        let c = b.add_atom(Atom::new(Element::C));
        let d = b.add_atom(Atom::new(Element::C));
        b.add_bond(a, bb, BondOrder::Single).unwrap(); // bond 0 (to be removed)
        b.add_bond(bb, c, BondOrder::Single).unwrap(); // bond 1
        let cd = b.add_bond(c, d, BondOrder::Single).unwrap(); // bond 2
        b.set_bond_direction(cd, BondOrder::Up);
        (b.build(), cd)
    }

    #[test]
    fn test_remove_bond_remaps_bond_direction_not_misattributes() {
        let (mut mol, _cd) = chain_with_direction_on_last_bond();
        assert_eq!(mol.bond_count(), 3);
        mol.remove_bond(BondIdx(0)); // remove A-B; C-D shifts from index 2 to 1
        assert_eq!(mol.bond_count(), 2);
        // The stash must have followed C-D to its new index...
        assert_eq!(mol.bond_direction(BondIdx(1)), Some(BondOrder::Up));
        // ...and must NOT have leaked onto the bond that shifted into the
        // old numeric slot 2 (which no longer exists) or onto B-C (index 0
        // after the shift), which never had a direction.
        assert_eq!(mol.bond_direction(BondIdx(0)), None);
        assert_eq!(mol.bond_opt(BondIdx(2)), None);
    }

    #[test]
    fn test_remove_bond_drops_direction_for_the_removed_bond_itself() {
        let (mut mol, _cd) = chain_with_direction_on_last_bond();
        mol.remove_bond(BondIdx(2)); // remove C-D itself — its stash must go with it
        assert_eq!(mol.bond_count(), 2);
        assert!(mol.bond_direction(BondIdx(0)).is_none());
        assert!(mol.bond_direction(BondIdx(1)).is_none());
    }

    #[test]
    fn test_with_atom_removed_remaps_bond_direction() {
        let (mol, _cd) = chain_with_direction_on_last_bond();
        // Remove atom A (index 0), unrelated to the C-D bond carrying the
        // stash. Bonds incident to A (A-B) disappear; B-C and C-D survive,
        // renumbered 0 and 1 respectively — the direction must follow C-D.
        let (updated, _atom_remap) = mol.with_atom_removed(AtomIdx(0));
        assert_eq!(updated.bond_count(), 2);
        // Find the surviving C-D bond by scanning for the stash directly,
        // rather than assuming a specific bond index, so this test doesn't
        // depend on internal re-numbering order.
        let has_direction = (0..updated.bond_count())
            .map(|i| BondIdx(i as u32))
            .any(|bidx| updated.bond_direction(bidx) == Some(BondOrder::Up));
        assert!(
            has_direction,
            "bond_direction on C-D must survive atom removal, remapped to its new bond index"
        );
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

    // --- safe Option-returning variants ---

    #[test]
    fn test_atom_opt_valid() {
        let mol = ethane();
        assert!(mol.atom_opt(AtomIdx(0)).is_some());
        assert!(mol.atom_opt(AtomIdx(1)).is_some());
        let atom = mol.atom_opt(AtomIdx(0)).unwrap();
        assert_eq!(atom.element.atomic_number(), 6);
    }

    #[test]
    fn test_atom_opt_invalid() {
        let mol = ethane();
        assert!(mol.atom_opt(AtomIdx(2)).is_none());
        assert!(mol.atom_opt(AtomIdx(1000)).is_none());
    }

    #[test]
    fn test_bond_opt_valid() {
        let mol = ethane();
        assert!(mol.bond_opt(BondIdx(0)).is_some());
        let bond = mol.bond_opt(BondIdx(0)).unwrap();
        assert_eq!(bond.order, BondOrder::Single);
    }

    #[test]
    fn test_bond_opt_invalid() {
        let mol = ethane();
        assert!(mol.bond_opt(BondIdx(1)).is_none());
        assert!(mol.bond_opt(BondIdx(1000)).is_none());
    }

    #[test]
    fn test_neighbors_opt_valid() {
        let mol = ethane();
        let neighbors = mol.neighbors_opt(AtomIdx(0)).unwrap();
        assert_eq!(neighbors.len(), 1);
        assert_eq!(neighbors[0].0, AtomIdx(1));
    }

    #[test]
    fn test_neighbors_opt_isolated_atom() {
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::new(Element::C));
        b.add_atom(Atom::new(Element::N));
        let mol = b.build();
        let neighbors = mol.neighbors_opt(AtomIdx(0)).unwrap();
        assert_eq!(neighbors.len(), 0);
    }

    #[test]
    fn test_neighbors_opt_invalid() {
        let mol = ethane();
        assert!(mol.neighbors_opt(AtomIdx(2)).is_none());
        assert!(mol.neighbors_opt(AtomIdx(1000)).is_none());
    }

    #[test]
    fn test_degree_opt_valid() {
        let mol = ethane();
        assert_eq!(mol.degree_opt(AtomIdx(0)), Some(1));
        assert_eq!(mol.degree_opt(AtomIdx(1)), Some(1));
    }

    #[test]
    fn test_degree_opt_isolated_atom() {
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::new(Element::C));
        b.add_atom(Atom::new(Element::N));
        let mol = b.build();
        assert_eq!(mol.degree_opt(AtomIdx(0)), Some(0));
        assert_eq!(mol.degree_opt(AtomIdx(1)), Some(0));
    }

    #[test]
    fn test_degree_opt_invalid() {
        let mol = ethane();
        assert!(mol.degree_opt(AtomIdx(2)).is_none());
        assert!(mol.degree_opt(AtomIdx(1000)).is_none());
    }

    #[test]
    fn test_degree_opt_multiple_bonds() {
        // Create a central atom with 3 neighbors
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(Element::C));
        let n1 = b.add_atom(Atom::new(Element::C));
        let n2 = b.add_atom(Atom::new(Element::N));
        let n3 = b.add_atom(Atom::new(Element::O));
        b.add_bond(center, n1, BondOrder::Single).unwrap();
        b.add_bond(center, n2, BondOrder::Double).unwrap();
        b.add_bond(center, n3, BondOrder::Single).unwrap();
        let mol = b.build();
        assert_eq!(mol.degree_opt(center), Some(3));
        assert_eq!(mol.degree_opt(n1), Some(1));
        assert_eq!(mol.degree_opt(n2), Some(1));
        assert_eq!(mol.degree_opt(n3), Some(1));
    }
}
