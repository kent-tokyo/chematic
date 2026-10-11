/// Builder for constructing a [`Molecule`] incrementally.
///
/// Usage: add atoms, add bonds, then call `build()`.
#[derive(Default)]
pub struct MoleculeBuilder {
    atoms: Vec<Atom>,
    atom_tags: Vec<Option<core::num::NonZeroU16>>,
    bonds: Vec<BondEntry>,
    adjacency: Vec<Vec<(AtomIdx, BondIdx)>>,
    r_groups: std::collections::HashMap<u32, RGroupLabel>,
    stereo_groups: Vec<StereoGroup>,
    stereo_neighbor_order: std::collections::HashMap<u32, Vec<u32>>,
    bond_directions: std::collections::HashMap<u32, BondOrder>,
    bond_direction_anchors: std::collections::HashMap<u32, AtomIdx>,
    smiles_ring_closure_keys: std::collections::HashMap<u32, u64>,
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    /// Create an empty builder with storage reserved for a known graph size.
    ///
    /// Parsers commonly learn atom and bond counts before constructing the
    /// graph. Reserving those primary arrays avoids repeated growth while
    /// preserving exactly the same builder semantics as [`Self::new`].
    pub fn with_capacity(atom_count: usize, bond_count: usize) -> Self {
        Self {
            atoms: Vec::with_capacity(atom_count),
            atom_tags: Vec::new(),
            bonds: Vec::with_capacity(bond_count),
            adjacency: Vec::with_capacity(atom_count),
            r_groups: std::collections::HashMap::new(),
            stereo_groups: Vec::new(),
            stereo_neighbor_order: std::collections::HashMap::new(),
            bond_directions: std::collections::HashMap::new(),
            bond_direction_anchors: std::collections::HashMap::new(),
            smiles_ring_closure_keys: std::collections::HashMap::new(),
        }
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
        b.r_groups = mol.r_groups.clone();
        b.atom_tags = mol.atom_tags.clone();
        b.stereo_neighbor_order = mol.stereo_neighbor_order.clone();
        b.bond_directions = mol.bond_directions.clone();
        b.bond_direction_anchors = mol.bond_direction_anchors.clone();
        b.smiles_ring_closure_keys = mol.smiles_ring_closure_keys.clone();
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

    /// Attach an R-group label to an already-added wildcard atom.
    pub fn set_r_group(&mut self, idx: AtomIdx, label: RGroupLabel) {
        let atom = &mut self.atoms[idx.0 as usize];
        atom.element = Element::C;
        atom.wildcard = true;
        atom.hydrogen_count = Some(0);
        atom.aromatic = false;
        self.r_groups.insert(idx.0, label);
    }

    /// Remove an R-group sidecar while leaving the atom itself in place.
    pub fn clear_r_group(&mut self, idx: AtomIdx) {
        self.r_groups.remove(&idx.0);
    }

    /// Copy R-group sidecars when atom indices are unchanged.
    pub fn copy_r_groups_from(&mut self, mol: &Molecule) {
        self.r_groups = mol.r_groups.clone();
    }

    /// Set a caller tag on an already-added atom. `None` and `Some(0)` clear it.
    /// Panics if `idx` is out of range.
    pub fn set_tag(&mut self, idx: AtomIdx, tag: Option<u16>) {
        assert!(
            (idx.0 as usize) < self.atom_count(),
            "atom index out of range"
        );
        let tag = tag.and_then(core::num::NonZeroU16::new);
        if tag.is_some() || !self.atom_tags.is_empty() {
            self.atom_tags.resize(self.atom_count(), None);
            self.atom_tags[idx.0 as usize] = tag;
        }
    }

    /// Copy caller tags when all source atoms retain their indices. Additional
    /// trailing atoms are untagged. Panics if this builder has fewer atoms.
    pub fn copy_atom_tags_from(&mut self, mol: &Molecule) {
        assert!(self.atom_count() >= mol.atom_count(), "atom count mismatch");
        self.atom_tags.clone_from(&mol.atom_tags);
        if !self.atom_tags.is_empty() {
            self.atom_tags.resize(self.atom_count(), None);
        }
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

    /// Set the parser-side endpoint that anchored a stashed direction.
    pub fn set_bond_direction_anchor(&mut self, idx: BondIdx, atom: AtomIdx) {
        self.bond_direction_anchors.insert(idx.0, atom);
    }

    /// Record that bond `idx` closed SMILES ring label `label` at its
    /// `occurrence`-th use (see [`Molecule::rdkit_bond_order`]).
    pub fn set_smiles_ring_closure(&mut self, idx: BondIdx, label: u32, occurrence: u32) {
        // Low two bits: whether the begin atom RDKit gives the bond is known
        // (bit 1) and is the ring-opening atom (bit 0); see
        // `set_smiles_ring_closure_begins_at_open`.
        self.smiles_ring_closure_keys.insert(
            idx.0,
            (u64::from(label) << 32) | ((u64::from(occurrence) << 2) & 0xFFFF_FFFF),
        );
    }

    /// Record, for ring-closure bond `idx` (already registered with
    /// [`MoleculeBuilder::set_smiles_ring_closure`]), whether RDKit's SMILES
    /// parser makes the ring-opening atom its begin atom (`CloseMolRings`
    /// keeps the opening partial bond only when it carried an explicit,
    /// non-directional bond symbol).
    pub fn set_smiles_ring_closure_begins_at_open(&mut self, idx: BondIdx, at_open: bool) {
        if let Some(k) = self.smiles_ring_closure_keys.get_mut(&idx.0) {
            *k = (*k & !3) | 2 | u64::from(at_open);
        }
    }

    /// Copy all bond-direction entries from `mol` into this builder verbatim.
    ///
    /// Only valid when bond indices are unchanged from `mol` (atoms/bonds
    /// re-added in the same order, none skipped) — e.g. a rebuild that only
    /// touches atom fields or promotes bond order to `Aromatic`. A rebuild
    /// that removes or reorders bonds must remap directions bond-by-bond
    /// instead (see `Molecule::with_bond_removed`).
    pub fn copy_bond_directions_from(&mut self, mol: &Molecule) {
        self.smiles_ring_closure_keys = mol.smiles_ring_closure_keys.clone();
        self.bond_directions = mol.bond_directions.clone();
        self.bond_direction_anchors = mol.bond_direction_anchors.clone();
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
        if !self.atom_tags.is_empty() {
            self.atom_tags.push(None);
        }
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
            atom_tags: self.atom_tags,
            bonds: self.bonds,
            adjacency: self.adjacency,
            r_groups: self.r_groups,
            stereo_groups: self.stereo_groups,
            stereo_neighbor_order: self.stereo_neighbor_order,
            bond_directions: self.bond_directions,
            bond_direction_anchors: self.bond_direction_anchors,
            smiles_ring_closure_keys: self.smiles_ring_closure_keys,
            derived: Default::default(),
        }
    }
}
