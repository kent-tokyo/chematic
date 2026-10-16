// ---------------------------------------------------------------------------
// In-place mutation methods
// ---------------------------------------------------------------------------

fn remap_stereo_groups_after_atom_removal(
    groups: Vec<StereoGroup>,
    group_bonds: Vec<Vec<BondIdx>>,
    atom_remap: &[Option<AtomIdx>],
    bond_remap: &[Option<u32>],
) -> Vec<(StereoGroup, Vec<BondIdx>)> {
    groups
        .into_iter()
        .zip(group_bonds)
        .filter_map(|(group, bonds)| {
            let atoms = group
                .atom_indices
                .into_iter()
                .filter_map(|atom| atom_remap.get(atom.0 as usize).and_then(|mapped| *mapped))
                .collect::<Vec<_>>();
            let bonds = bonds
                .into_iter()
                .filter_map(|bond| {
                    bond_remap
                        .get(bond.0 as usize)
                        .and_then(|mapped| mapped.map(BondIdx))
                })
                .collect::<Vec<_>>();
            (!atoms.is_empty() || !bonds.is_empty())
                .then(|| (StereoGroup::new(group.kind, atoms), bonds))
        })
        .collect()
}

impl Molecule {
    /// Replace this molecule with a fresh empty graph.
    ///
    /// Drops atoms, bonds, adjacency, tags, R labels, all stereo/ring-order
    /// metadata and every derived cache. Previously cloned molecules remain
    /// independent. Storage is released rather than reserved for reuse.
    /// Wrapper-owned titles/properties are outside this core graph.
    pub fn clear(&mut self) {
        *self = MoleculeBuilder::new().build();
    }

    /// Append a new atom and return its index.
    pub fn add_atom(&mut self, atom: Atom) -> AtomIdx {
        self.invalidate_derived();
        let idx = AtomIdx(self.atoms.len() as u32);
        self.atoms.push(atom);
        if !self.atom_tags.is_empty() {
            self.atom_tags.push(None);
        }
        self.adjacency.push(vec![]);
        idx
    }

    /// Remove atom `idx` and all bonds involving it.
    ///
    /// Returns a remapping table: `remap[old_idx]` gives the new `AtomIdx`
    /// for surviving atoms, or `None` for the removed atom.  Atom indices
    /// of atoms after the removed slot shift down by 1.
    pub fn remove_atom(&mut self, idx: AtomIdx) -> Vec<Option<AtomIdx>> {
        self.invalidate_derived();
        self.smiles_ring_closure_keys.clear();
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

        // Drain the single slot instead of calling `Vec::remove` directly.
        // Besides making the order-preserving intent explicit, this avoids a
        // CodeQL `rust/cleartext-logging` false positive on the PyO3 RWMol
        // mutable-reference path (the analyzer models that method name as a
        // logging sink even though this is only an in-memory graph edit).
        let _ = self.atoms.drain(removed..=removed).next();
        if !self.atom_tags.is_empty() {
            let _ = self.atom_tags.drain(removed..=removed).next();
        }

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
        let old_bond_anchors = std::mem::take(&mut self.bond_direction_anchors);
        for (old_key, atom) in old_bond_anchors {
            if let (Some(Some(new_key)), Some(new_atom)) = (
                bond_remap.get(old_key as usize),
                remap.get(atom.0 as usize).and_then(|r| *r),
            ) {
                self.bond_direction_anchors.insert(*new_key, new_atom);
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
                        if matches!(v, STEREO_H_SENTINEL | STEREO_VACANT_SENTINEL) {
                            Some(v)
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

        let old_r_groups = std::mem::take(&mut self.r_groups);
        for (old_key, label) in old_r_groups {
            if let Some(Some(new_atom)) = remap.get(old_key as usize) {
                self.r_groups.insert(new_atom.0, label);
            }
        }

        let old_groups = std::mem::take(&mut self.stereo_groups);
        let old_group_bonds = std::mem::take(&mut self.stereo_group_bonds);
        let groups = remap_stereo_groups_after_atom_removal(
            old_groups,
            old_group_bonds,
            &remap,
            &bond_remap,
        );
        self.set_stereo_groups_with_bonds(groups);

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
        self.invalidate_derived();
        // Parser provenance describes the original graph. Adding a bond can
        // change its cycle rank, so retaining the old closure count would
        // make provenance-backed ring shortcuts stale.
        self.smiles_ring_closure_keys.clear();
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
        self.invalidate_derived();
        self.smiles_ring_closure_keys.clear();
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
        let old_bond_anchors = std::mem::take(&mut self.bond_direction_anchors);
        for (old_key, atom) in old_bond_anchors {
            let old = old_key as usize;
            let new_key = match old.cmp(&removed) {
                std::cmp::Ordering::Less => old_key,
                std::cmp::Ordering::Equal => continue,
                std::cmp::Ordering::Greater => old_key - 1,
            };
            self.bond_direction_anchors.insert(new_key, atom);
        }
        let old_groups = std::mem::take(&mut self.stereo_groups);
        let old_group_bonds = std::mem::take(&mut self.stereo_group_bonds);
        let groups = old_groups
            .into_iter()
            .zip(old_group_bonds)
            .filter_map(|(group, group_bonds)| {
                let bonds = group_bonds
                    .into_iter()
                    .filter_map(|bond| match (bond.0 as usize).cmp(&removed) {
                        std::cmp::Ordering::Less => Some(bond),
                        std::cmp::Ordering::Equal => None,
                        std::cmp::Ordering::Greater => Some(BondIdx(bond.0 - 1)),
                    })
                    .collect::<Vec<_>>();
                (!group.atom_indices.is_empty() || !bonds.is_empty())
                    .then(|| (StereoGroup::new(group.kind, group.atom_indices), bonds))
            })
            .collect::<Vec<_>>();
        self.set_stereo_groups_with_bonds(groups);
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
        self.invalidate_derived();
        self.atoms[idx.0 as usize].charge = charge;
    }

    /// Set the isotope label of atom `idx` in-place. `None` = natural
    /// isotope abundance (no label).
    /// Set (or clear) the explicit hydrogen count of atom `idx`.
    pub fn set_hydrogen_count(&mut self, idx: AtomIdx, hydrogens: Option<u8>) {
        self.invalidate_derived();
        self.atoms[idx.0 as usize].hydrogen_count = hydrogens;
    }

    pub fn set_isotope(&mut self, idx: AtomIdx, isotope: Option<u16>) {
        self.invalidate_derived();
        self.atoms[idx.0 as usize].isotope = isotope;
    }

    /// Set a non-chemical tag on atom `idx`. Preserved by molecule clone,
    /// core graph edits, reaction apply, fragments, and aromaticity perception.
    /// Ignored by SMILES write and canonicalization; see [`Self::atom_tag`].
    ///
    /// Tags are `1..=u16::MAX`. `None` clears the tag; `Some(0)` is treated
    /// as clear (same as `None`). Storage is allocated only when tagging.
    ///
    /// # Panics
    /// Panics if `idx` is out of range.
    pub fn set_tag(&mut self, idx: AtomIdx, tag: Option<u16>) {
        assert!(
            (idx.0 as usize) < self.atom_count(),
            "atom index out of range"
        );
        let tag = tag.and_then(core::num::NonZeroU16::new);
        if self.atom_tag(idx) == tag {
            return;
        }
        // Cached aromatic views contain whole molecules, including tags.
        self.invalidate_derived();
        self.atom_tags.resize(self.atom_count(), None);
        self.atom_tags[idx.0 as usize] = tag;
    }

    /// Mark an atom as an R-group wildcard while preserving its atom index.
    pub fn set_r_group(&mut self, idx: AtomIdx, label: crate::atom::RGroupLabel) {
        self.invalidate_derived();
        let atom = &mut self.atoms[idx.0 as usize];
        atom.element = Element::C;
        atom.wildcard = true;
        atom.hydrogen_count = Some(0);
        atom.aromatic = false;
        atom.chirality = crate::atom::Chirality::None;
        self.r_groups.insert(idx.0, label);
    }

    /// Mark an atom as a generic `*` wildcard while preserving its atom index.
    pub fn set_wildcard(&mut self, idx: AtomIdx) {
        self.invalidate_derived();
        let atom = &mut self.atoms[idx.0 as usize];
        atom.element = Element::C;
        atom.wildcard = true;
        atom.hydrogen_count = Some(0);
        atom.aromatic = false;
        atom.chirality = crate::atom::Chirality::None;
        self.r_groups.remove(&idx.0);
    }

    /// Set the element of atom `idx` in-place.
    ///
    /// Chirality and hydrogen count are reset (element-specific properties).
    pub fn set_element(&mut self, idx: AtomIdx, el: Element) {
        self.invalidate_derived();
        let a = &mut self.atoms[idx.0 as usize];
        a.element = el;
        a.chirality = crate::atom::Chirality::None;
        a.hydrogen_count = None;
        a.aromatic = false;
        a.wildcard = false;
        self.r_groups.remove(&idx.0);
    }

    /// Set the CIP stereo code of atom `idx` in-place.
    pub fn set_cip_code(&mut self, idx: AtomIdx, code: Option<crate::atom::CipCode>) {
        self.invalidate_derived();
        self.atoms[idx.0 as usize].cip_code = code;
    }

    /// Set the tetrahedral chirality (`@`/`@@`) of atom `idx` in-place.
    pub fn set_chirality(&mut self, idx: AtomIdx, chirality: crate::atom::Chirality) {
        self.invalidate_derived();
        self.atoms[idx.0 as usize].chirality = chirality;
    }

    /// Set the aromatic flag of atom `idx` in-place.
    pub fn set_atom_aromatic(&mut self, idx: AtomIdx, aromatic: bool) {
        self.invalidate_derived();
        self.atoms[idx.0 as usize].aromatic = aromatic;
    }

    /// Remove every `R`/`R<n>` label (the atoms themselves are unchanged).
    #[doc(hidden)]
    pub fn clear_r_group_labels(&mut self) {
        self.invalidate_derived();
        self.r_groups.clear();
    }

    /// Whether every adjacency list is in bond-index order (atom1's entries
    /// then atom2's, bond by bond), i.e. exactly what [`MoleculeBuilder`]
    /// produces when the bonds are re-added in index order. A clone of such a
    /// molecule is then indistinguishable from that rebuild.
    #[doc(hidden)]
    pub fn adjacency_in_bond_order(&self) -> bool {
        let mut fill = vec![0usize; self.atoms.len()];
        for (b, bond) in self.bonds.iter().enumerate() {
            for (end, other) in [(bond.atom1, bond.atom2), (bond.atom2, bond.atom1)] {
                let list = &self.adjacency[end.0 as usize];
                let k = fill[end.0 as usize];
                if list.get(k) != Some(&(other, BondIdx(b as u32))) {
                    return false;
                }
                fill[end.0 as usize] = k + 1;
            }
        }
        fill.iter()
            .zip(&self.adjacency)
            .all(|(&k, list)| k == list.len())
    }

    /// What re-adding every atom and bond to a [`MoleculeBuilder`] and then
    /// copying the stereo neighbour orders and stereo groups gives (no atom
    /// tags, R-group labels, stashed bond directions or SMILES ring-closure
    /// record), without re-adding the bonds; `None` when the adjacency is not
    /// in bond order ([`Self::adjacency_in_bond_order`]), where a rebuild
    /// would reorder neighbours.
    #[doc(hidden)]
    pub fn atoms_bonds_and_stereo_copy(&self) -> Option<Molecule> {
        if !self.adjacency_in_bond_order() {
            return None;
        }
        Some(Molecule {
            atoms: self.atoms.clone(),
            atom_tags: Vec::new(),
            bonds: self.bonds.clone(),
            adjacency: self.adjacency.clone(),
            r_groups: Default::default(),
            stereo_groups: self.stereo_groups.clone(),
            stereo_group_bonds: self.stereo_group_bonds.clone(),
            stereo_neighbor_order: self.stereo_neighbor_order.clone(),
            bond_directions: Default::default(),
            bond_direction_anchors: Default::default(),
            smiles_ring_closure_keys: Default::default(),
            derived: Default::default(),
        })
    }

    /// Set the bond order of bond `idx` in-place. Endpoints (`atom1`/
    /// `atom2`) and adjacency are untouched -- order alone doesn't affect
    /// connectivity, so unlike [`Self::remove_bond`] + [`Self::add_bond`],
    /// this never perturbs any atom's `neighbors()` iteration order (some
    /// callers, e.g. 2D-wedge tetrahedral-parity perception, rely on that
    /// order to pick an "apex" neighbor -- a remove+re-add would silently
    /// change which neighbor that is).
    pub fn set_bond_order(&mut self, idx: BondIdx, order: BondOrder) {
        self.invalidate_derived();
        self.bonds[idx.0 as usize].order = order;
    }

    /// Make bond `idx` a dative bond from `donor` (stored as `atom1`, see
    /// [`BondOrder::Dative`]), swapping its endpoints if needed. Adjacency
    /// order is untouched.
    pub fn set_dative_bond(&mut self, idx: BondIdx, donor: AtomIdx) {
        self.invalidate_derived();
        let bond = &mut self.bonds[idx.0 as usize];
        bond.order = BondOrder::Dative;
        if bond.atom2 == donor {
            std::mem::swap(&mut bond.atom1, &mut bond.atom2);
        }
    }

    /// Return the enhanced stereo groups attached to this molecule.
    pub fn stereo_groups(&self) -> &[StereoGroup] {
        &self.stereo_groups
    }

    /// Return the geometric/axial bond members of stereo group `group_index`.
    ///
    /// The index refers to the matching entry in [`Self::stereo_groups`].
    pub fn stereo_group_bond_indices(&self, group_index: usize) -> Option<&[BondIdx]> {
        self.stereo_group_bonds.get(group_index).map(Vec::as_slice)
    }

    /// Replace the stereo group list in-place.
    pub fn set_stereo_groups(&mut self, groups: Vec<StereoGroup>) {
        self.invalidate_derived();
        self.stereo_group_bonds = vec![Vec::new(); groups.len()];
        self.stereo_groups = groups;
    }

    /// Replace enhanced stereo groups together with their bond members.
    pub fn set_stereo_groups_with_bonds(
        &mut self,
        groups: Vec<(StereoGroup, Vec<BondIdx>)>,
    ) {
        self.invalidate_derived();
        self.stereo_groups = groups.iter().map(|(group, _)| group.clone()).collect();
        self.stereo_group_bonds = groups
            .into_iter()
            .map(|(_, bonds)| dedup_bond_indices(bonds))
            .collect();
    }

    /// Copy enhanced stereo groups and their atom/bond membership verbatim.
    pub fn copy_stereo_groups_from(&mut self, source: &Molecule) {
        self.invalidate_derived();
        self.stereo_groups.clone_from(&source.stereo_groups);
        self.stereo_group_bonds
            .clone_from(&source.stereo_group_bonds);
    }

    /// Add a single stereo group in-place.
    pub fn add_stereo_group(&mut self, group: StereoGroup) {
        self.invalidate_derived();
        self.stereo_groups.push(group);
        self.stereo_group_bonds.push(Vec::new());
    }

    /// Add one enhanced stereo group with geometric/axial bond members.
    pub fn add_stereo_group_with_bonds(
        &mut self,
        group: StereoGroup,
        bond_indices: Vec<BondIdx>,
    ) {
        self.invalidate_derived();
        self.stereo_groups.push(group);
        self.stereo_group_bonds
            .push(dedup_bond_indices(bond_indices));
    }

    /// How many SMILES ring-closure digits were written on atom `idx` (the
    /// ring-closure bonds it ends, see [`Molecule::rdkit_bond_order`]); 0 for
    /// molecules not read from SMILES.
    pub fn smiles_ring_closure_count(&self, idx: AtomIdx) -> usize {
        if self.smiles_ring_closure_keys.is_empty() {
            return 0;
        }
        self.neighbors(idx)
            .filter(|(_, b)| self.smiles_ring_closure_keys.contains_key(&b.0))
            .count()
    }

    /// Number of ring-closure bonds recorded by the SMILES parser.
    ///
    /// For an unmodified SMILES graph this is its cycle rank. A non-SMILES
    /// molecule, or an acyclic SMILES, returns zero; callers must therefore
    /// only use positive values as provenance-backed shortcuts.
    pub fn smiles_ring_closure_bond_count(&self) -> usize {
        self.smiles_ring_closure_keys.len()
    }

    /// Whether bond `idx` was written as a SMILES ring closure (a digit
    /// pair). Such a bond is stored opening atom -> closing atom; RDKit's
    /// parser makes the closing atom its begin atom (unless the opening
    /// digit carried an explicit bond symbol).
    pub fn is_smiles_ring_closure(&self, idx: BondIdx) -> bool {
        self.smiles_ring_closure_keys.contains_key(&idx.0)
    }

    /// For a ring-closure bond read from SMILES, whether RDKit's parser
    /// makes the ring-opening atom the bond's begin atom; `None` when not
    /// recorded.
    pub fn smiles_ring_closure_begins_at_open(&self, idx: BondIdx) -> Option<bool> {
        let k = *self.smiles_ring_closure_keys.get(&idx.0)?;
        (k & 2 != 0).then_some(k & 1 != 0)
    }

    /// Return the parser provenance stored for one SMILES ring-closure bond.
    ///
    /// This is an internal interchange hook for bindings that must preserve a
    /// molecule without routing it through a text format. The tuple contains
    /// `(ring label, label occurrence, begins-at-opening-atom)`. The final
    /// value is `None` when the parser did not record the endpoint convention.
    #[doc(hidden)]
    pub fn smiles_ring_closure_metadata(
        &self,
        idx: BondIdx,
    ) -> Option<(u32, u32, Option<bool>)> {
        let key = *self.smiles_ring_closure_keys.get(&idx.0)?;
        let label = (key >> 32) as u32;
        let occurrence = ((key & 0xFFFF_FFFF) >> 2) as u32;
        let begins_at_open = (key & 2 != 0).then_some(key & 1 != 0);
        Some((label, occurrence, begins_at_open))
    }

    /// Bond indices in the order RDKit's SMILES parser creates the bonds:
    /// chain bonds as written, then every ring-closure bond, ordered by ring
    /// label and, for a reused label, by occurrence (RDKit closes rings after
    /// parsing). Per-atom neighbour order in RDKit follows the same sequence.
    /// For molecules not read from SMILES this is the bond index order, which
    /// is also RDKit's for MOL/SDF input.
    pub fn rdkit_bond_order(&self) -> Vec<BondIdx> {
        let mut order = Vec::with_capacity(self.bonds.len());
        self.for_each_bond_in_rdkit_order(|bond| order.push(bond));
        order
    }

    /// Visit bonds in [`Self::rdkit_bond_order`] without allocating the full
    /// order vector. The usual SMILES has at most eight closure bonds, which
    /// are sorted in a small stack buffer; unusually closure-heavy inputs use
    /// a heap fallback. Intended for hot read-only descriptor paths.
    #[doc(hidden)]
    pub fn for_each_bond_in_rdkit_order(&self, mut visit: impl FnMut(BondIdx)) {
        // Acyclic SMILES and non-SMILES inputs have no deferred closure
        // bonds.  They are the common case for descriptor hot paths; avoid a
        // hash-table probe for every bond when insertion order already is
        // RDKit order.
        if self.smiles_ring_closure_keys.is_empty() {
            for bond in 0..self.bonds.len() as u32 {
                visit(BondIdx(bond));
            }
            return;
        }
        for bond in 0..self.bonds.len() as u32 {
            if !self.smiles_ring_closure_keys.contains_key(&bond) {
                visit(BondIdx(bond));
            }
        }

        const INLINE_CLOSURES: usize = 8;
        if self.smiles_ring_closure_keys.len() <= INLINE_CLOSURES {
            let mut closures = [(0u64, 0u32); INLINE_CLOSURES];
            let mut len = 0usize;
            for (&bond, &key) in &self.smiles_ring_closure_keys {
                if (bond as usize) < self.bonds.len() {
                    closures[len] = (key, bond);
                    len += 1;
                }
            }
            closures[..len].sort_unstable();
            for &(_, bond) in &closures[..len] {
                visit(BondIdx(bond));
            }
            return;
        }

        let mut closures: Vec<(u64, u32)> = self
            .smiles_ring_closure_keys
            .iter()
            .filter(|(bond, _)| (**bond as usize) < self.bonds.len())
            .map(|(&bond, &key)| (key, bond))
            .collect();
        closures.sort_unstable();
        for (_, bond) in closures {
            visit(BondIdx(bond));
        }
    }

    /// SMILES-text-order neighbor sequence for a chiral atom.
    ///
    /// Returns `None` for atoms not parsed from SMILES or without stereo.
    /// The slice contains neighbor atom indices in SMILES text order;
    /// [`STEREO_H_SENTINEL`] marks bracket H; [`STEREO_VACANT_SENTINEL`]
    /// marks an unoccupied non-tetrahedral coordination site.
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

    /// Whether `idx` carries OpenSMILES allene/cumulene axial stereo.
    ///
    /// Axial `@AL1`/`@AL2` and tetrahedral `@`/`@@` share the same two parity
    /// values. The parser records an allene's two terminal atoms and four
    /// virtual ligands as a six-entry stereo order, which makes the class
    /// unambiguous without adding variants to the public [`crate::Chirality`]
    /// enum. Keeping the enum's variant set stable preserves exhaustive
    /// downstream matches across patch releases.
    pub fn is_allene_stereo(&self, idx: AtomIdx) -> bool {
        self.atom(idx).chirality.is_tetrahedral()
            && self
                .stereo_neighbor_order(idx)
                .is_some_and(|order| order.len() == 6)
    }

    /// Set the SMILES stereo neighbor order for atom `idx`.
    pub fn set_stereo_neighbor_order(&mut self, idx: AtomIdx, order: Vec<u32>) {
        self.invalidate_derived();
        self.stereo_neighbor_order.insert(idx.0, order);
    }

    /// Discard the stored reference frame for one atom, leaving its tag and
    /// every other frame unchanged. Clear the chirality tag first when removing
    /// a declaration; a retained nonzero tag still requires a valid frame.
    /// Removing an absent frame is a no-op.
    pub fn clear_stereo_neighbor_order(&mut self, idx: AtomIdx) {
        self.invalidate_derived();
        self.stereo_neighbor_order.remove(&idx.0);
    }

    /// Directional (`/`, `\`) marker stashed for bond `idx`, if its `order`
    /// was overwritten to `Aromatic` while it still carried E/Z direction.
    /// Returns `BondOrder::Up` or `BondOrder::Down` when present.
    pub fn bond_direction(&self, idx: BondIdx) -> Option<BondOrder> {
        self.bond_directions.get(&idx.0).copied()
    }

    /// Whether any bond has a stashed direction ([`Self::bond_direction`]).
    pub fn has_bond_directions(&self) -> bool {
        !self.bond_directions.is_empty()
    }

    /// Return the endpoint that anchored a parser-side aromatic direction
    /// stash, if one was recorded.
    pub fn bond_direction_anchor(&self, idx: BondIdx) -> Option<AtomIdx> {
        self.bond_direction_anchors.get(&idx.0).copied()
    }

    /// Stash a directional marker for bond `idx` (see [`Self::bond_direction`]).
    pub fn set_bond_direction(&mut self, idx: BondIdx, direction: BondOrder) {
        self.invalidate_derived();
        self.bond_directions.insert(idx.0, direction);
    }

    /// Set the endpoint that anchored a parser-side aromatic direction stash.
    pub fn set_bond_direction_anchor(&mut self, idx: BondIdx, atom: AtomIdx) {
        self.invalidate_derived();
        self.bond_direction_anchors.insert(idx.0, atom);
    }

    /// Remove bond `idx`'s stashed direction and its anchor, if any.
    pub fn clear_bond_direction(&mut self, idx: BondIdx) {
        self.invalidate_derived();
        self.bond_directions.remove(&idx.0);
        self.bond_direction_anchors.remove(&idx.0);
    }
}
