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
        builder.copy_atom_tags_from(self);
        builder.copy_r_groups_from(self);
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
                a.wildcard = false;
            }
            builder.add_atom(a);
        }
        for (_, bond) in self.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        builder.copy_stereo_from(self);
        builder.copy_atom_tags_from(self);
        builder.copy_r_groups_from(self);
        builder.clear_r_group(idx);
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
            let new_idx = builder.add_atom(atom.clone());
            builder.set_tag(new_idx, self.atom_tag(aidx).map(core::num::NonZeroU16::get));
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
                        if matches!(v, STEREO_H_SENTINEL | STEREO_VACANT_SENTINEL) {
                            Some(v)
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
        // Anchors follow their atoms (a stash is read from its anchor).
        for (old_bidx, anchor) in &self.bond_direction_anchors {
            if let (Some(Some(new_bidx)), Some(Some(new_anchor))) = (
                bond_remap.get(*old_bidx as usize),
                remap.get(anchor.0 as usize),
            ) {
                builder.set_bond_direction_anchor(*new_bidx, *new_anchor);
            }
        }
        for (old_atom, label) in &self.r_groups {
            if let Some(Some(new_atom)) = remap.get(*old_atom as usize) {
                builder.set_r_group(*new_atom, *label);
            }
        }
        for (group_index, group) in self.stereo_groups.iter().enumerate() {
            let atoms = group
                .atom_indices
                .iter()
                .filter_map(|atom| remap.get(atom.0 as usize).and_then(|mapped| *mapped))
                .collect::<Vec<_>>();
            let bonds = self.stereo_group_bonds[group_index]
                .iter()
                .filter_map(|bond| bond_remap.get(bond.0 as usize).and_then(|mapped| *mapped))
                .collect::<Vec<_>>();
            if !atoms.is_empty() || !bonds.is_empty() {
                builder.add_stereo_group_with_bonds(
                    StereoGroup::new(group.kind.clone(), atoms),
                    bonds,
                );
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
        builder.copy_atom_tags_from(self);
        builder.copy_r_groups_from(self);
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
        builder.copy_atom_tags_from(self);
        builder.copy_r_groups_from(self);
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
        let mut bond_remap = vec![None; self.bond_count()];
        for (_, atom) in self.atoms() {
            builder.add_atom(atom.clone());
        }
        for (bidx, bond) in self.bonds() {
            if bidx == idx {
                continue;
            }
            if let Ok(new_bidx) = builder.add_bond(bond.atom1, bond.atom2, bond.order) {
                bond_remap[bidx.0 as usize] = Some(new_bidx);
                if let Some(direction) = self.bond_direction(bidx) {
                    builder.set_bond_direction(new_bidx, direction);
                    if let Some(anchor) = self.bond_direction_anchor(bidx) {
                        builder.set_bond_direction_anchor(new_bidx, anchor);
                    }
                }
            }
        }
        builder.copy_stereo_from(self);
        builder.copy_atom_tags_from(self);
        builder.copy_r_groups_from(self);
        for (group_index, group) in self.stereo_groups.iter().enumerate() {
            let bonds = self.stereo_group_bonds[group_index]
                .iter()
                .filter_map(|bond| bond_remap.get(bond.0 as usize).and_then(|mapped| *mapped))
                .collect::<Vec<_>>();
            if !group.atom_indices.is_empty() || !bonds.is_empty() {
                builder.add_stereo_group_with_bonds(group.clone(), bonds);
            }
        }
        builder.build()
    }
}
