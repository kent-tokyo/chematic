impl Molecule {
    /// Return the memoized RDKit-compatible ring count.
    ///
    /// This scalar has a dedicated cache instead of a type-erased derived
    /// slot because it is commonly read immediately after parsing and its
    /// lookup cost is measurable through language bindings.
    #[doc(hidden)]
    #[inline(always)]
    pub fn derived_rdkit_num_rings(&self, compute: impl FnOnce() -> usize) -> usize {
        self.derived.rdkit_num_rings(compute)
    }

    /// Return the memoized value for `slot`, computing it with `compute` on
    /// first use. Intended for perception crates in this workspace; the slot
    /// set and its types are not a stable public API.
    #[doc(hidden)]
    pub fn derived<T, F>(
        &self,
        slot: crate::derived_cache::DerivedSlot,
        compute: F,
    ) -> std::sync::Arc<T>
    where
        T: std::any::Any + Send + Sync + std::panic::RefUnwindSafe + std::panic::UnwindSafe,
        F: FnOnce() -> T,
    {
        self.derived.get_or_compute(slot, compute)
    }

    /// Seed `slot` with a value computed on another molecule. Only valid when
    /// that value is exactly what the slot's own computation would produce on
    /// `self` (for example ring-bond flags of an index-aligned copy with the
    /// same graph and ring-eligible bonds).
    #[doc(hidden)]
    pub fn seed_derived<T>(&self, slot: crate::derived_cache::DerivedSlot, value: std::sync::Arc<T>)
    where
        T: std::any::Any + Send + Sync + std::panic::RefUnwindSafe + std::panic::UnwindSafe,
    {
        self.derived.seed(slot, value);
    }

    /// The value memoized for `slot`, if it has already been computed.
    #[doc(hidden)]
    pub fn derived_if_computed<T>(
        &self,
        slot: crate::derived_cache::DerivedSlot,
    ) -> Option<std::sync::Arc<T>>
    where
        T: std::any::Any + Send + Sync + std::panic::RefUnwindSafe + std::panic::UnwindSafe,
    {
        self.derived.peek(slot)
    }

    #[inline]
    fn invalidate_derived(&mut self) {
        self.derived = crate::derived_cache::DerivedCache::default();
    }

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

    /// Non-chemical caller tag for `idx`, or `None` when untagged or out of range.
    ///
    /// Tags are caller-managed labels, not unique identities: duplicates are
    /// allowed, including when combining independently tagged reactants. They
    /// do not affect chemical equality, canonicalization, or SMILES output.
    /// A SMILES write/parse round trip requires explicit remapping using the
    /// writer's atom-order result. Cloning an `Atom` alone does not copy a tag.
    /// Custom rebuilds must use [`MoleculeBuilder::set_tag`] or, when indices
    /// are unchanged, [`MoleculeBuilder::copy_atom_tags_from`]. Tag preservation
    /// is guaranteed for the operations listed on [`Self::set_tag`]; other
    /// molecule-rebuilding APIs are not covered by that contract.
    ///
    /// ```
    /// use chematic_core::{Atom, Element, MoleculeBuilder};
    /// let mut builder = MoleculeBuilder::new();
    /// let carbon = builder.add_atom(Atom::new(Element::C));
    /// let mut mol = builder.build();
    /// mol.set_tag(carbon, Some(42));
    /// assert_eq!(mol.clone().atom_tag(carbon).unwrap().get(), 42);
    /// mol.set_tag(carbon, None);
    /// assert_eq!(mol.atom_tag(carbon), None);
    /// ```
    pub fn atom_tag(&self, idx: AtomIdx) -> Option<core::num::NonZeroU16> {
        self.atom_tags.get(idx.0 as usize).copied().flatten()
    }

    /// R-group label attached to `idx`, if this wildcard came from `R`/`R<n>`.
    pub fn r_group_label(&self, idx: AtomIdx) -> Option<RGroupLabel> {
        self.r_groups.get(&idx.0).copied()
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

    /// The `(neighbor, bond)` adjacency list of atom `idx` as a slice, in the
    /// same order as [`Self::neighbors`].
    ///
    /// # Panics
    /// Panics if `idx` is out of range.
    #[inline]
    pub fn neighbor_slice(&self, idx: AtomIdx) -> &[(AtomIdx, BondIdx)] {
        &self.adjacency[idx.0 as usize]
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
