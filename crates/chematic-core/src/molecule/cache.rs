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
    /// [`Self::fragments_with_source_atoms`] also returns, per fragment, the
    /// original index of each renumbered atom.
    pub fn fragments(&self) -> Vec<Molecule> {
        self.fragments_with_source_atoms()
            .into_iter()
            .map(|(fragment, _)| fragment)
            .collect()
    }

    /// [`Self::fragments`], also returning for each fragment the index in
    /// `self` of every fragment atom: `source[i]` is the original index of
    /// fragment atom `AtomIdx(i)`.
    ///
    /// Fragments are ordered by their lowest original atom index and keep
    /// the original relative atom order, so each `source` vector is strictly
    /// increasing and the vectors together partition `0..self.atom_count()`.
    pub fn fragments_with_source_atoms(&self) -> Vec<(Molecule, Vec<AtomIdx>)> {
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
                let mut source: Vec<AtomIdx> = Vec::new();
                for (aidx, atom) in self.atoms() {
                    if component[aidx.0 as usize] == cid {
                        let new_idx = builder.add_atom(atom.clone());
                        builder
                            .set_tag(new_idx, self.atom_tag(aidx).map(core::num::NonZeroU16::get));
                        if let Some(label) = self.r_group_label(aidx) {
                            builder.set_r_group(new_idx, label);
                        }
                        old_to_new.insert(aidx, new_idx);
                        source.push(aidx);
                    }
                }
                for (_, bond) in self.bonds() {
                    if let (Some(&a1), Some(&a2)) =
                        (old_to_new.get(&bond.atom1), old_to_new.get(&bond.atom2))
                    {
                        let _ = builder.add_bond(a1, a2, bond.order);
                    }
                }
                (builder.build(), source)
            })
            .collect()
    }
}
