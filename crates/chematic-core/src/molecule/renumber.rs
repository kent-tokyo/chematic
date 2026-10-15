/// Invalid atom permutation or an invalid reference in caller-set metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RenumberError {
    Length { expected: usize, actual: usize },
    InvalidAtom(AtomIdx),
    DuplicateAtom(AtomIdx),
    InvalidMetadata { field: &'static str, index: u32 },
}

impl core::fmt::Display for RenumberError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Length { expected, actual } => {
                write!(f, "atom order needs {expected} indices, got {actual}")
            }
            Self::InvalidAtom(idx) => write!(f, "atom order index {} is out of range", idx.0),
            Self::DuplicateAtom(idx) => write!(f, "atom order repeats index {}", idx.0),
            Self::InvalidMetadata { field, index } => {
                write!(f, "invalid {field} reference {index}")
            }
        }
    }
}
impl std::error::Error for RenumberError {}

fn mapped_atom(
    inverse: &[AtomIdx],
    old: u32,
    field: &'static str,
) -> Result<AtomIdx, RenumberError> {
    inverse
        .get(old as usize)
        .copied()
        .ok_or(RenumberError::InvalidMetadata { field, index: old })
}

impl Molecule {
    /// Copy with `new_order[new_index] = old_index`, a full atom permutation.
    ///
    /// Bond indices and each atom's ordered neighbour slots remain unchanged.
    /// Atom references in stereo frames, groups and direction anchors are
    /// remapped, so chirality tags retain their meaning without parity flips.
    /// Atom data, tags, R labels and bond metadata survive; all derived caches
    /// start empty. No perception/sanitization is performed. Coordinates and
    /// external per-atom arrays must be reordered separately by the caller.
    pub fn renumber_atoms(&self, new_order: &[AtomIdx]) -> Result<Self, RenumberError> {
        let n = self.atoms.len();
        if new_order.len() != n {
            return Err(RenumberError::Length {
                expected: n,
                actual: new_order.len(),
            });
        }
        let mut inverse = vec![AtomIdx(STEREO_H_SENTINEL); n];
        for (new, &old) in new_order.iter().enumerate() {
            let slot = inverse
                .get_mut(old.0 as usize)
                .ok_or(RenumberError::InvalidAtom(old))?;
            if slot.0 != STEREO_H_SENTINEL {
                return Err(RenumberError::DuplicateAtom(old));
            }
            *slot = AtomIdx(new as u32);
        }
        // Molecule::clone deliberately clears every derived perception cache.
        let mut output = self.clone();
        output.atoms = new_order
            .iter()
            .map(|old| self.atoms[old.0 as usize].clone())
            .collect();
        if !self.atom_tags.is_empty() {
            output.atom_tags = new_order
                .iter()
                .map(|old| self.atom_tags.get(old.0 as usize).copied().flatten())
                .collect();
        }
        output.adjacency = new_order
            .iter()
            .map(|old| {
                self.adjacency[old.0 as usize]
                    .iter()
                    .map(|(nb, bond)| Ok((mapped_atom(&inverse, nb.0, "adjacency")?, *bond)))
                    .collect::<Result<Vec<_>, RenumberError>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        for bond in &mut output.bonds {
            bond.atom1 = mapped_atom(&inverse, bond.atom1.0, "bond endpoint")?;
            bond.atom2 = mapped_atom(&inverse, bond.atom2.0, "bond endpoint")?;
        }
        output.remap_renumber_metadata(&inverse)?;
        Ok(output)
    }

    fn remap_renumber_metadata(&mut self, inverse: &[AtomIdx]) -> Result<(), RenumberError> {
        self.validate_renumber_bond_metadata()?;
        self.r_groups = std::mem::take(&mut self.r_groups)
            .into_iter()
            .map(|(old, label)| Ok((mapped_atom(inverse, old, "R group")?.0, label)))
            .collect::<Result<_, RenumberError>>()?;
        for group in &mut self.stereo_groups {
            for idx in &mut group.atom_indices {
                *idx = mapped_atom(inverse, idx.0, "stereo group")?;
            }
        }
        self.stereo_neighbor_order = std::mem::take(&mut self.stereo_neighbor_order)
            .into_iter()
            .map(|(old, frame)| {
                let center = mapped_atom(inverse, old, "stereo center")?;
                let frame = frame
                    .into_iter()
                    .map(|nb| {
                        if nb == STEREO_H_SENTINEL {
                            Ok(nb)
                        } else {
                            mapped_atom(inverse, nb, "stereo neighbour").map(|a| a.0)
                        }
                    })
                    .collect::<Result<_, _>>()?;
                Ok((center.0, frame))
            })
            .collect::<Result<_, RenumberError>>()?;
        for anchor in self.bond_direction_anchors.values_mut() {
            *anchor = mapped_atom(inverse, anchor.0, "direction anchor")?;
        }
        Ok(())
    }

    fn validate_renumber_bond_metadata(&self) -> Result<(), RenumberError> {
        for (field, indices) in [
            (
                "stereo group bond",
                self.stereo_group_bonds
                    .iter()
                    .flatten()
                    .map(|b| b.0)
                    .collect::<Vec<_>>(),
            ),
            (
                "bond direction",
                self.bond_directions.keys().copied().collect(),
            ),
            (
                "direction anchor bond",
                self.bond_direction_anchors.keys().copied().collect(),
            ),
            (
                "ring closure bond",
                self.smiles_ring_closure_keys.keys().copied().collect(),
            ),
        ] {
            for index in indices {
                if index as usize >= self.bonds.len() {
                    return Err(RenumberError::InvalidMetadata { field, index });
                }
            }
        }
        Ok(())
    }
}
