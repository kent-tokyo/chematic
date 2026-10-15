/// Structural errors in an enhanced stereo grouping. Indices are zero-based.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoGroupError {
    EmptyGroup {
        group: usize,
    },
    InvalidGroupNumber {
        group: usize,
    },
    DuplicateGroup {
        first: usize,
        second: usize,
    },
    InvalidAtom {
        group: usize,
        atom: AtomIdx,
    },
    InvalidBond {
        group: usize,
        bond: BondIdx,
    },
    AtomOverlap {
        first: usize,
        second: usize,
        atom: AtomIdx,
    },
    BondOverlap {
        first: usize,
        second: usize,
        bond: BondIdx,
    },
}
impl core::fmt::Display for StereoGroupError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyGroup { group } => {
                write!(f, "stereo group {group} has no atom or bond members")
            }
            Self::InvalidGroupNumber { group } => write!(
                f,
                "stereo group {group} requires a positive AND/OR group number"
            ),
            Self::DuplicateGroup { first, second } => write!(
                f,
                "stereo groups {first} and {second} repeat the same kind and number"
            ),
            Self::InvalidAtom { group, atom } => write!(
                f,
                "stereo group {group} references out-of-range atom {}",
                atom.0
            ),
            Self::InvalidBond { group, bond } => write!(
                f,
                "stereo group {group} references out-of-range bond {}",
                bond.0
            ),
            Self::AtomOverlap {
                first,
                second,
                atom,
            } => write!(
                f,
                "atom {} belongs to stereo groups {first} and {second}",
                atom.0
            ),
            Self::BondOverlap {
                first,
                second,
                bond,
            } => write!(
                f,
                "bond {} belongs to stereo groups {first} and {second}",
                bond.0
            ),
        }
    }
}
impl std::error::Error for StereoGroupError {}

impl Molecule {
    /// Replace enhanced atom/bond groups on an independent copy, after validating
    /// references, positive AND/OR numbers, unique group keys and disjoint
    /// membership. An empty list clears groups; an empty group is rejected.
    /// Duplicate members within one group are deduplicated in first-seen order.
    ///
    /// This checks metadata structure, not stereogenicity or 3D feasibility.
    /// It preserves the indexed graph, stored stereo frames and all metadata;
    /// derived caches start empty. No configuration is inverted or canonicalized.
    pub fn with_checked_stereo_groups(
        &self,
        groups: Vec<(StereoGroup, Vec<BondIdx>)>,
    ) -> Result<Self, StereoGroupError> {
        let mut kinds = std::collections::HashMap::new();
        let mut atoms = std::collections::HashMap::new();
        let mut bonds = std::collections::HashMap::new();
        for (i, (group, members)) in groups.iter().enumerate() {
            if group.atom_indices.is_empty() && members.is_empty() {
                return Err(StereoGroupError::EmptyGroup { group: i });
            }
            if matches!(
                group.kind,
                crate::StereoGroupKind::And(0) | crate::StereoGroupKind::Or(0)
            ) {
                return Err(StereoGroupError::InvalidGroupNumber { group: i });
            }
            if let Some(first) = kinds.insert(group.kind.clone(), i) {
                return Err(StereoGroupError::DuplicateGroup { first, second: i });
            }
            for &atom in &group.atom_indices {
                if atom.0 as usize >= self.atom_count() {
                    return Err(StereoGroupError::InvalidAtom { group: i, atom });
                }
                if let Some(first) = atoms.insert(atom, i)
                    && first != i
                {
                    return Err(StereoGroupError::AtomOverlap {
                        first,
                        second: i,
                        atom,
                    });
                }
            }
            for &bond in members {
                if bond.0 as usize >= self.bond_count() {
                    return Err(StereoGroupError::InvalidBond { group: i, bond });
                }
                if let Some(first) = bonds.insert(bond, i)
                    && first != i
                {
                    return Err(StereoGroupError::BondOverlap {
                        first,
                        second: i,
                        bond,
                    });
                }
            }
        }
        let groups = groups
            .into_iter()
            .map(|(g, bonds)| (StereoGroup::new(g.kind, g.atom_indices), bonds))
            .collect();
        let mut result = self.clone();
        result.set_stereo_groups_with_bonds(groups);
        Ok(result)
    }
}
