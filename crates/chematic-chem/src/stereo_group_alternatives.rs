//! Bounded, correlated full alternatives of independent enhanced stereo groups.
use crate::{StereoGroupInversionError, invert_stereo_group};
use chematic_core::{Molecule, StereoGroupError, StereoGroupKind};

/// Hard allocation/search boundaries for the coordinate-free family API.
pub const MAX_STEREO_GROUP_ALTERNATIVES: usize = 64;
pub const MAX_STEREO_GROUP_ALTERNATIVE_ATOMS: usize = 1024;
pub const MAX_STEREO_GROUP_ALTERNATIVE_BONDS: usize = 4096;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoGroupAlternativesError {
    InvalidLimit,
    AtomLimit {
        actual: usize,
        limit: usize,
    },
    BondLimit {
        actual: usize,
        limit: usize,
    },
    TooManyRelativeGroups {
        actual: usize,
        limit: usize,
    },
    AlternativeLimit {
        required: usize,
        limit: usize,
    },
    InvalidGroups(StereoGroupError),
    Inversion {
        group: usize,
        cause: StereoGroupInversionError,
    },
    AllocationFailed,
}
impl core::fmt::Display for StereoGroupAlternativesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidLimit => write!(f, "max_alternatives must be between 1 and 64"),
            Self::AtomLimit { actual, limit } => write!(
                f,
                "stereo-group alternatives have {actual} atoms, exceeding {limit}"
            ),
            Self::BondLimit { actual, limit } => write!(
                f,
                "stereo-group alternatives have {actual} bonds, exceeding {limit}"
            ),
            Self::TooManyRelativeGroups { actual, limit } => {
                write!(f, "{actual} relative groups exceed the limit of {limit}")
            }
            Self::AlternativeLimit { required, limit } => write!(
                f,
                "{required} group alternatives exceed max_alternatives={limit}"
            ),
            Self::InvalidGroups(e) => e.fmt(f),
            Self::Inversion { group, cause } => write!(f, "stereo group {group}: {cause}"),
            Self::AllocationFailed => write!(f, "could not allocate stereo-group alternatives"),
        }
    }
}
impl std::error::Error for StereoGroupAlternativesError {}

/// Correlated drawn representatives, not deduplicated pure stereoisomers.
pub struct StereoGroupAlternatives {
    /// Original stored group indices for the AND/OR groups, in stored order.
    pub relative_group_indices: Vec<usize>,
    /// Output index is a bit mask: bit k selects full inversion of stored group
    /// `relative_group_indices[k]`. ABS and ungrouped configurations stay fixed.
    /// Every output retains the original enhanced-group metadata.
    pub molecules: Vec<Molecule>,
}

fn checked_source(mol: &Molecule) -> Result<Molecule, StereoGroupAlternativesError> {
    let groups = mol
        .stereo_groups()
        .iter()
        .enumerate()
        .map(|(i, group)| {
            (
                group.clone(),
                mol.stereo_group_bond_indices(i)
                    .unwrap_or_default()
                    .to_vec(),
            )
        })
        .collect();
    let mut result = mol
        .with_checked_stereo_groups(groups)
        .map_err(StereoGroupAlternativesError::InvalidGroups)?;
    // Retain legacy raw membership order/duplicates after structural validation.
    result.copy_stereo_groups_from(mol);
    // Use the same cold-cache/CIP policy for mask zero and changed representatives.
    for (idx, atom) in mol.atoms() {
        if atom.cip_code.is_some() {
            result.set_cip_code(idx, None);
        }
    }
    Ok(result)
}

fn check_limits(
    mol: &Molecule,
    max_alternatives: usize,
) -> Result<(), StereoGroupAlternativesError> {
    if !(1..=MAX_STEREO_GROUP_ALTERNATIVES).contains(&max_alternatives) {
        return Err(StereoGroupAlternativesError::InvalidLimit);
    }
    for (actual, limit, atoms) in [
        (mol.atom_count(), MAX_STEREO_GROUP_ALTERNATIVE_ATOMS, true),
        (mol.bond_count(), MAX_STEREO_GROUP_ALTERNATIVE_BONDS, false),
    ] {
        if actual > limit {
            return Err(if atoms {
                StereoGroupAlternativesError::AtomLimit { actual, limit }
            } else {
                StereoGroupAlternativesError::BondLimit { actual, limit }
            });
        }
    }
    Ok(())
}

/// Expand every AND/OR group by its full alternative, flipping all of that
/// group's tetrahedral and geometric double-bond members together. ABS groups
/// and ungrouped assignments remain fixed. Unspecified centers are not assigned.
///
/// At most six relative groups / 64 representatives / 1024 atoms / 4096 bonds.
/// The caller may lower `max_alternatives` (1..=64); incomplete results are never
/// returned. Unsupported members or any failed inversion return an error with
/// the original stored group index, leaving the input untouched.
///
/// Each result retains the enhanced record's groups, tags, indexed graph and
/// stereo frames; all stored CIP annotations are cleared and caches start cold.
/// There is no canonicalization, meso reduction, coordinate generation or
/// compound counting: symmetry-equivalent masks intentionally remain separate.
pub fn enumerate_stereo_group_alternatives(
    mol: &Molecule,
    max_alternatives: usize,
) -> Result<StereoGroupAlternatives, StereoGroupAlternativesError> {
    check_limits(mol, max_alternatives)?;
    let relative_group_indices: Vec<_> = mol
        .stereo_groups()
        .iter()
        .enumerate()
        .filter(|(_, group)| group.kind != StereoGroupKind::Absolute)
        .map(|(i, _)| i)
        .collect();
    if relative_group_indices.len() > 6 {
        return Err(StereoGroupAlternativesError::TooManyRelativeGroups {
            actual: relative_group_indices.len(),
            limit: 6,
        });
    }
    // The explicit six-group bound precedes this shift on every platform.
    let required = 1usize << relative_group_indices.len();
    if required > max_alternatives {
        return Err(StereoGroupAlternativesError::AlternativeLimit {
            required,
            limit: max_alternatives,
        });
    }
    let base = checked_source(mol)?;
    let mut molecules = Vec::new();
    molecules
        .try_reserve_exact(required)
        .map_err(|_| StereoGroupAlternativesError::AllocationFailed)?;
    molecules.push(base);
    for &group in &relative_group_indices {
        let prefix = molecules.len();
        for mask in 0..prefix {
            let inverted = invert_stereo_group(&molecules[mask], group)
                .map_err(|cause| StereoGroupAlternativesError::Inversion { group, cause })?;
            molecules.push(inverted);
        }
    }
    Ok(StereoGroupAlternatives {
        relative_group_indices,
        molecules,
    })
}
