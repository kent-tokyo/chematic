//! Conservative proof-backed reduction of a relative stereo group's full alternative.
//! Exact colored graph mappings, not canonical-string equality, justify promotion.
use crate::{StereoGroupAlternativesError, enumerate_stereo_group_alternatives};
use chematic_core::{AtomIdx, BondIdx, Molecule, StereoGroup, StereoGroupKind};
mod graph;
mod search;

pub const MAX_MESO_GROUP_ATOMS: usize = 64;
pub const MAX_MESO_GROUP_BONDS: usize = 128;
pub const MAX_MESO_GROUP_VISITS: u64 = 10_000_000;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MesoStereoGroupError {
    GroupIndex { index: usize, groups: usize },
    AbsoluteGroup,
    GraphLimit,
    InvalidVisitLimit,
    BudgetExceeded { visits: u64, limit: u64 },
    Alternatives(StereoGroupAlternativesError),
    UnsupportedAtom(AtomIdx),
    IncompleteFrame(AtomIdx),
    UnspecifiedTetrahedralAtom(AtomIdx),
    UnsupportedBond(BondIdx),
    UnsupportedDirections,
}
impl core::fmt::Display for MesoStereoGroupError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::GroupIndex { index, groups } => {
                write!(f, "group index {index} is outside {groups} groups")
            }
            Self::AbsoluteGroup => write!(f, "ABS is already absolute"),
            Self::GraphLimit => write!(
                f,
                "meso-group proof supports at most 64 atoms and 128 bonds"
            ),
            Self::InvalidVisitLimit => write!(f, "max_visits must be between 1 and 10000000"),
            Self::BudgetExceeded { visits, limit } => write!(
                f,
                "meso-group proof exhausted {visits} candidate visits (limit {limit}); no reduction returned"
            ),
            Self::Alternatives(e) => e.fmt(f),
            Self::UnsupportedAtom(a) => write!(
                f,
                "atom {} has unsupported wildcard/axial/non-tetrahedral stereo",
                a.0
            ),
            Self::IncompleteFrame(a) => write!(
                f,
                "atom {} requires a complete unique four-ligand frame",
                a.0
            ),
            Self::UnspecifiedTetrahedralAtom(a) => write!(
                f,
                "atom {} has an unspecified locally tetrahedral configuration",
                a.0
            ),
            Self::UnsupportedBond(b) => write!(
                f,
                "bond {} has unsupported query/dative/zero order, a loop or parallel edge",
                b.0
            ),
            Self::UnsupportedDirections => write!(
                f,
                "unsupported, conflicting or unaccounted geometric directions"
            ),
        }
    }
}
impl std::error::Error for MesoStereoGroupError {}

/// A checked mapping from the lower-mask configuration to its full alternative.
pub struct MesoStereoWitness {
    pub mask: usize,
    /// Source atom index -> target atom index; preserves groups and source labels.
    pub atom_permutation: Vec<AtomIdx>,
}
pub struct MesoStereoGroupReduction {
    pub molecule: Molecule,
    /// True only after every other-group assignment has a checked witness.
    pub reduced: bool,
    pub checked_pairs: usize,
    pub search_visits: u64,
    pub witnesses: Vec<MesoStereoWitness>,
}

fn promote(mol: &Molecule, selected: usize) -> Molecule {
    let mut groups: Vec<(StereoGroup, Vec<BondIdx>)> = mol
        .stereo_groups()
        .iter()
        .enumerate()
        .map(|(i, g)| {
            (
                g.clone(),
                mol.stereo_group_bond_indices(i)
                    .unwrap_or_default()
                    .to_vec(),
            )
        })
        .collect();
    if let Some(abs) = groups
        .iter()
        .position(|(g, _)| g.kind == StereoGroupKind::Absolute)
    {
        let (relative, bonds) = groups[selected].clone();
        groups[abs].0.atom_indices.extend(relative.atom_indices);
        groups[abs].1.extend(bonds);
        groups.remove(selected);
    } else {
        groups[selected].0.kind = StereoGroupKind::Absolute;
    }
    // Existing structural validation established disjoint membership. The only
    // change merges the selected members into the unique ABS collection.
    let mut result = mol.clone();
    result.set_stereo_groups_with_bonds(groups);
    for (idx, atom) in mol.atoms() {
        if atom.cip_code.is_some() {
            result.set_cip_code(idx, None);
        }
    }
    result
}

/// Promote one AND/OR group to ABS only when its full alternative is equivalent
/// under a configuration-preserving graph mapping for every assignment of all
/// other relative groups. Group membership, atom maps and atom tags are colors:
/// a proof cannot erase their distinctions or exchange independent groups.
///
/// Returns an unchanged cold-cache copy when no mapping is found for a pair;
/// this is conservative retention, not a general "non-meso" classification.
/// Budget exhaustion is an error, never a negative result or partial promotion.
/// No drawing, indexed graph or tetrahedral/geometric configuration is changed.
pub fn reduce_meso_stereo_group(
    mol: &Molecule,
    index: usize,
    max_visits: u64,
) -> Result<MesoStereoGroupReduction, MesoStereoGroupError> {
    if !(1..=MAX_MESO_GROUP_VISITS).contains(&max_visits) {
        return Err(MesoStereoGroupError::InvalidVisitLimit);
    }
    if mol.atom_count() > MAX_MESO_GROUP_ATOMS || mol.bond_count() > MAX_MESO_GROUP_BONDS {
        return Err(MesoStereoGroupError::GraphLimit);
    }
    let group = mol
        .stereo_groups()
        .get(index)
        .ok_or(MesoStereoGroupError::GroupIndex {
            index,
            groups: mol.stereo_groups().len(),
        })?;
    if group.kind == StereoGroupKind::Absolute {
        return Err(MesoStereoGroupError::AbsoluteGroup);
    }
    let alternatives =
        enumerate_stereo_group_alternatives(mol, 64).map_err(MesoStereoGroupError::Alternatives)?;
    let bit = alternatives
        .relative_group_indices
        .iter()
        .position(|&i| i == index)
        .expect("validated relative group");
    let mut budget = search::Budget {
        visits: 0,
        limit: max_visits,
    };
    let mut witnesses = Vec::new();
    let mut checked_pairs = 0;
    for mask in 0..alternatives.molecules.len() {
        if mask & (1 << bit) != 0 {
            continue;
        }
        let a = graph::Graph::new(&alternatives.molecules[mask])?;
        let b = graph::Graph::new(&alternatives.molecules[mask | (1 << bit)])?;
        checked_pairs += 1;
        let Some(mapping) = search::find_mapping(&a, &b, &mut budget)? else {
            return Ok(MesoStereoGroupReduction {
                molecule: alternatives.molecules[0].clone(),
                reduced: false,
                checked_pairs,
                search_visits: budget.visits,
                witnesses,
            });
        };
        witnesses.push(MesoStereoWitness {
            mask,
            atom_permutation: mapping.into_iter().map(|i| AtomIdx(i as u32)).collect(),
        });
    }
    Ok(MesoStereoGroupReduction {
        molecule: promote(mol, index),
        reduced: true,
        checked_pairs,
        search_visits: budget.visits,
        witnesses,
    })
}
