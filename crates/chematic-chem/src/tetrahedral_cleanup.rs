//! Explicit conservative stereo cleanup motivated by RDKit issue #9656.
use crate::meso_stereo_group::{self, MesoStereoGroupError};
use chematic_core::{AtomIdx, Chirality, Molecule};

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum TetrahedralCleanupError {
    InvalidAtom { atom: AtomIdx, atoms: usize },
    UnspecifiedAtom(AtomIdx),
    EnhancedGroups,
    GraphLimit,
    InvalidVisitLimit,
    Proof(MesoStereoGroupError),
}
impl core::fmt::Display for TetrahedralCleanupError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidAtom { atom, atoms } => {
                write!(f, "cleanup atom {} is outside {atoms} atoms", atom.0)
            }
            Self::UnspecifiedAtom(a) => write!(
                f,
                "cleanup atom {} requires specified tetrahedral chirality",
                a.0
            ),
            Self::EnhancedGroups => write!(
                f,
                "tetrahedral cleanup does not support enhanced stereo groups"
            ),
            Self::GraphLimit => write!(
                f,
                "tetrahedral cleanup supports at most 64 atoms and 128 bonds"
            ),
            Self::InvalidVisitLimit => write!(f, "max_visits must be between 1 and 10000000"),
            Self::Proof(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for TetrahedralCleanupError {}

pub struct TetrahedralStereoCleanup {
    pub molecule: Molecule,
    pub atom: AtomIdx,
    /// True only after an exact anchored inversion witness has been found.
    pub removed: bool,
    pub search_visits: u64,
    /// Original -> inverted molecule atom mapping, with the target fixed.
    pub witness: Option<Vec<AtomIdx>>,
}

/// Remove a specified tetrahedral declaration only after proving it redundant.
///
/// Invert this one center in an independent copy and search for an exact graph
/// mapping back to the original that FIXES the selected atom and preserves all
/// other stored tetrahedral/geometric configurations, element/isotope/charge/H,
/// aromaticity, bond orders, atom maps and tags. A successful mapping necessarily
/// permutes the selected center's four-ligand frame oddly. No CIP-label or
/// canonical-string equality is used as a proof. On success only this tag/frame
/// is removed, and all atom CIP annotations are cleared for recomputation.
///
/// No witness means conservative retention, not a stereogenicity classification.
/// Errors/budget exhaustion never return a partial cleanup. Source is unchanged.
/// The strict proof profile rejects enhanced groups, wildcard/axial/non-tetrahedral
/// stereo, incomplete frames, unspecified locally tetrahedral atoms and unsupported
/// bond directions/types. Inputs are capped at 64 atoms/128 bonds before search.
/// This opt-in API does not change parsing/writing defaults or generate coordinates.
pub fn remove_redundant_tetrahedral_stereo(
    mol: &Molecule,
    atom: AtomIdx,
    max_visits: u64,
) -> Result<TetrahedralStereoCleanup, TetrahedralCleanupError> {
    if !(1..=meso_stereo_group::MAX_MESO_GROUP_VISITS).contains(&max_visits) {
        return Err(TetrahedralCleanupError::InvalidVisitLimit);
    }
    if mol.atom_count() > meso_stereo_group::MAX_MESO_GROUP_ATOMS
        || mol.bond_count() > meso_stereo_group::MAX_MESO_GROUP_BONDS
    {
        return Err(TetrahedralCleanupError::GraphLimit);
    }
    if atom.0 as usize >= mol.atom_count() {
        return Err(TetrahedralCleanupError::InvalidAtom {
            atom,
            atoms: mol.atom_count(),
        });
    }
    if !mol.stereo_groups().is_empty() {
        return Err(TetrahedralCleanupError::EnhancedGroups);
    }
    let opposite = match mol.atom(atom).chirality {
        Chirality::Clockwise => Chirality::CounterClockwise,
        Chirality::CounterClockwise => Chirality::Clockwise,
        _ => return Err(TetrahedralCleanupError::UnspecifiedAtom(atom)),
    };
    let mut inverted = mol.clone();
    inverted.set_chirality(atom, opposite);
    let (witness, search_visits) =
        meso_stereo_group::anchored_inversion_witness(mol, &inverted, atom, max_visits)
            .map_err(TetrahedralCleanupError::Proof)?;
    let mut molecule = mol.clone();
    if witness.is_some() {
        molecule.set_chirality(atom, Chirality::None);
        molecule.clear_stereo_neighbor_order(atom);
        for (idx, a) in mol.atoms() {
            if a.cip_code.is_some() {
                molecule.set_cip_code(idx, None);
            }
        }
    }
    Ok(TetrahedralStereoCleanup {
        molecule,
        atom,
        removed: witness.is_some(),
        search_visits,
        witness,
    })
}
