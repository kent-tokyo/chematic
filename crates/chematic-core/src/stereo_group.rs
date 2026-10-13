//! Enhanced stereochemistry group data model.
//!
//! Corresponds to the ChemDraw V3000 `M V30 STEREOGROUP` / `MDLV30/STE*`
//! collection entries, which distinguish three types of stereo specification:
//!
//! - [`StereoGroupKind::Absolute`] — the absolute configuration is known.
//! - [`StereoGroupKind::Or`] — the compound is one of two configurations (or
//!   a mixture), but it is not known which; all atoms in the same Or-group
//!   share the same ambiguity.
//! - [`StereoGroupKind::And`] — the compound is a mixture of the drawn
//!   configuration and its full inversion; all members vary together.

use crate::molecule::{AtomIdx, BondIdx};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Kind of an enhanced stereo group.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum StereoGroupKind {
    /// Absolute configuration — the drawn stereocenters are in the single
    /// stated configuration (no uncertainty).
    Absolute,
    /// "Or" group — the drawn configuration OR its inverse; group number
    /// identifies which Or-group this is (1, 2, ...).
    Or(u32),
    /// "And" group — independent mixture; group number identifies the
    /// And-group (1, 2, ...).
    And(u32),
}

/// A named collection of stereocenters and stereogenic bonds with a shared
/// ambiguity classification.
///
/// Stored in [`Molecule::stereo_groups`]. The member lists use 0-based core
/// indices. Bond members cover geometric double bonds and atropisomer axes as
/// represented by the V3000 `STEB*` collections introduced by CTfile 2026.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    pub kind: StereoGroupKind,
    pub atom_indices: Vec<AtomIdx>,
    pub bond_indices: Vec<BondIdx>,
}

impl StereoGroup {
    /// Construct a new stereo group.
    ///
    /// Duplicate atom indices are silently removed (RDKit PR #9258: duplicate
    /// indices in StereoGroups caused heap-use-after-free during removeHs).
    pub fn new(kind: StereoGroupKind, atom_indices: Vec<AtomIdx>) -> Self {
        Self::with_bonds(kind, atom_indices, Vec::new())
    }

    /// Construct a stereo group containing atom and/or bond members.
    ///
    /// Duplicate indices are removed independently while preserving the first
    /// occurrence. At least one member should be present; format parsers reject
    /// empty groups before calling this constructor.
    pub fn with_bonds(
        kind: StereoGroupKind,
        atom_indices: Vec<AtomIdx>,
        bond_indices: Vec<BondIdx>,
    ) -> Self {
        let mut seen = std::collections::HashSet::new();
        let atom_indices: Vec<AtomIdx> = atom_indices
            .into_iter()
            .filter(|idx| seen.insert(*idx))
            .collect();
        let mut seen = std::collections::HashSet::new();
        let bond_indices: Vec<BondIdx> = bond_indices
            .into_iter()
            .filter(|idx| seen.insert(*idx))
            .collect();
        StereoGroup {
            kind,
            atom_indices,
            bond_indices,
        }
    }
}
