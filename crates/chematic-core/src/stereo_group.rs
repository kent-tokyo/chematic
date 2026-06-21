//! Enhanced stereochemistry group data model.
//!
//! Corresponds to the ChemDraw V3000 `M V30 STEREOGROUP` / `MDLV30/STE*`
//! collection entries, which distinguish three types of stereo specification:
//!
//! - [`StereoGroupKind::Absolute`] — the absolute configuration is known.
//! - [`StereoGroupKind::Or`] — the compound is one of two configurations (or
//!   a mixture), but it is not known which; all atoms in the same Or-group
//!   share the same ambiguity.
//! - [`StereoGroupKind::And`] — the compound is a mixture of configurations;
//!   each atom may independently be R or S.

use crate::molecule::AtomIdx;

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

/// A named collection of stereocenters with a shared ambiguity classification.
///
/// Stored in [`Molecule::stereo_groups`].  The `atom_indices` list the
/// 0-based internal `AtomIdx` values of the stereocenters that belong to
/// this group.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    pub kind: StereoGroupKind,
    pub atom_indices: Vec<AtomIdx>,
}

impl StereoGroup {
    /// Construct a new stereo group.
    ///
    /// Duplicate atom indices are silently removed (RDKit PR #9258: duplicate
    /// indices in StereoGroups caused heap-use-after-free during removeHs).
    pub fn new(kind: StereoGroupKind, atom_indices: Vec<AtomIdx>) -> Self {
        let mut seen = std::collections::HashSet::new();
        let atom_indices: Vec<AtomIdx> = atom_indices
            .into_iter()
            .filter(|idx| seen.insert(*idx))
            .collect();
        StereoGroup { kind, atom_indices }
    }
}
