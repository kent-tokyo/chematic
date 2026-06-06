//! `chematic-perception` — molecular perception algorithms.
//!
//! Provides:
//! - [`sssr`]: Smallest Set of Smallest Rings (SSSR) via Balducci-Pearlman algorithm.
//! - [`aromaticity`]: Hückel aromaticity perception for kekulized molecules.

#![forbid(unsafe_code)]

pub mod aromaticity;
pub mod sssr;

pub mod stereo2d;

pub use aromaticity::{AromaticityModel, apply_aromaticity, assign_aromaticity};
pub use chematic_core::{ValenceError, validate_valence};
pub use sssr::{RingSet, find_sssr};
pub use stereo2d::{StereoAssignment2D, apply_stereo_from_2d, assign_stereo_from_2d};

use chematic_core::Molecule;

/// Apply aromaticity to `mol` in-place (wrapper for [`apply_aromaticity`]).
pub fn aromatize(mol: &mut Molecule) {
    *mol = apply_aromaticity(mol);
}

/// Convert `mol` to Kekulé form in-place (wrapper for `kekulize` + `apply_kekule`).
///
/// Returns `Err` if kekulization fails (e.g. invalid aromatic system).
pub fn kekulize_inplace(mol: &mut Molecule) -> Result<(), chematic_core::KekuleError> {
    use chematic_core::{apply_kekule, kekulize};
    let result = kekulize(mol)?;
    *mol = apply_kekule(mol, &result);
    Ok(())
}
