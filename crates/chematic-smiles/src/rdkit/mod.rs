//! RDKit-compatible canonical SMILES.
//!
//! [`rdkit_canonical_smiles`] writes the string RDKit 2026.03.1's
//! `Chem.MolToSmiles(Chem.MolFromSmiles(s))` writes (isomeric, canonical,
//! default parameters, legacy stereo perception) for a molecule chematic
//! read from the SMILES `s`. It is separate from chematic's own canonical
//! SMILES ([`crate::canonical_smiles`]), which it does not change.
//!
//! The pipeline is a port of the RDKit C++ code, operation for operation,
//! on a model of RDKit's molecule ([`mol`]):
//!
//! 1. the parser's molecule (`toMol`: ring-closure bonds last, chiral tags
//!    relative to RDKit's bond order) — [`parse`];
//! 2. `MolOps::removeHs` and `MolOps::sanitizeMol` (`cleanUp`,
//!    `cleanUpOrganometallics`, valences, symmetrized SSSR, `Kekulize`,
//!    `assignRadicals`, `setAromaticity`, `setConjugation`,
//!    `setHybridization`, `cleanupChirality`, `adjustHs`) — [`sanitize`],
//!    [`kekulize`], [`aromaticity`];
//! 3. legacy `assignStereochemistry(cleanIt=true, force=true,
//!    flagPossibleStereoCenters=true)` — [`stereo`];
//! 4. `MolToSmiles`: per fragment `Canon::rankMolAtoms`,
//!    `Canon::canonicalizeFragment` and `FragmentSmilesConstruct`, fragments
//!    sorted and joined with `.` — [`rank`], [`canon`], [`write`].
//!
//! Inputs outside what the port models (non-tetrahedral chirality, query
//! bonds, molecules not read from SMILES, ...) are refused with
//! [`RdkitSmilesError::Unsupported`]; molecules RDKit itself would reject
//! while sanitizing are refused with [`RdkitSmilesError::Sanitization`].
//! No string is returned for them.

// The port keeps RDKit's index-based loops and its if/else-if chains (some
// arms share a body) so it reads side by side with the C++.
#![allow(clippy::needless_range_loop, clippy::if_same_then_else)]

mod aromaticity;
mod canon;
mod kekulize;
mod mol;
mod parse;
mod periodic;
mod rank;
mod sanitize;
mod stereo;
mod write;

use chematic_core::Molecule;

/// Why [`rdkit_canonical_smiles`] produced no string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitSmilesError {
    /// The molecule uses a feature the RDKit port does not model.
    Unsupported(String),
    /// RDKit's sanitization would reject the molecule (valence or
    /// kekulization failure), so `Chem.MolFromSmiles` returns `None`.
    Sanitization(String),
}

impl core::fmt::Display for RdkitSmilesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Unsupported(what) => {
                write!(f, "RDKit-compatible SMILES: unsupported input: {what}")
            }
            Self::Sanitization(what) => {
                write!(f, "RDKit-compatible SMILES: sanitization failed: {what}")
            }
        }
    }
}

impl std::error::Error for RdkitSmilesError {}

/// The canonical SMILES RDKit 2026.03.1 writes for `mol`:
/// `Chem.MolToSmiles(Chem.MolFromSmiles(s))` for the SMILES `s` chematic
/// parsed `mol` from (see the module documentation).
///
/// ```
/// let mol = chematic_smiles::parse("OC(=O)[C@@H]1CCCN1").unwrap();
/// assert_eq!(
///     chematic_smiles::rdkit_canonical_smiles(&mol).unwrap(),
///     "O=C(O)[C@@H]1CCCN1"
/// );
/// ```
pub fn rdkit_canonical_smiles(mol: &Molecule) -> Result<String, RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    sanitize::remove_hs_and_sanitize(&mut m)?;
    stereo::legacy_stereo_perception(&mut m, true, true);
    write::mol_to_smiles(&m)
}

#[cfg(test)]
mod tests;
