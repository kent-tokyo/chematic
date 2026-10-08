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
mod molblock;
mod parse;
mod periodic;
mod rank;
mod sanitize;
mod stereo;
mod write;

use chematic_core::Molecule;

pub use molblock::{RdkitMolBlock, RdkitMolBlockAtom, RdkitMolBlockBond, rdkit_mol_block};

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

/// RDKit's `CalcNumAtomStereoCenters` and
/// `CalcNumUnspecifiedAtomStereoCenters` for `mol` as `MolFromSmiles`
/// leaves it: atoms flagged `_ChiralityPossible` by RDKit's legacy stereo
/// perception (RDKit 2026.03's default), and those of them without a
/// chiral tag.
pub fn rdkit_atom_stereocenter_counts(mol: &Molecule) -> Result<(usize, usize), RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    sanitize::remove_hs_and_sanitize(&mut m)?;
    stereo::legacy_stereo_perception(&mut m, true, true);
    let possible = m.atoms.iter().filter(|a| a.chirality_possible);
    let total = possible.clone().count();
    let unspecified = possible
        .filter(|a| a.chiral == mol::ChiralTag::Unspecified)
        .count();
    Ok((total, unspecified))
}

/// RDKit's legacy stereo perception (`MolFromSmiles`, RDKit 2026.03's
/// default) of `mol`, indexed like `mol`: per atom whether it keeps a chiral
/// tag and its `_CIPCode` (`b'R'`/`b'S'`), per bond its `Bond::BondStereo`
/// value (0 none, 1 any, 2 Z, 3 E).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitLegacyStereo {
    pub atom_tagged: Vec<bool>,
    pub atom_cip: Vec<Option<u8>>,
    pub bond_stereo: Vec<u8>,
}

/// [`RdkitLegacyStereo`] for `mol`; `Err` where the port cannot model the
/// molecule or RDKit's hydrogen removal would renumber its atoms.
pub fn rdkit_legacy_stereo(mol: &Molecule) -> Result<RdkitLegacyStereo, RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    let n_atoms = m.atoms.len();
    sanitize::remove_hs_and_sanitize(&mut m)?;
    if m.atoms.len() != n_atoms || m.bonds.len() != mol.bond_count() {
        return Err(RdkitSmilesError::Unsupported(
            "hydrogen removal renumbers atoms".into(),
        ));
    }
    stereo::legacy_stereo_perception(&mut m, true, true);
    let order = mol.rdkit_bond_order();
    let mut bond_stereo = vec![0u8; mol.bond_count()];
    for (k, b) in order.iter().enumerate() {
        bond_stereo[b.0 as usize] = m.bonds[k].stereo as u8;
    }
    Ok(RdkitLegacyStereo {
        atom_tagged: m
            .atoms
            .iter()
            .map(|a| a.chiral != mol::ChiralTag::Unspecified)
            .collect(),
        atom_cip: m.atoms.iter().map(|a| a.cip_code).collect(),
        bond_stereo,
    })
}

#[cfg(test)]
mod tests;
