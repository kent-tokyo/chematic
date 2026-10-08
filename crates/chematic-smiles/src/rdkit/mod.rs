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
mod depict;
mod embed_view;
mod enumerate;
mod inchi_read;
mod kekulize;
mod mol;
mod molblock;
mod molblock2d;
mod parse;
mod periodic;
mod pyrandom;
mod rank;
mod sanitize;
mod stereo;
mod write;

use chematic_core::Molecule;

pub use embed_view::{RdkitMolView, RdkitViewAtom, RdkitViewBond, rdkit_mol_view};
pub use inchi_read::{InchiOutputAtom, InchiOutputStereo0D, rdkit_molecule_from_inchi_output};
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
    rdkit_smiles(mol, &RdkitSmilesParams::default())
}

/// RDKit's `SmilesWriteParams` as `Chem.MolToSmiles` exposes them
/// (`doRandom` is always false). The default is `MolToSmiles`'s.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RdkitSmilesParams {
    /// `isomericSmiles`: write chirality, `/` `\` and isotopes.
    pub isomeric: bool,
    /// `kekuleSmiles`: write the (canonical) Kekulé form.
    pub kekule: bool,
    /// `canonical`: canonical atom ranks and sorted fragments.
    pub canonical: bool,
    /// `allBondsExplicit`: write every bond symbol.
    pub all_bonds_explicit: bool,
    /// `allHsExplicit`: write every atom in brackets with its H count.
    pub all_hs_explicit: bool,
    /// `rootedAtAtom`: start the fragment holding this atom there.
    pub rooted_at_atom: Option<usize>,
}

impl Default for RdkitSmilesParams {
    fn default() -> Self {
        Self {
            isomeric: true,
            kekule: false,
            canonical: true,
            all_bonds_explicit: false,
            all_hs_explicit: false,
            rooted_at_atom: None,
        }
    }
}

/// `Chem.MolToSmiles(Chem.MolFromSmiles(s), **params)` for the SMILES `s`
/// chematic parsed `mol` from (RDKit 2026.03.1).
pub fn rdkit_smiles(
    mol: &Molecule,
    params: &RdkitSmilesParams,
) -> Result<String, RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    if has_added_hydrogens(mol) {
        sanitize::sanitize_keeping_hs(&mut m)?;
    } else {
        sanitize::remove_hs_and_sanitize(&mut m)?;
    }
    stereo::legacy_stereo_perception(&mut m, true, true);
    write::mol_to_smiles(&m, params)
}

/// Whether `mol` carries hydrogen atoms the SMILES parser cannot produce:
/// an H outside brackets (no H count), as `add_hydrogens` adds them. Such a
/// molecule stands for RDKit's molecule after `Chem.AddHs`, whose hydrogen
/// atoms are graph atoms that `MolToSmiles` writes; a molecule read from
/// SMILES only has bracket `[H]` atoms, which `MolFromSmiles` removes by
/// its `removeHs` rules.
fn has_added_hydrogens(mol: &Molecule) -> bool {
    mol.atoms().any(|(_, a)| {
        !a.wildcard && a.element == chematic_core::Element::H && a.hydrogen_count.is_none()
    })
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

/// The stereoisomers `rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers`
/// yields for `mol` (as `MolFromSmiles` reads it) with options
/// `onlyUnassigned=True, unique=True, maxIsomers=max_isomers,
/// tryEmbedding=False`: their sorted, distinct RDKit canonical SMILES.
///
/// `Err` where the port cannot model the molecule (including enhanced
/// stereo groups), or where RDKit would pick a random sample because there
/// are more flip combinations than `max_isomers`.
pub fn rdkit_stereoisomer_smiles(
    mol: &Molecule,
    max_isomers: usize,
) -> Result<Vec<String>, RdkitSmilesError> {
    if !mol.stereo_groups().is_empty() {
        return Err(RdkitSmilesError::Unsupported(
            "enhanced stereo groups".into(),
        ));
    }
    let mut m = parse::from_chematic(mol)?;
    sanitize::remove_hs_and_sanitize(&mut m)?;
    stereo::legacy_stereo_perception(&mut m, true, true);
    for a in &mut m.atoms {
        a.cip_code = None;
    }
    enumerate::enumerate(&m, max_isomers)
}

/// The 2D coordinates RDKit 2026.03.1's default depiction gives the atoms
/// of `Chem.MolFromSmiles(s)`: `rdDepictor.Compute2DCoords(mol)` (RDKit's
/// own depictor, not CoordGen; `canonOrient=True`, no coordinate map, no
/// random sampling, no ring templates), for the SMILES `s` chematic parsed
/// `mol` from. One `[x, y]` per atom of RDKit's molecule, in RDKit's atom
/// order (chematic's atom order when `MolFromSmiles` removes no hydrogen
/// atoms).
///
/// `MolFromSmiles` is modelled as in [`rdkit_canonical_smiles`]; the
/// depiction is a port of RDKit's `compute2DCoords` that reproduces its
/// floating-point operations in order. Exact bits are platform-dependent
/// because RDKit and chematic call the host libm; portable validation uses an
/// absolute coordinate tolerance of `1e-12`.
///
/// ```
/// let mol = chematic_smiles::parse("CCO").unwrap();
/// let xy = chematic_smiles::rdkit_2d_coords(&mol).unwrap();
/// assert_eq!(xy.len(), 3);
/// assert_eq!(xy[1], [0.0, 0.5000000000000001]);
/// ```
pub fn rdkit_2d_coords(mol: &Molecule) -> Result<Vec<[f64; 2]>, RdkitSmilesError> {
    let (m, cip) = rdkit_mol_from_smiles(mol)?;
    depict::compute_2d_coords(&m, &cip)
}

/// `Chem.MolToMolBlock(m)` after `rdDepictor.Compute2DCoords(m)` for
/// `m = Chem.MolFromSmiles(s)`, the SMILES `s` chematic parsed `mol` from
/// (RDKit 2026.03.1): the V2000 MOL block with the coordinates of
/// [`rdkit_2d_coords`], RDKit's kekulization, wedge/hash bonds chosen by
/// `pickBondsToWedge`, crossed double bonds and `M  CHG`/`RAD`/`ISO` lines.
///
/// Molecules RDKit would write as V3000 (dative bonds, more than 999 atoms
/// or bonds) are refused. A dummy atom without map number, isotope, charge
/// or hydrogens is written as RDKit writes a bare `*` (chematic does not
/// keep whether it was bracketed).
///
/// ```
/// let mol = chematic_smiles::parse("C[C@H](O)F").unwrap();
/// let block = chematic_smiles::rdkit_mol_block_2d(&mol).unwrap();
/// assert!(block.starts_with("\n     RDKit          2D\n\n  4  3  0"));
/// assert!(block.contains("  2  1  1  1\n"));
/// ```
pub fn rdkit_mol_block_2d(mol: &Molecule) -> Result<String, RdkitSmilesError> {
    let (m, cip) = rdkit_mol_from_smiles(mol)?;
    let xy = depict::compute_2d_coords(&m, &cip)?;
    let bare_dummies: Vec<bool> = m
        .atoms
        .iter()
        .map(|a| {
            a.anum == 0
                && a.map.is_none()
                && a.isotope == 0
                && a.charge == 0
                && a.num_explicit_hs == 0
        })
        .collect();
    molblock2d::mol_block_2d(&m, &cip, &xy, &bare_dummies)
}

/// `Chem.MolFromSmiles(s)` for the SMILES `s` chematic parsed `mol` from,
/// with the atoms' `_CIPRank` values (empty when RDKit sets none).
fn rdkit_mol_from_smiles(mol: &Molecule) -> Result<(mol::Mol, Vec<u32>), RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    if has_added_hydrogens(mol) {
        sanitize::sanitize_keeping_hs(&mut m)?;
    } else {
        sanitize::remove_hs_and_sanitize(&mut m)?;
    }
    let cip = stereo::legacy_stereo_perception(&mut m, true, true);
    Ok((m, cip))
}

#[cfg(test)]
mod depict_tests;
#[cfg(test)]
mod tests;
