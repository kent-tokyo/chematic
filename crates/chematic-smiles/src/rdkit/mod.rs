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

mod align;
mod aromaticity;
mod canon;
mod crippen_params;
mod depict;
mod embed_view;
mod enumerate;
mod findstereo;
mod inchi_read;
mod kekulize;
mod matrices;
mod mol;
mod mol2_read;
mod molblock;
mod molblock2d;
mod molhash;
mod murcko;
mod parse;
mod pdb;
mod pdb_read;
mod periodic;
mod pyrandom;
mod rank;
mod rxn_smarts;
mod sanitize;
mod smarts_match;
mod smarts_write;
mod stereo;
mod substruct;
mod tautomer;
mod write;
mod xyz_read;

use chematic_core::Molecule;

pub use embed_view::{RdkitMolView, RdkitViewAtom, RdkitViewBond, rdkit_mol_view};
pub use inchi_read::{InchiOutputAtom, InchiOutputStereo0D, rdkit_molecule_from_inchi_output};
pub use molblock::{RdkitMolBlock, RdkitMolBlockAtom, RdkitMolBlockBond, rdkit_mol_block};
pub use molhash::RdkitHashFunction;
pub use tautomer::RdkitTautomerStatus;

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
    stereo::legacy_stereo_perception_unflagged(&mut m);
    write::mol_to_smiles_owned(m, params)
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

/// What `MolFromSmiles`' sanitization changes in the bonds and charges of
/// a molecule before perception (`cleanUp`: nitro-like nitrogens, P(=O)=C,
/// halogen oxides; `cleanUpOrganometallics`: dative bonds from hypervalent
/// atoms to metals), indexed like `mol`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RdkitCleanupEdits {
    /// `(atom, new formal charge)`.
    pub charges: Vec<(usize, i32)>,
    /// `(bond, new order, dative donor atom)`.
    pub bonds: Vec<(usize, chematic_core::BondOrder, Option<usize>)>,
}

/// [`RdkitCleanupEdits`] for `mol` (`None`: nothing changes); `Err` where the
/// port cannot model the molecule.
pub fn rdkit_cleanup_edits(mol: &Molecule) -> Result<Option<RdkitCleanupEdits>, RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    for a in 0..m.atoms.len() {
        m.update_atom_property_cache(a, false)?;
    }
    let before = m.clone();
    sanitize::clean_up(&mut m)?;
    sanitize::clean_up_organometallics(&mut m)?;
    let mut edits = RdkitCleanupEdits::default();
    for (a, (x, y)) in before.atoms.iter().zip(&m.atoms).enumerate() {
        if x.charge != y.charge {
            edits.charges.push((a, y.charge));
        }
    }
    let order = mol.rdkit_bond_order();
    for (k, (x, y)) in before.bonds.iter().zip(&m.bonds).enumerate() {
        if x.bt == y.bt {
            continue;
        }
        let bo = match y.bt {
            mol::BondType::Single => chematic_core::BondOrder::Single,
            mol::BondType::Double => chematic_core::BondOrder::Double,
            mol::BondType::Triple => chematic_core::BondOrder::Triple,
            mol::BondType::Dative => chematic_core::BondOrder::Dative,
            _ => return Err(RdkitSmilesError::Unsupported("clean-up bond type".into())),
        };
        let donor = (y.bt == mol::BondType::Dative).then_some(y.begin);
        edits.bonds.push((order[k].0 as usize, bo, donor));
    }
    Ok((!edits.charges.is_empty() || !edits.bonds.is_empty()).then_some(edits))
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

/// `Chem.MolToPDBBlock(m)` for `m = Chem.MolFromSmiles(s)` (RDKit
/// 2026.03.1, default flavor): `HETATM` lines named per element (`C1`,
/// `C2`, ...) in residue `UNL`, then `CONECT` records of the canonical
/// Kekulé structure (bond orders as repeated entries). `coords` (one per
/// atom of RDKit's molecule) stands for a conformer; without it RDKit
/// writes zero coordinates.
pub fn rdkit_pdb_block(
    mol: &Molecule,
    coords: Option<&[[f64; 3]]>,
) -> Result<String, RdkitSmilesError> {
    let m = rdkit_mol_for_writing(mol)?;
    pdb::mol_to_pdb_block(&m, coords)
}

/// `Chem.MolToSmarts(m, isomericSmiles, rootedAtAtom)` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1): atoms as `[#n]`/`[Sym]`
/// with isotope, chirality (and `H` on chiral atoms with one explicit H),
/// charge and map number, every bond explicit, input atom order (no
/// canonicalization).
pub fn rdkit_smarts(
    mol: &Molecule,
    isomeric: bool,
    rooted_at_atom: Option<usize>,
) -> Result<String, RdkitSmilesError> {
    let m = rdkit_mol_for_writing(mol)?;
    Ok(smarts_write::mol_to_smarts_owned(m, isomeric, rooted_at_atom, true)?.0)
}

/// `Chem.MolToCXSmiles(m, params)` for `m = Chem.MolFromSmiles(s)` (RDKit
/// 2026.03.1, `CXSmilesFields.CX_ALL`): the SMILES written with dative
/// bonds as single bonds, then the CXSMILES extension of a molecule read
/// from SMILES (radicals `^n:`, ring double bonds of unknown geometry
/// `ctu:`, coordinate bonds `C:`). With `params.kekule` the whole molecule
/// is kekulized first, as `MolToCXSmiles` does. Molecules with enhanced
/// stereo groups are refused.
///
/// ```
/// let mol = chematic_smiles::parse("C[CH2]").unwrap();
/// assert_eq!(
///     chematic_smiles::rdkit_cx_smiles(&mol, &Default::default()).unwrap(),
///     "[CH2]C |^1:0|"
/// );
/// ```
pub fn rdkit_cx_smiles(
    mol: &Molecule,
    params: &RdkitSmilesParams,
) -> Result<String, RdkitSmilesError> {
    if !mol.stereo_groups().is_empty() {
        return Err(RdkitSmilesError::Unsupported(
            "enhanced stereo groups".into(),
        ));
    }
    let mut m = rdkit_mol_for_writing(mol)?;
    let mut p = *params;
    if p.kekule {
        kekulize::kekulize(&mut m)?;
        p.kekule = false;
    }
    write::mol_to_cx_smiles(m, &p, false)
}

/// `Chem.MolFragmentToSmiles(m, atomsToUse, bondsToUse, **params)` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1, no atom or bond symbols):
/// the atoms `atoms` (indices into RDKit's molecule) with the bonds
/// `bonds` (default: every bond between two of them), ranked by
/// `Canon::rankFragmentAtoms`; chiral atoms with a bond outside the
/// fragment lose their tag, disconnected pieces are joined with `.`.
/// `params.kekule` is refused for fragments with aromatic atoms.
///
/// ```
/// let mol = chematic_smiles::parse("OC(=O)c1ccccc1").unwrap();
/// let p = Default::default();
/// assert_eq!(chematic_smiles::rdkit_fragment_smiles(&mol, &[0, 1, 2], None, &p).unwrap(), "O=CO");
/// ```
pub fn rdkit_fragment_smiles(
    mol: &Molecule,
    atoms: &[usize],
    bonds: Option<&[usize]>,
    params: &RdkitSmilesParams,
) -> Result<String, RdkitSmilesError> {
    let m = rdkit_mol_for_writing(mol)?;
    write::mol_fragment_to_smiles(m, atoms, bonds, params)
}

/// `Chem.MolFromSmiles(s).GetNumAtoms()` for the SMILES `s` chematic
/// parsed `mol` from: the atom count after RDKit's hydrogen removal (the
/// index range of RDKit's atom numbering).
pub fn rdkit_num_atoms(mol: &Molecule) -> Result<usize, RdkitSmilesError> {
    Ok(rdkit_mol_for_writing(mol)?.atoms.len())
}

/// `Chem.GetDistanceMatrix(m, useBO, useAtomWts)` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1; the `Chem.AddHs` molecule
/// when `mol` carries added hydrogens): topological distances by
/// Floyd-Warshall, `1e8` between fragments; `use_bo` weights bonds by
/// `1 / bond order` (aromatic `2/3`), `use_atom_wts` sets the diagonal to
/// `6 / atomic number`. One row per atom in RDKit's atom order.
pub fn rdkit_distance_matrix(
    mol: &Molecule,
    use_bo: bool,
    use_atom_wts: bool,
) -> Result<Vec<Vec<f64>>, RdkitSmilesError> {
    let m = rdkit_mol_for_writing(mol)?;
    let n = m.atoms.len();
    let d = matrices::distance_mat(&m, use_bo, use_atom_wts);
    Ok(d.chunks(n.max(1)).map(<[f64]>::to_vec).take(n).collect())
}

/// `Chem.Get3DDistanceMatrix(m, useAtomWts=...)` for RDKit's molecule of
/// `mol` (as in [`rdkit_distance_matrix`]) with the conformer `coords`
/// (one `[x, y, z]` per atom in RDKit's atom order).
pub fn rdkit_distance_matrix_3d(
    mol: &Molecule,
    coords: &[[f64; 3]],
    use_atom_wts: bool,
) -> Result<Vec<Vec<f64>>, RdkitSmilesError> {
    let m = rdkit_mol_for_writing(mol)?;
    let n = m.atoms.len();
    if coords.len() != n {
        return Err(RdkitSmilesError::Unsupported(format!(
            "{} coordinates for {n} atoms",
            coords.len()
        )));
    }
    let d = matrices::distance_mat_3d(&m, coords, use_atom_wts);
    Ok(d.chunks(n.max(1)).map(<[f64]>::to_vec).take(n).collect())
}

/// One problem `MolOps::detectChemistryProblems` reports (see
/// [`rdkit_detect_chemistry_problems`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitChemistryProblem {
    /// `MolSanitizeException::getType()`: `"AtomValenceException"`,
    /// `"AtomKekulizeException"` or `"KekulizeException"`.
    pub kind: String,
    /// The atom (`getAtomIdx()`) or atoms (`getAtomIndices()`) involved,
    /// in the parser's atom order.
    pub atoms: Vec<usize>,
    /// RDKit's message.
    pub message: String,
}

/// `Chem.DetectChemistryProblems(Chem.MolFromSmiles(s, sanitize=False))`
/// (RDKit 2026.03.1, `SANITIZE_ALL`) for the SMILES `s` chematic parsed
/// `mol` from: `cleanUp`, a strict `updatePropertyCache` per atom (each
/// failure an `AtomValenceException`), then a non-canonical `Kekulize` on
/// SSSR rings (`AtomKekulizeException` for a non-ring aromatic atom,
/// `KekulizeException` with the unkekulized atoms). Problems come in
/// RDKit's order. Parse SMILES with impossible O/F valences with
/// [`crate::parse_template`], which [`crate::parse`] rejects.
///
/// ```
/// let mol = chematic_smiles::parse("c1cccc1").unwrap();
/// let p = chematic_smiles::rdkit_detect_chemistry_problems(&mol).unwrap();
/// assert_eq!(p[0].kind, "KekulizeException");
/// assert_eq!(p[0].atoms, vec![0, 1, 2, 3, 4]);
/// ```
pub fn rdkit_detect_chemistry_problems(
    mol: &Molecule,
) -> Result<Vec<RdkitChemistryProblem>, RdkitSmilesError> {
    detect_chemistry_problems(parse::from_chematic(mol)?)
}

/// `Chem.DetectChemistryProblems(Chem.MolFromSmiles(s))`: as
/// [`rdkit_detect_chemistry_problems`] on RDKit's sanitized molecule (atom
/// indices after its hydrogen removal); `Err` where `MolFromSmiles` itself
/// fails.
pub fn rdkit_detect_chemistry_problems_sanitized(
    mol: &Molecule,
) -> Result<Vec<RdkitChemistryProblem>, RdkitSmilesError> {
    // `clearComputedProps` leaves the ring information (symmetrized SSSR);
    // `Kekulize` runs `findSSSR` only without it.
    detect_chemistry_problems(rdkit_mol_for_writing(mol)?)
}

fn detect_chemistry_problems(
    mut m: mol::Mol,
) -> Result<Vec<RdkitChemistryProblem>, RdkitSmilesError> {
    let mut res = Vec::new();
    sanitize::clean_up(&mut m)?;
    for a in 0..m.atoms.len() {
        if let Err(e) = m.update_atom_property_cache(a, true) {
            let message = match e {
                RdkitSmilesError::Sanitization(msg) => msg,
                other => return Err(other),
            };
            res.push(RdkitChemistryProblem {
                kind: "AtomValenceException".into(),
                atoms: vec![a],
                message,
            });
        }
    }
    if m.rings.is_none()
        && (m.bonds.iter().any(|b| b.aromatic) || m.atoms.iter().any(|a| a.aromatic))
    {
        m.find_sssr()?;
    }
    if let Err(e) = kekulize::kekulize(&mut m) {
        let message = match e {
            RdkitSmilesError::Sanitization(msg) => msg,
            other => return Err(other),
        };
        let numbers = |s: &str| -> Vec<usize> {
            s.split(|c: char| !c.is_ascii_digit())
                .filter_map(|t| t.parse().ok())
                .collect()
        };
        let problem =
            if let Some(rest) = message.strip_prefix("Can't kekulize mol.  Unkekulized atoms:") {
                RdkitChemistryProblem {
                    kind: "KekulizeException".into(),
                    atoms: numbers(rest),
                    message,
                }
            } else if let Some(rest) = message
                .strip_prefix("non-ring atom ")
                .or_else(|| message.strip_prefix("Kekulization somehow screwed up valence on "))
            {
                RdkitChemistryProblem {
                    kind: "AtomKekulizeException".into(),
                    atoms: numbers(rest).into_iter().take(1).collect(),
                    message,
                }
            } else {
                RdkitChemistryProblem {
                    kind: "AtomValenceException".into(),
                    atoms: numbers(&message).into_iter().take(1).collect(),
                    message,
                }
            };
        res.push(problem);
    }
    Ok(res)
}

/// `Chem.MolToSmarts(Chem.MolFromSmarts(smarts))` (RDKit 2026.03.1): the
/// query re-written as RDKit writes query molecules (`[OH]` becomes
/// `[O&H1]`, `;` becomes `&` where no `,` needs it, recursive queries
/// re-written, chirality relative to the output order). Directional and
/// dative bonds and chirality classes are refused.
///
/// ```
/// assert_eq!(chematic_smiles::rdkit_smarts_to_smarts("[CX3](=O)[OX2H1]").unwrap(),
///            "[C&X3](=O)[O&X2&H1]");
/// ```
pub fn rdkit_smarts_to_smarts(smarts: &str) -> Result<String, RdkitSmilesError> {
    rxn_smarts::smarts_to_smarts(smarts)
}

/// `rdChemReactions.ReactionToSmarts(rdChemReactions.ReactionFromSmarts(s))`
/// (RDKit 2026.03.1): each reactant, agent and product template written by
/// `MolToSmarts` (see [`rdkit_smarts_to_smarts`]) in input order, a
/// template with several components (`(A.B)` grouping) in parentheses.
///
/// ```
/// assert_eq!(
///     chematic_smiles::rdkit_reaction_to_smarts("[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[N:3]").unwrap(),
///     "[C:1](=[O:2])[O&H1].[N&H2:3]>>[C:1](=[O:2])[N:3]"
/// );
/// ```
pub fn rdkit_reaction_to_smarts(reaction_smarts: &str) -> Result<String, RdkitSmilesError> {
    rxn_smarts::reaction_to_smarts(reaction_smarts)
}

/// RDKit's process-wide random generator (`getRandomGenerator()`, seeded
/// with 42 at load), as the random SMILES writer draws from it.
static RDKIT_RANDOM_GENERATOR: std::sync::Mutex<Option<write::MinstdRand>> =
    std::sync::Mutex::new(None);

/// `Chem.MolToRandomSmilesVect(m, n, randomSeed, isomericSmiles,
/// kekuleSmiles, allBondsExplicit, allHsExplicit)` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1): `n` non-canonical SMILES
/// with a random root per fragment and random DFS branch order, drawn
/// from RDKit's `boost::minstd_rand` generator. `random_seed > 0` reseeds
/// the generator first; `0` continues a process-wide generator that starts
/// as RDKit's does (seed 42). `canonical` and `rooted_at_atom` of `params`
/// are ignored.
///
/// ```
/// let mol = chematic_smiles::parse("CCO").unwrap();
/// let v = chematic_smiles::rdkit_random_smiles(&mol, 2, 42, &Default::default()).unwrap();
/// assert_eq!(v.len(), 2);
/// ```
pub fn rdkit_random_smiles(
    mol: &Molecule,
    n: usize,
    random_seed: u32,
    params: &RdkitSmilesParams,
) -> Result<Vec<String>, RdkitSmilesError> {
    let p = RdkitSmilesParams {
        canonical: false,
        rooted_at_atom: None,
        ..*params
    };
    rdkit_random_smiles_with(mol, n, random_seed, &p)
}

/// `n` calls of `Chem.MolToSmiles(m, doRandom=True, **params)` after
/// seeding as in [`rdkit_random_smiles`] (with `params.canonical` the
/// fragments are sorted; a `rooted_at_atom` fixes the root of its
/// fragment).
pub fn rdkit_random_smiles_with(
    mol: &Molecule,
    n: usize,
    random_seed: u32,
    params: &RdkitSmilesParams,
) -> Result<Vec<String>, RdkitSmilesError> {
    let m = rdkit_mol_for_writing(mol)?;
    let mut guard = RDKIT_RANDOM_GENERATOR
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    if random_seed > 0 || guard.is_none() {
        *guard = Some(write::MinstdRand::new(if random_seed > 0 {
            random_seed
        } else {
            42
        }));
    }
    let rng = guard.as_mut().expect("seeded");
    write::mol_to_random_smiles(&m, n, params, rng)
}

/// `Chem.MolToCXSmarts(m)` for `m = Chem.MolFromSmiles(s)` (RDKit
/// 2026.03.1): [`rdkit_smarts`] with dative bonds written as `-` plus the
/// CXSMILES extension (radicals `^n:`, coordinate bonds `C:`). Molecules
/// with enhanced stereo groups are refused.
pub fn rdkit_cx_smarts(mol: &Molecule) -> Result<String, RdkitSmilesError> {
    if !mol.stereo_groups().is_empty() {
        return Err(RdkitSmilesError::Unsupported(
            "enhanced stereo groups".into(),
        ));
    }
    let m = rdkit_mol_for_writing(mol)?;
    let (mut res, atoms, bonds) = smarts_write::mol_to_smarts(&m, true, None, false)?;
    if !res.is_empty() {
        let ext = smarts_write::cx_extensions(&m, &atoms, &bonds);
        if !ext.is_empty() {
            res.push(' ');
            res.push_str(&ext);
        }
    }
    Ok(res)
}

/// `Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1): `MurckoDecompose` keeps
/// ring atoms, linkers between ring systems and atoms doubly bonded to
/// them, fixes the hydrogens of the atoms that lose a neighbour as RDKit
/// does (aromatic heteroatoms and aromatic carbocations get one explicit H;
/// bracket or chiral atoms get their implicit hydrogens back and lose their
/// chiral tag), then the scaffold's stereo is re-perceived by `MolToSmiles`.
pub fn rdkit_murcko_scaffold(mol: &Molecule) -> Result<String, RdkitSmilesError> {
    let m = rdkit_mol_for_writing(mol)?;
    let mut scaffold = murcko::murcko_decompose(m)?;
    stereo::legacy_stereo_perception(&mut scaffold, true, false);
    write::mol_to_smiles_owned(scaffold, &RdkitSmilesParams::default())
}

/// `rdMolHash.MolHash(m, function, useCXSmiles)` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1), every
/// `rdMolHash.HashFunction`: graph, element-graph, scaffold, tautomer,
/// mesomer and regioisomer SMILES hashes and the formula and count hashes.
pub fn rdkit_mol_hash(
    mol: &Molecule,
    function: RdkitHashFunction,
    use_cx_smiles: bool,
) -> Result<String, RdkitSmilesError> {
    if use_cx_smiles && !mol.stereo_groups().is_empty() {
        return Err(RdkitSmilesError::Unsupported(
            "enhanced stereo groups".into(),
        ));
    }
    let m = rdkit_mol_for_writing(mol)?;
    molhash::mol_hash(m, function, use_cx_smiles)
}

/// `EnumerateStereoisomers.GetStereoisomerCount(m)` (default options) for
/// `m = Chem.MolFromSmiles(s)`: `2 ** len(flippers)`. `Err` where the port
/// cannot model the molecule (including enhanced stereo groups).
pub fn rdkit_stereoisomer_count(mol: &Molecule) -> Result<u128, RdkitSmilesError> {
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
    let (atoms, bonds) = enumerate::flippers(&m);
    let n = atoms.len() + bonds.len();
    if n >= 128 {
        return Err(RdkitSmilesError::Unsupported(format!(
            "2^{n} stereoisomers"
        )));
    }
    Ok(1u128 << n)
}

/// `Chem.FindMolChiralCenters(m, force=True, includeUnassigned=...)` for
/// `m = Chem.MolFromSmiles(s)` with RDKit 2026.03's default (legacy) stereo
/// perception: `(atom index, "R" | "S" | "?")` in atom order, indices in
/// RDKit's atom numbering.
pub fn rdkit_chiral_centers(
    mol: &Molecule,
    include_unassigned: bool,
) -> Result<Vec<(usize, String)>, RdkitSmilesError> {
    let (mut m, _) = rdkit_mol_from_smiles(mol)?;
    for a in &mut m.atoms {
        a.chirality_possible = false;
    }
    stereo::legacy_stereo_perception(&mut m, true, include_unassigned);
    Ok(m.atoms
        .iter()
        .enumerate()
        .filter_map(|(i, a)| match a.cip_code {
            Some(c) => Some((i, (c as char).to_string())),
            None if include_unassigned && a.chirality_possible => Some((i, "?".to_string())),
            None => None,
        })
        .collect())
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
/// Where RDKit switches to V3000 (dative bonds, more than 999 atoms or
/// bonds, coordinates outside the V2000 fields) the V3000 CTAB is written
/// as RDKit writes it. A dummy atom without map number, isotope, charge
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
    // A block whose coordinates contradict a double bond's E/Z would be
    // read back with the other configuration: refuse it (the depiction
    // port does not reproduce RDKit's layout of some macrocycles).
    for bond in &m.bonds {
        if !matches!(bond.stereo, mol::BondStereo::E | mol::BondStereo::Z)
            || bond.stereo_atoms.len() != 2
        {
            continue;
        }
        let (b, e) = (xy[bond.begin], xy[bond.end]);
        let side = |p: [f64; 2]| (e[0] - b[0]) * (p[1] - b[1]) - (e[1] - b[1]) * (p[0] - b[0]);
        let (s0, s1) = (
            side(xy[bond.stereo_atoms[0]]),
            side(xy[bond.stereo_atoms[1]]),
        );
        let cis = s0 * s1 > 0.0;
        if s0 * s1 == 0.0 || cis != (bond.stereo == mol::BondStereo::Z) {
            return Err(RdkitSmilesError::Unsupported(
                "2D coordinates do not reproduce a double bond's E/Z configuration".into(),
            ));
        }
    }
    // An unspecified double bond RDKit does not cross because a neighbouring
    // single bond carries a direction (conjugated partial stereo) reads back
    // with the configuration its coordinates happen to show: refuse it too.
    for (b, bond) in m.bonds.iter().enumerate() {
        if bond.bt != mol::BondType::Double
            || bond.stereo != mol::BondStereo::None
            || m.num_bond_rings(b) != 0
            || !molblock2d::is_bond_potential_stereo_bond(&m, b)
        {
            continue;
        }
        let directed_nbr = [bond.begin, bond.end].iter().any(|&a| {
            m.atom_bonds[a].iter().any(|&nb| {
                nb != b && m.bonds[nb].bt == mol::BondType::Single && m.bonds[nb].dir.is_set()
            })
        });
        if directed_nbr {
            return Err(RdkitSmilesError::Unsupported(
                "an unspecified double bond next to directed bonds would read back with stereo"
                    .into(),
            ));
        }
    }
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
pub use align::{RdkitAlignment, Transform3D};

/// The sanitized RDKit molecule `Chem.MolFromSmiles` builds for `mol`
/// (no stereo perception: alignment does not use it).
fn rdkit_mol_sanitized(mol: &Molecule) -> Result<mol::Mol, RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    if has_added_hydrogens(mol) {
        sanitize::sanitize_keeping_hs(&mut m)?;
    } else {
        sanitize::remove_hs_and_sanitize(&mut m)?;
    }
    Ok(m)
}

/// `rdMolAlign.GetAlignmentTransform` / `AlignMol(prb, ref, atomMap,
/// weights, reflect, maxIters)` for two conformers of `mol` (as
/// `Chem.MolFromSmiles` numbers its atoms): the RMSD after the best fit of
/// `probe` onto `reference` and its transform. Without an atom map
/// (`(probe atom, reference atom)` pairs) the first substructure match is
/// used, as RDKit does.
pub fn rdkit_align_mol(
    mol: &Molecule,
    probe: &[[f64; 3]],
    reference: &[[f64; 3]],
    atom_map: Option<&[(usize, usize)]>,
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
) -> Result<RdkitAlignment, RdkitSmilesError> {
    let m = rdkit_mol_sanitized(mol)?;
    align::align_mol(
        &m,
        probe,
        reference,
        atom_map,
        weights,
        reflect,
        max_iterations,
    )
}

/// `rdMolAlign.GetBestRMS(prb, ref, maxMatches=max_matches,
/// symmetrizeConjugatedTerminalGroups=symmetrize, weights)` for two
/// conformers of `mol`: the smallest RMSD over the best fits of every
/// substructure match (RDKit's enumeration order, `uniquify=False`), with
/// that fit's transform and match.
pub fn rdkit_best_rms(
    mol: &Molecule,
    probe: &[[f64; 3]],
    reference: &[[f64; 3]],
    max_matches: usize,
    symmetrize: bool,
    weights: Option<&[f64]>,
) -> Result<RdkitAlignment, RdkitSmilesError> {
    let m = rdkit_mol_sanitized(mol)?;
    align::best_rms(&m, probe, reference, max_matches, symmetrize, weights)
}

/// `rdMolAlign.CalcRMS(prb, ref, maxMatches=max_matches,
/// symmetrizeConjugatedTerminalGroups=symmetrize, weights)`: the smallest
/// RMSD over all matches without moving the probe.
pub fn rdkit_calc_rms(
    mol: &Molecule,
    probe: &[[f64; 3]],
    reference: &[[f64; 3]],
    max_matches: usize,
    symmetrize: bool,
    weights: Option<&[f64]>,
) -> Result<f64, RdkitSmilesError> {
    let m = rdkit_mol_sanitized(mol)?;
    align::calc_rms(&m, probe, reference, max_matches, symmetrize, weights)
}

/// `RDNumeric::Alignments::AlignPoints(refPoints, probePoints, weights,
/// reflect, maxIterations)`: the sum of squared residuals of the best fit
/// and its transform.
pub fn rdkit_align_points(
    reference: &[[f64; 3]],
    probe: &[[f64; 3]],
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
) -> Result<(f64, Transform3D), RdkitSmilesError> {
    align::align_points(reference, probe, weights, reflect, max_iterations)
}

/// A molecule read by [`rdkit_mol_from_pdb_block`].
#[derive(Clone)]
pub struct RdkitPdbMolecule {
    /// The molecule, in RDKit's atom order.
    pub molecule: Molecule,
    /// One position per atom (the first model).
    pub coords: Vec<[f64; 3]>,
    /// `Chem.MolToSmiles` of RDKit's molecule.
    pub smiles: String,
}

/// `Chem.MolFromPDBBlock(text, sanitize, removeHs, flavor,
/// proximityBonding)` (RDKit 2026.03.1): `None` where RDKit returns no
/// molecule. `Err` where RDKit would fail (or the port cannot model the
/// result, such as zero-order bonds).
pub fn rdkit_mol_from_pdb_block(
    text: &str,
    sanitize: bool,
    remove_hs: bool,
    flavor: u32,
    proximity_bonding: bool,
) -> Result<Option<RdkitPdbMolecule>, RdkitSmilesError> {
    let Some(pdb) =
        pdb_read::mol_from_pdb_block(text, sanitize, remove_hs, flavor, proximity_bonding)?
    else {
        return Ok(None);
    };
    let molecule = pdb_read::to_chematic(&pdb.mol, sanitize)?;
    let mut m = pdb.mol;
    let smiles = if sanitize {
        // `MolToSmiles` runs `assignStereochemistry(cleanIt=true)` first.
        stereo::legacy_stereo_perception(&mut m, true, false);
        write::mol_to_smiles_owned(m, &RdkitSmilesParams::default())?
    } else {
        String::new()
    };
    Ok(Some(RdkitPdbMolecule {
        molecule,
        coords: pdb.coords,
        smiles,
    }))
}

/// `Chem.MolFromXYZBlock(text)` (RDKit 2026.03.1): the atoms (no bonds)
/// and their positions. `Err` with
/// [`RdkitSmilesError::Sanitization`] where RDKit's parser fails.
pub fn rdkit_mol_from_xyz_block(text: &str) -> Result<(Molecule, Vec<[f64; 3]>), RdkitSmilesError> {
    let (m, coords) = xyz_read::mol_from_xyz_block(text)?;
    Ok((pdb_read::to_chematic(&m, false)?, coords))
}

/// A molecule read by [`rdkit_mol_from_mol2_block`].
#[derive(Clone)]
pub struct RdkitMol2Molecule {
    /// The molecule, in RDKit's atom order.
    pub molecule: Molecule,
    /// One position per atom.
    pub coords: Vec<[f64; 3]>,
    /// `Chem.MolToSmiles` of RDKit's molecule (empty when not sanitized).
    pub smiles: String,
}

/// `Chem.MolFromMol2Block(text, sanitize, removeHs, cleanupSubstructures)`
/// (RDKit 2026.03.1): Tripos atom types, Corina-style substructure cleanup
/// and formal-charge guessing (or `UNITY_ATOM_ATTR` charges), chirality and
/// double-bond stereo from the 3D coordinates. `Err` with
/// [`RdkitSmilesError::Sanitization`] where RDKit returns no molecule;
/// [`RdkitSmilesError::Unsupported`] for features the port does not model
/// (Tripos query atoms, `du`/`un` bonds).
pub fn rdkit_mol_from_mol2_block(
    text: &str,
    sanitize: bool,
    remove_hs: bool,
    cleanup_substructures: bool,
) -> Result<RdkitMol2Molecule, RdkitSmilesError> {
    let (m, coords) =
        mol2_read::mol_from_mol2_block(text, sanitize, remove_hs, cleanup_substructures)?;
    let molecule = pdb_read::to_chematic(&m, sanitize)?;
    let smiles = if sanitize {
        // Stereochemistry is already assigned (`_StereochemDone`).
        write::mol_to_smiles_owned(m, &RdkitSmilesParams::default())?
    } else {
        String::new()
    };
    Ok(RdkitMol2Molecule {
        molecule,
        coords,
        smiles,
    })
}

/// [`rdkit_mol_from_smiles`]'s molecule for writers, which read neither
/// `chirality_possible` nor the CIP ranks.
/// The result of [`rdkit_enumerate_tautomers`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitTautomerEnumeration {
    /// `Chem.MolToSmiles(t)` of every tautomer, in the result's order (its
    /// map keys' order).
    pub smiles: Vec<String>,
    /// The result's map keys (`TautomerEnumeratorResult.smiles`).
    pub keys: Vec<String>,
    /// `TautomerEnumeratorResult.status`.
    pub status: RdkitTautomerStatus,
}

/// `Chem.MolToSmiles(rdMolStandardize.TautomerEnumerator().Canonicalize(m))`
/// for `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1, default
/// `CleanupParameters`): the port of `TautomerEnumerator::canonicalize`
/// (transforms, enumeration order and limits, `scoreTautomer`, ties broken
/// by the smaller canonical SMILES).
///
/// ```
/// let mol = chematic_smiles::parse("Oc1ccccn1").unwrap();
/// assert_eq!(chematic_smiles::rdkit_canonical_tautomer(&mol).unwrap(), "O=c1cccc[nH]1");
/// ```
pub fn rdkit_canonical_tautomer(mol: &Molecule) -> Result<String, RdkitSmilesError> {
    let (m, _) = rdkit_mol_from_smiles(mol)?;
    let (canon, _) = tautomer::canonicalize(&m)?;
    write::mol_to_smiles_owned(canon, &RdkitSmilesParams::default())
}

/// `rdMolStandardize.TautomerEnumerator().Enumerate(m)` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1, default
/// `CleanupParameters`): the tautomers' SMILES and the result status.
pub fn rdkit_enumerate_tautomers(
    mol: &Molecule,
) -> Result<RdkitTautomerEnumeration, RdkitSmilesError> {
    let (m, _) = rdkit_mol_from_smiles(mol)?;
    let res = tautomer::enumerate(&m, &tautomer::Settings::default())?;
    let mut smiles = Vec::with_capacity(res.tautomers.len());
    for t in res.tautomers.values() {
        smiles.push(write::mol_to_smiles(&t.mol, &RdkitSmilesParams::default())?);
    }
    Ok(RdkitTautomerEnumeration {
        smiles,
        keys: res.tautomers.keys().cloned().collect(),
        status: res.status,
    })
}

/// `rdMolStandardize.TautomerEnumerator.ScoreTautomer(m)` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1): ring, substructure and
/// hetero-H terms of `TautomerScoringFunctions::scoreTautomer`.
pub fn rdkit_tautomer_score(mol: &Molecule) -> Result<i32, RdkitSmilesError> {
    let (m, _) = rdkit_mol_from_smiles(mol)?;
    Ok(tautomer::score_tautomer(&m))
}

fn rdkit_mol_for_writing(mol: &Molecule) -> Result<mol::Mol, RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    if has_added_hydrogens(mol) {
        sanitize::sanitize_keeping_hs(&mut m)?;
    } else {
        sanitize::remove_hs_and_sanitize(&mut m)?;
    }
    stereo::legacy_stereo_perception_unflagged(&mut m);
    Ok(m)
}

/// RDKit's sanitized model of `mol` (`Chem.MolFromSmiles`) projected onto
/// chematic's graph: what `removeHs` + `sanitizeMol` (cleanup, Kekulize,
/// aromaticity) leave on each atom and bond, by chematic index. Atoms RDKit
/// removes (explicit hydrogens) and their bonds are `None`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitSanitizedModel {
    /// Per chematic atom: `(is_aromatic, formal_charge, total H count
    /// counting removed hydrogen neighbours as graph atoms)`.
    pub atoms: Vec<Option<(bool, i8, u8)>>,
    /// Per chematic atom: number of hydrogen neighbours RDKit removed.
    pub removed_h_neighbors: Vec<u8>,
    /// Per chematic bond: RDKit's bond as a chematic order (`Aromatic` for
    /// aromatic bonds; `Single`, `Double`, `Triple`, `Quadruple`, `Dative`
    /// otherwise).
    pub bonds: Vec<Option<chematic_core::BondOrder>>,
}

/// [`RdkitSanitizedModel`] for `mol` (RDKit 2026.03.1). `Err` where the
/// port cannot model the molecule or RDKit's sanitization rejects it.
pub fn rdkit_sanitized_model(mol: &Molecule) -> Result<RdkitSanitizedModel, RdkitSmilesError> {
    use chematic_core::BondOrder;
    let order = mol.rdkit_bond_order();
    let mut m = parse::from_chematic_ordered(mol, &order)?;
    let n = m.atoms.len();
    // Bonds before H removal, by RDKit index: (begin, end, chematic bond).
    let bonds_before: Vec<(usize, usize, usize)> = m
        .bonds
        .iter()
        .zip(order.iter())
        .map(|(b, cb)| (b.begin, b.end, cb.0 as usize))
        .collect();
    let kept: Vec<usize> = if has_added_hydrogens(mol) {
        sanitize::sanitize_keeping_hs(&mut m)?;
        (0..n).collect()
    } else {
        sanitize::remove_hs(&mut m, true)?
    };
    let mut new_idx = vec![usize::MAX; n];
    for (k, &old) in kept.iter().enumerate() {
        new_idx[old] = k;
    }
    let mut atoms = vec![None; n];
    let mut removed_h_neighbors = vec![0u8; n];
    for (b0, e0, _) in &bonds_before {
        if new_idx[*b0] == usize::MAX && new_idx[*e0] != usize::MAX {
            removed_h_neighbors[*e0] += 1;
        } else if new_idx[*e0] == usize::MAX && new_idx[*b0] != usize::MAX {
            removed_h_neighbors[*b0] += 1;
        }
    }
    for (old, &k) in new_idx.iter().enumerate() {
        if k == usize::MAX {
            continue;
        }
        let a = &m.atoms[k];
        atoms[old] = Some((a.aromatic, a.charge as i8, m.total_num_hs(k).min(255) as u8));
    }
    let mut bonds = vec![None; mol.bond_count()];
    let mut next = 0usize;
    for (b0, e0, cb) in bonds_before {
        if new_idx[b0] == usize::MAX || new_idx[e0] == usize::MAX {
            continue;
        }
        let bond = &m.bonds[next];
        next += 1;
        // RDKit flags a ring triple (or double) bond it could not make
        // aromatic-typed as aromatic but keeps its type (`c1ccccc#1`); the
        // type is what chematic's bond order can carry.
        bonds[cb] = Some(
            if bond.aromatic && !matches!(bond.bt, mol::BondType::Triple) {
                BondOrder::Aromatic
            } else {
                match bond.bt {
                    mol::BondType::Single => BondOrder::Single,
                    mol::BondType::Double => BondOrder::Double,
                    mol::BondType::Triple => BondOrder::Triple,
                    mol::BondType::Quadruple => BondOrder::Quadruple,
                    mol::BondType::Aromatic => BondOrder::Aromatic,
                    mol::BondType::Dative => BondOrder::Dative,
                }
            },
        );
    }
    Ok(RdkitSanitizedModel {
        atoms,
        removed_h_neighbors,
        bonds,
    })
}

/// `Chem.MolFromSmiles(s)` as a chematic molecule when RDKit's `removeHs`
/// takes hydrogen graph atoms of `mol` away (`None` otherwise, or where the
/// port cannot model `mol`): RDKit's atoms and bonds in RDKit's order, each
/// heavy atom with its total hydrogen count. RDKit-compatible fingerprints
/// read this molecule so that `[H]C([H])([H])[H]` is methane, as in RDKit.
pub fn rdkit_hydrogen_suppressed(mol: &Molecule) -> Option<Molecule> {
    if !mol
        .atoms()
        .any(|(_, a)| !a.wildcard && a.element == chematic_core::Element::H)
    {
        return None;
    }
    let (m, _) = rdkit_mol_from_smiles(mol).ok()?;
    if m.atoms.len() == mol.atom_count() {
        return None;
    }
    pdb_read::to_chematic(&m, true).ok()
}

/// RDKit's hybridization of every atom of `mol` (`Chem.MolFromSmiles`):
/// 1 SP, 2 SP2, 3 SP3, 0 otherwise (S, SP2D, SP3D, SP3D2, unspecified).
/// `None` where RDKit removes atoms (explicit hydrogens) or the port cannot
/// model `mol`.
pub fn rdkit_hybridizations(mol: &Molecule) -> Option<Vec<u8>> {
    let (m, _) = rdkit_mol_from_smiles(mol).ok()?;
    if m.atoms.len() != mol.atom_count() {
        return None;
    }
    Some(
        m.atoms
            .iter()
            .map(|a| match a.hybrid {
                mol::Hybridization::Sp => 1,
                mol::Hybridization::Sp2 => 2,
                mol::Hybridization::Sp3 => 3,
                _ => 0,
            })
            .collect(),
    )
}

/// RDKit's `MolOps::addHs` on a sanitized molecule (no coordinates): every
/// atom's hydrogens become graph atoms appended in atom order, the atom
/// keeping none.
fn add_hs_graph(m: &mut mol::Mol) -> Result<(), RdkitSmilesError> {
    let n = m.atoms.len();
    for a in 0..n {
        let hs = m.total_num_hs(a);
        for _ in 0..hs {
            let mut h = mol::Atom::new(1);
            h.no_implicit = false;
            let hi = m.add_atom(h);
            m.add_bond(mol::Bond::new(a, hi, mol::BondType::Single));
        }
        m.atoms[a].num_explicit_hs = 0;
        m.atoms[a].no_implicit = true;
    }
    if let Some(ri) = &m.rings {
        let rings = ri.atom_rings.clone();
        m.set_rings(rings);
    }
    m.update_property_cache(false)
}

/// For RDKit's `Chem.AddHs(Chem.MolFromSmiles(s))` molecule: per atom, the
/// index of the first of `patterns` (SMARTS) with a match starting at that
/// atom (`SubstructMatch(..., uniquify=false)`), as RDKit's Crippen atom
/// typing assigns types; and the number of heavy (original) atoms, which
/// come first. `None` where the port cannot model `mol`.
pub fn rdkit_addhs_first_pattern(
    mol: &Molecule,
    patterns: &[&str],
) -> Option<(Vec<Option<usize>>, usize)> {
    thread_local! {
        static CACHE: std::cell::RefCell<std::collections::HashMap<String, Option<smarts_match::Query>>> =
            std::cell::RefCell::new(std::collections::HashMap::new());
    }
    let (mut m, _) = rdkit_mol_from_smiles(mol).ok()?;
    let heavy = m.atoms.len();
    add_hs_graph(&mut m).ok()?;
    let mut types: Vec<Option<usize>> = vec![None; m.atoms.len()];
    let mut left = types.len();
    for (pi, pat) in patterns.iter().enumerate() {
        if left == 0 {
            break;
        }
        let hits = CACHE.with(|c| {
            let mut c = c.borrow_mut();
            let q = c
                .entry(pat.to_string())
                .or_insert_with(|| std::panic::catch_unwind(|| smarts_match::parse(pat)).ok());
            q.as_ref().map(|q| smarts_match::first_atoms(q, &m))
        });
        let Some(hits) = hits else { continue };
        for (a, hit) in hits.into_iter().enumerate() {
            if hit && types[a].is_none() {
                types[a] = Some(pi);
                left -= 1;
            }
        }
    }
    Some((types, heavy))
}

/// RDKit's `getCrippenAtomContribs` on `Chem.AddHs(Chem.MolFromSmiles(s))`
/// (RDKit 2026.03.1 `Crippen.cpp` parameters): per atom of that molecule
/// (heavy atoms first, then the added hydrogens) its `(logP, MR)`
/// contribution, and the number of heavy atoms. `None` where the port
/// cannot model `mol`.
pub fn rdkit_crippen_contribs(mol: &Molecule) -> Option<(Vec<(f64, f64)>, usize)> {
    crippen_contribs(mol, true)
}

/// `getCrippenAtomContribs(Chem.MolFromSmiles(s))` (no added hydrogens), as
/// RDKit's SlogP/SMR VSA descriptors call it: per atom of RDKit's molecule.
pub fn rdkit_crippen_contribs_no_hs(mol: &Molecule) -> Option<Vec<(f64, f64)>> {
    crippen_contribs(mol, false).map(|(c, _)| c)
}

fn crippen_contribs(mol: &Molecule, add_hs: bool) -> Option<(Vec<(f64, f64)>, usize)> {
    thread_local! {
        static QUERIES: Vec<Option<smarts_match::Query>> = crippen_params::CRIPPEN_PARAMS
            .iter()
            .map(|p| std::panic::catch_unwind(|| smarts_match::parse(p.1)).ok())
            .collect();
    }
    let (mut m, _) = rdkit_mol_from_smiles(mol).ok()?;
    let heavy = m.atoms.len();
    if add_hs {
        add_hs_graph(&mut m).ok()?;
    }
    let mut contribs = vec![(0.0, 0.0); m.atoms.len()];
    let mut needed = vec![true; m.atoms.len()];
    let mut left = needed.len();
    QUERIES.with(|qs| {
        for (pi, q) in qs.iter().enumerate() {
            if left == 0 {
                break;
            }
            let Some(q) = q else { continue };
            for (a, hit) in smarts_match::first_atoms(q, &m).into_iter().enumerate() {
                if hit && needed[a] {
                    needed[a] = false;
                    left -= 1;
                    let p = crippen_params::CRIPPEN_PARAMS[pi];
                    contribs[a] = (p.2, p.3);
                }
            }
        }
    });
    Some((contribs, heavy))
}

/// `Crippen.MolLogP(m)` and `Crippen.MolMR(m)` for `m = MolFromSmiles(s)`
/// (contributions summed in atom order).
pub fn rdkit_crippen_logp_mr(mol: &Molecule) -> Option<(f64, f64)> {
    let (c, _) = rdkit_crippen_contribs(mol)?;
    let mut logp = 0.0f64;
    let mut mr = 0.0f64;
    for (l, _) in &c {
        logp += l;
    }
    for (_, r) in &c {
        mr += r;
    }
    Some((logp, mr))
}

/// `Chem.MolFromSmiles(s)` as a chematic molecule (RDKit's atoms, bonds,
/// aromaticity, hydrogen counts and kept stereo, in RDKit's order), for a
/// molecule read from `s` without added hydrogens; `None` where the port
/// cannot model it.
pub fn rdkit_parsed_molecule(mol: &Molecule) -> Option<Molecule> {
    if has_added_hydrogens(mol) {
        return None;
    }
    let (m, _) = rdkit_mol_from_smiles(mol).ok()?;
    pdb_read::to_chematic(&m, true).ok()
}

/// `Chem.AddHs(Chem.MolFromSmiles(s))` as a chematic molecule built by the
/// port (heavy atoms in RDKit's order, then every hydrogen in atom order),
/// for a molecule read from `s`; `None` where the port cannot model it.
pub fn rdkit_added_hs_molecule(mol: &Molecule) -> Option<Molecule> {
    if has_added_hydrogens(mol) {
        return None;
    }
    let (mut m, _) = rdkit_mol_from_smiles(mol).ok()?;
    add_hs_graph(&mut m).ok()?;
    pdb_read::to_chematic(&m, true).ok()
}

/// [`rdkit_hydrogen_suppressed`] with, per atom of the returned molecule,
/// the index of the atom of `mol` it is (RDKit's `removeHs` keeps the
/// other atoms in order).
pub fn rdkit_hydrogen_suppressed_with_map(mol: &Molecule) -> Option<(Molecule, Vec<usize>)> {
    if has_added_hydrogens(mol)
        || !mol
            .atoms()
            .any(|(_, a)| !a.wildcard && a.element == chematic_core::Element::H)
    {
        return None;
    }
    let mut m = parse::from_chematic(mol).ok()?;
    let kept = sanitize::remove_hs(&mut m, true).ok()?;
    if kept.len() == mol.atom_count() {
        return None;
    }
    stereo::legacy_stereo_perception(&mut m, true, true);
    let out = pdb_read::to_chematic(&m, true).ok()?;
    (out.atom_count() == kept.len()).then_some((out, kept))
}

thread_local! {
    static IN_RDKIT_MODEL: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Runs `f` on RDKit's own molecule for `mol` (`Chem.MolFromSmiles(s)`:
/// hydrogen graph atoms removed, RDKit's aromaticity, bond orders, charges
/// and per-atom hydrogen counts, as [`rdkit_hydrogen_suppressed`] builds it)
/// when chematic's reading of `mol` may differ from RDKit's
/// ([`chematic_perception::rdkit_model_may_disagree`], or hydrogen graph
/// atoms), and on `mol` itself otherwise. Nested calls run on the molecule
/// they are given.
pub fn with_rdkit_model_molecule<R>(mol: &Molecule, f: impl FnOnce(&Molecule) -> R) -> R {
    if IN_RDKIT_MODEL.with(|c| c.get()) {
        return f(mol);
    }
    let needed = chematic_perception::rdkit_model_may_disagree(mol)
        || mol
            .atoms()
            .any(|(_, a)| !a.wildcard && a.element == chematic_core::Element::H);
    let model = if needed {
        rdkit_mol_from_smiles(mol)
            .ok()
            .and_then(|(m, _)| pdb_read::to_chematic(&m, true).ok())
    } else {
        None
    };
    match model {
        Some(m) => {
            struct Reset;
            impl Drop for Reset {
                fn drop(&mut self) {
                    IN_RDKIT_MODEL.with(|c| c.set(false));
                }
            }
            IN_RDKIT_MODEL.with(|c| c.set(true));
            let _reset = Reset;
            f(&m)
        }
        None => f(mol),
    }
}

/// RDKit's `numPiElectrons` of every atom of `mol` (`Chem.MolFromSmiles`):
/// 1 for an aromatic atom, 0 for an SP3 atom, otherwise the explicit
/// valence beyond one per physical bond. `None` where RDKit removes atoms
/// or the port cannot model `mol`.
pub fn rdkit_num_pi_electrons(mol: &Molecule) -> Option<Vec<u32>> {
    let (m, _) = rdkit_mol_from_smiles(mol).ok()?;
    if m.atoms.len() != mol.atom_count() {
        return None;
    }
    Some(
        (0..m.atoms.len())
            .map(|a| {
                let atom = &m.atoms[a];
                if atom.aromatic {
                    1
                } else if atom.hybrid != mol::Hybridization::Sp3 {
                    let val = atom.explicit_valence.max(0) as u32;
                    let physical = atom.num_explicit_hs
                        + m.atom_bonds[a]
                            .iter()
                            .filter(|&&b| m.bonds[b].valence_contrib(a) != 0.0)
                            .count() as u32;
                    val.saturating_sub(physical)
                } else {
                    0
                }
            })
            .collect(),
    )
}

/// Registers [`rdkit_model_correct_view`] as chematic-perception's
/// RDKit-model hook, so the shared RDKit-parity aromatic view (used by the
/// RDKit-compatible fingerprints, descriptors and SMARTS matching) follows
/// this crate's port of RDKit's sanitization on molecules where the parity
/// view may disagree with RDKit. Called by every SMILES parse; idempotent.
pub fn register_rdkit_model_hook() {
    chematic_perception::set_rdkit_model_hook(rdkit_model_correct_view);
}

/// The bond-order class RDKit's model distinguishes.
fn order_class(order: chematic_core::BondOrder) -> u8 {
    use chematic_core::BondOrder;
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 12,
        BondOrder::Dative => 17,
        _ => 0,
    }
}

/// `view` (an RDKit-parity view of `mol`, on `mol`'s graph) with the
/// aromaticity, bond orders, charges and hydrogen counts of
/// [`rdkit_sanitized_model`] where they disagree; `None` when they agree or
/// the port cannot model `mol`.
pub fn rdkit_model_correct_view(mol: &Molecule, view: &Molecule) -> Option<Molecule> {
    use chematic_core::{AtomIdx, BondIdx};
    if view.atom_count() != mol.atom_count() || view.bond_count() != mol.bond_count() {
        return None;
    }
    let model = rdkit_sanitized_model(mol).ok()?;
    let atom_differs = |i: usize| -> bool {
        match model.atoms[i] {
            Some((arom, charge, _)) => {
                let a = view.atom(AtomIdx(i as u32));
                a.aromatic != arom || a.charge != charge
            }
            None => false,
        }
    };
    let bond_differs = |b: usize| -> bool {
        match model.bonds[b] {
            Some(order) => order_class(view.bond(BondIdx(b as u32)).order) != order_class(order),
            None => false,
        }
    };
    if !(0..mol.atom_count()).any(atom_differs) && !(0..mol.bond_count()).any(bond_differs) {
        return None;
    }
    let mut out = view.clone();
    let mut touched = vec![false; mol.atom_count()];
    for b in 0..mol.bond_count() {
        if bond_differs(b) {
            let order = model.bonds[b].expect("modelled bond");
            let bond = view.bond(BondIdx(b as u32));
            touched[bond.atom1.0 as usize] = true;
            touched[bond.atom2.0 as usize] = true;
            out.set_bond_order(BondIdx(b as u32), order);
        }
    }
    for i in 0..mol.atom_count() {
        let Some((arom, charge, hs)) = model.atoms[i] else {
            continue;
        };
        if atom_differs(i) || touched[i] {
            let idx = AtomIdx(i as u32);
            out.set_atom_aromatic(idx, arom);
            out.set_charge(idx, charge);
            if !view.atom(idx).wildcard {
                out.set_hydrogen_count(idx, Some(hs.saturating_sub(model.removed_h_neighbors[i])));
            }
        }
    }
    Some(out)
}

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
mod profile_tests;
#[cfg(test)]
mod tautomer_tests;
#[cfg(test)]
mod tests;
