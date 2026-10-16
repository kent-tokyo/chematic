//! Stereo-checked rooted serialization, motivated by RDKit #9646.
use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule, RenumberError, STEREO_H_SENTINEL};

/// Maximum graph size for the recursive rooted writer's checked profile.
pub const ROOTED_SMILES_MAX_ATOMS: usize = 256;

/// Refusal to serialize a rooted molecule without a verified stereo frame.
#[derive(Debug)]
#[non_exhaustive]
pub enum RootedSmilesError {
    InvalidRoot { root: usize, atoms: usize },
    ResourceLimit { atoms: usize, limit: usize },
    MissingStereoFrame { atom: usize },
    InvalidStereoFrame { atom: usize },
    UnsupportedStereo { atom: Option<usize> },
    UnsupportedBond { bond: usize, order: BondOrder },
    InvalidMetadata(RenumberError),
}

impl core::fmt::Display for RootedSmilesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidRoot { root, atoms } => {
                write!(f, "root atom {root} is outside {atoms} atoms")
            }
            Self::ResourceLimit { atoms, limit } => {
                write!(f, "rooted SMILES atom limit: {atoms} > {limit}")
            }
            Self::MissingStereoFrame { atom } => {
                write!(f, "atom {atom} has no recorded stereo frame")
            }
            Self::InvalidStereoFrame { atom } => {
                write!(f, "atom {atom} has an invalid tetrahedral stereo frame")
            }
            Self::UnsupportedStereo { atom } => write!(
                f,
                "rooted SMILES supports ordinary tetrahedral stereo only (atom {atom:?}); enhanced groups are unsupported"
            ),
            Self::UnsupportedBond { bond, order } => write!(
                f,
                "rooted SMILES bond {bond} has unsupported order {order:?}"
            ),
            Self::InvalidMetadata(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for RootedSmilesError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidMetadata(error) => Some(error),
            _ => None,
        }
    }
}

/// Write non-canonical SMILES starting at `root`, with original atom indices.
///
/// Returns `(text, order)`, where `order[parsed_index]` identifies the original
/// graph atom. Disconnected components follow the root component. No chemical
/// perception or source mutation occurs. Four-ligand tetrahedral frames must
/// be recorded explicitly: three graph neighbors plus one inferred/bracket H,
/// or four graph neighbors without H. For a constructed graph, the caller
/// specifies the frame to which its `@`/`@@` tag refers; adjacency is never
/// guessed. Implicit H need not be converted to an explicit bracket-H count.
///
/// Rejects invalid roots, graphs above [`ROOTED_SMILES_MAX_ATOMS`], missing or
/// malformed frames, axial/non-tetrahedral/lone-pair stereo, enhanced groups,
/// and zero/query bond orders. Other SMILES chemistry uses the existing writer.
/// This projects SMILES chemistry, not arbitrary source/depiction metadata.
/// Existing native and RDKit-compatibility writers retain their contracts.
pub fn write_rooted_checked(
    mol: &Molecule,
    root: usize,
) -> Result<(String, Vec<AtomIdx>), RootedSmilesError> {
    let n = mol.atom_count();
    if n > ROOTED_SMILES_MAX_ATOMS {
        return Err(RootedSmilesError::ResourceLimit {
            atoms: n,
            limit: ROOTED_SMILES_MAX_ATOMS,
        });
    }
    if root >= n {
        return Err(RootedSmilesError::InvalidRoot { root, atoms: n });
    }
    validate_stereo(mol)?;
    validate_bonds(mol)?;
    let permutation: Vec<_> = std::iter::once(AtomIdx(root as u32))
        .chain((0..n).filter(|&i| i != root).map(|i| AtomIdx(i as u32)))
        .collect();
    let reordered = mol
        .renumber_atoms(&permutation)
        .map_err(RootedSmilesError::InvalidMetadata)?;
    let (text, order) = crate::write_with_atom_order(&reordered);
    Ok((
        text,
        order
            .into_iter()
            .map(|i| permutation[i.0 as usize])
            .collect(),
    ))
}

fn validate_stereo(mol: &Molecule) -> Result<(), RootedSmilesError> {
    if !mol.stereo_groups().is_empty() {
        return Err(RootedSmilesError::UnsupportedStereo { atom: None });
    }
    for (idx, atom) in mol.atoms() {
        if atom.chirality == Chirality::None {
            continue;
        }
        if !atom.chirality.is_tetrahedral() || mol.is_allene_stereo(idx) {
            return Err(RootedSmilesError::UnsupportedStereo {
                atom: Some(idx.0 as usize),
            });
        }
        validate_frame(mol, idx)?;
    }
    Ok(())
}

fn validate_frame(mol: &Molecule, idx: AtomIdx) -> Result<(), RootedSmilesError> {
    let atom = idx.0 as usize;
    let order = mol
        .stereo_neighbor_order(idx)
        .ok_or(RootedSmilesError::MissingStereoFrame { atom })?;
    if order.len() != 4 {
        return Err(RootedSmilesError::InvalidStereoFrame { atom });
    }
    let h = chematic_core::implicit_hcount(mol, idx);
    let mut expected: Vec<_> = mol.neighbors(idx).map(|(nb, _)| nb.0).collect();
    if expected.len() == 3 && h == 1 {
        expected.push(STEREO_H_SENTINEL);
    } else if expected.len() != 4 || h != 0 {
        return Err(RootedSmilesError::InvalidStereoFrame { atom });
    }
    expected.sort_unstable();
    let mut actual = order.to_vec();
    actual.sort_unstable();
    if actual != expected {
        return Err(RootedSmilesError::InvalidStereoFrame { atom });
    }
    Ok(())
}

fn validate_bonds(mol: &Molecule) -> Result<(), RootedSmilesError> {
    for (idx, bond) in mol.bonds() {
        if matches!(
            bond.order,
            BondOrder::Zero
                | BondOrder::QueryAny
                | BondOrder::QuerySingleOrDouble
                | BondOrder::QuerySingleOrAromatic
                | BondOrder::QueryDoubleOrAromatic
        ) {
            return Err(RootedSmilesError::UnsupportedBond {
                bond: idx.0 as usize,
                order: bond.order,
            });
        }
    }
    Ok(())
}
