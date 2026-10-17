//! Checked selection from supplied candidates; no tautomer enumeration.

use chematic_core::Molecule;

use super::{RdkitSmilesError, RdkitSmilesParams, rdkit_mol_from_smiles, tautomer, write};

/// Hard input-count bound, including duplicate candidates.
pub const RDKIT_TAUTOMER_PICK_MAX_CANDIDATES: usize = 4096;

/// Selection in the caller's original iteration order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitTautomerPick {
    /// Zero-based position of the winning input. Identical ties keep the first.
    pub index: usize,
    /// Existing RDKit-port default tautomer score.
    pub score: i32,
    /// RDKit-port canonical SMILES used to break equal-score ties.
    pub smiles: String,
}

/// Invalid selection input or a candidate outside the existing RDKit port.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitTautomerPickError {
    InvalidLimit,
    Empty,
    TooManyCandidates {
        limit: usize,
    },
    CandidateTooLarge {
        index: usize,
    },
    Candidate {
        index: usize,
        error: RdkitSmilesError,
    },
}

impl std::fmt::Display for RdkitTautomerPickError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidLimit => write!(f, "max_candidates must be in 1..=4096"),
            Self::Empty => write!(f, "candidates must contain at least one molecule"),
            Self::TooManyCandidates { limit } => {
                write!(f, "candidates exceed max_candidates ({limit})")
            }
            Self::CandidateTooLarge { index } => {
                write!(f, "candidate {index} exceeds 1024 atoms or 4096 bonds")
            }
            Self::Candidate { index, error } => write!(f, "candidate {index}: {error}"),
        }
    }
}

impl std::error::Error for RdkitTautomerPickError {}

/// Pick the highest-scoring supplied molecule, breaking ties by the smallest
/// RDKit-port canonical SMILES, then by first occurrence. Returns its input
/// index, score and key; it does not transform or copy the winning molecule.
///
/// Accepts up to `max_candidates` inputs (1..=4096), each with at most 1024
/// atoms and 4096 bonds. Duplicates count toward the limit. At most one extra
/// input is requested to detect overflow; overflow never returns a partial
/// winner. Empty collections, including empty iterators, return a typed error.
///
/// Uses the existing RDKit 2026.03.1 port's default scoring and lexical tie
/// policy. This is not general RDKit `PickCanonical` emulation: no custom
/// scorer, newer stereo-count tie policy, stereo cleanup of the selected
/// source, or validation that the inputs belong to one tautomer family.
pub fn rdkit_pick_tautomer<'a>(
    candidates: impl IntoIterator<Item = &'a Molecule>,
    max_candidates: usize,
) -> Result<RdkitTautomerPick, RdkitTautomerPickError> {
    if !(1..=RDKIT_TAUTOMER_PICK_MAX_CANDIDATES).contains(&max_candidates) {
        return Err(RdkitTautomerPickError::InvalidLimit);
    }
    let mut best: Option<RdkitTautomerPick> = None;
    for (index, candidate) in candidates.into_iter().enumerate() {
        if index == max_candidates {
            return Err(RdkitTautomerPickError::TooManyCandidates {
                limit: max_candidates,
            });
        }
        if candidate.atom_count() > 1024 || candidate.bond_count() > 4096 {
            return Err(RdkitTautomerPickError::CandidateTooLarge { index });
        }
        let (model, _) = rdkit_mol_from_smiles(candidate)
            .map_err(|error| RdkitTautomerPickError::Candidate { index, error })?;
        let score = tautomer::score_tautomer(&model);
        let smiles = write::mol_to_smiles_owned(model, &RdkitSmilesParams::default())
            .map_err(|error| RdkitTautomerPickError::Candidate { index, error })?;
        if best
            .as_ref()
            .is_none_or(|b| score > b.score || (score == b.score && smiles < b.smiles))
        {
            best = Some(RdkitTautomerPick {
                index,
                score,
                smiles,
            });
        }
    }
    best.ok_or(RdkitTautomerPickError::Empty)
}
