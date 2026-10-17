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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse;

    #[test]
    fn checked_selection_and_limits_cover_the_public_contract() {
        let enol = parse("C=C(C)O").unwrap();
        let keto = parse("CC(C)=O").unwrap();
        let pick = rdkit_pick_tautomer([&enol, &keto], 2).unwrap();
        assert_eq!(
            (pick.index, pick.score, pick.smiles.as_str()),
            (1, 5, "CC(C)=O")
        );

        assert_eq!(
            rdkit_pick_tautomer([], 1),
            Err(RdkitTautomerPickError::Empty)
        );
        assert_eq!(
            rdkit_pick_tautomer([], 0),
            Err(RdkitTautomerPickError::InvalidLimit)
        );
        assert_eq!(
            rdkit_pick_tautomer([&enol, &keto], 1),
            Err(RdkitTautomerPickError::TooManyCandidates { limit: 1 })
        );

        let large = parse(&"C".repeat(1025)).unwrap();
        assert_eq!(
            rdkit_pick_tautomer([&large], 1),
            Err(RdkitTautomerPickError::CandidateTooLarge { index: 0 })
        );
        let invalid = parse("[CH5]").unwrap();
        assert!(matches!(
            rdkit_pick_tautomer([&invalid], 1),
            Err(RdkitTautomerPickError::Candidate { index: 0, .. })
        ));
    }

    #[test]
    fn typed_errors_have_actionable_messages() {
        assert_eq!(
            RdkitTautomerPickError::InvalidLimit.to_string(),
            "max_candidates must be in 1..=4096"
        );
        assert_eq!(
            RdkitTautomerPickError::Empty.to_string(),
            "candidates must contain at least one molecule"
        );
        assert_eq!(
            RdkitTautomerPickError::TooManyCandidates { limit: 7 }.to_string(),
            "candidates exceed max_candidates (7)"
        );
        assert_eq!(
            RdkitTautomerPickError::CandidateTooLarge { index: 3 }.to_string(),
            "candidate 3 exceeds 1024 atoms or 4096 bonds"
        );
    }
}
