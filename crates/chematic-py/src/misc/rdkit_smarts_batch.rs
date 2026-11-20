//! Prepared hit indices for the opt-in RDKit SMARTS profiles.

use super::{PyMoleculeBatch, PySmartsQuery};
use chematic_smarts::{
    AtomPrimitive, AtomQuery, MatchOutcome, RdkitParityConfig, RdkitParityError,
    RdkitRingCountModel, has_match_rdkit_parity_bounded,
    rdkit_parity_match::rdkit_relevant_ring_count_presence,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

type Refusals = Vec<(usize, &'static str)>;
type BatchMatches = (Vec<usize>, Refusals);

#[derive(Default)]
pub(super) struct RelevantRingCounts {
    masks: Vec<[u64; 4]>,
    refusals: Refusals,
}

impl RelevantRingCounts {
    fn prepare(batch: &PyMoleculeBatch, config: &RdkitParityConfig) -> Self {
        let mut counts = Self {
            masks: Vec::with_capacity(batch.molecules.len()),
            refusals: Vec::new(),
        };
        for (index, mol) in batch.molecules.iter().enumerate() {
            let suppressed = chematic_smiles::rdkit_hydrogen_suppressed(mol);
            let target = suppressed.as_ref().unwrap_or(mol);
            match rdkit_relevant_ring_count_presence(target, config) {
                Ok(mask) => counts.masks.push(mask),
                Err(error) => {
                    counts.refusals.push((index, parity_error_status(&error).1));
                    counts.masks.push([0; 4]);
                }
            }
        }
        counts
    }

    fn matches(&self, count: u8) -> BatchMatches {
        let word = usize::from(count / 64);
        let bit = 1 << (count % 64);
        let hits = self
            .masks
            .iter()
            .enumerate()
            .filter_map(|(i, mask)| (mask[word] & bit != 0).then_some(i))
            .collect();
        (hits, self.refusals.clone())
    }
}

pub(crate) fn parity_profile(profile: &str) -> PyResult<RdkitParityConfig> {
    let ring_count_model = match profile {
        "2026.03.6" => RdkitRingCountModel::SymmetrizedSssr,
        "2026.09.1" => RdkitRingCountModel::RelevantCycles,
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown profile '{other}': expected '2026.03.6' or '2026.09.1'"
            )));
        }
    };
    Ok(RdkitParityConfig {
        use_rdkit_parity_aromaticity: true,
        ring_count_model,
        ..Default::default()
    })
}

pub(crate) fn parity_error_status(error: &RdkitParityError) -> (&'static str, &'static str) {
    match error {
        RdkitParityError::RingModelAmbiguous { .. } => {
            ("typed_unsupported", "ring_model_ambiguous")
        }
        RdkitParityError::RingModelBudgetExceeded { .. } => {
            ("typed_refusal", "ring_model_budget_exceeded")
        }
        RdkitParityError::Aromaticity(_) => ("typed_refusal", "aromaticity_failure"),
    }
}

#[pymethods]
impl PySmartsQuery {
    /// Sorted hit indices and (index, reason) refusals for the RDKit profile.
    /// For isolated positive ring-count predicates in the 2026.09.1 profile,
    /// prepare the target ring-count properties on first use and reuse them.
    #[pyo3(signature = (batch, profile = "2026.03.6", use_chirality = false))]
    fn match_indices_rdkit_parity_batch(
        &self,
        py: Python<'_>,
        batch: &PyMoleculeBatch,
        profile: &str,
        use_chirality: bool,
    ) -> PyResult<BatchMatches> {
        let mut config = parity_profile(profile)?;
        config.base.use_chirality = use_chirality;
        Ok(py.detach(|| {
            if config.ring_count_model == RdkitRingCountModel::RelevantCycles
                && let [atom] = self.query.atoms.as_slice()
                && self.query.bonds.is_empty()
                && let AtomQuery::Primitive(AtomPrimitive::RingCount(count @ 1..=u8::MAX)) =
                    atom.query
            {
                return batch
                    .relevant_ring_counts
                    .get_or_init(|| RelevantRingCounts::prepare(batch, &config))
                    .matches(count);
            }
            let mut hits = Vec::new();
            let mut refusals = Vec::new();
            for (index, mol) in batch.molecules.iter().enumerate() {
                let suppressed = chematic_smiles::rdkit_hydrogen_suppressed(mol);
                let target = suppressed.as_ref().unwrap_or(mol);
                match has_match_rdkit_parity_bounded(&self.query, target, &config) {
                    Ok(MatchOutcome::Found) => hits.push(index),
                    Ok(MatchOutcome::NotFound) => {}
                    Ok(MatchOutcome::BudgetExhausted) => {
                        refusals.push((index, "matcher_budget_exceeded"));
                    }
                    Err(error) => refusals.push((index, parity_error_status(&error).1)),
                }
            }
            (hits, refusals)
        }))
    }
}
