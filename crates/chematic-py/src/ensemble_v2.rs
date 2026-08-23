//! Python binding for `chematic_3d::embed_ensemble_v2` (A2.1).
//!
//! Calls `embed_ensemble_v2` directly on the caller's own `Mol.inner` --
//! never canonicalizes/reparses first, same convention as
//! `embed_pipeline_v2` (see `pipeline_v2.rs`'s module doc, issue #172).
//!
//! # Error asymmetry with `embed_pipeline_v2`
//!
//! `embed_pipeline_v2` fails the whole call when a single embed doesn't
//! work out. `embed_ensemble_v2` is different: every per-attempt outcome --
//! including an ensemble where every attempt failed and zero conformers
//! were kept -- is a normal, fully-diagnosable `Ok(EnsembleV2Result)`, with
//! per-attempt detail in `attempts`. Its `Result` only rejects a config
//! that could never succeed regardless of the molecule or how many
//! attempts run (currently: an invalid `rmsd_threshold`). So there is no
//! `EnsembleV2Error` exception type here -- a rejected config just raises a
//! plain `ValueError` (mirrors `parse_stereo_policy` et al. in
//! `pipeline_v2.rs`, which do the same for bad policy strings). A Python
//! caller must not assume "no exception raised" implies "got at least one
//! conformer" -- check `len(result["conformers"])`.
//!
//! # Never cross-comparing MMFF94 vs UFF energy
//!
//! `embed_ensemble_v2` itself already scopes energy ranking and duplicate
//! pruning to within each `actual_force_field_used` group, and reports
//! `mixed_force_field` when kept conformers span more than one group (see
//! that function's own module doc in `chematic-3d` for the full
//! rationale). This binding deliberately adds no top-level "best
//! conformer" or flattened, globally energy-sorted field on top of that --
//! doing so would silently reintroduce the exact cross-group comparison
//! `embed_ensemble_v2` was built to avoid. `result["conformers"]` is
//! already correctly ordered (group-then-energy); per-conformer energy and
//! force-field context lives in the parallel `result["conformer_provenance"]`
//! list, and the full untrimmed record lives in `result["attempts"]`.

use chematic_3d::{ConformerAttempt, ConformerDisposition, EnsembleV2Result, embed_ensemble_v2};
use chematic_core::Molecule;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::pipeline_v2::{
    PyPipelineV2Config, coords_to_vec, failure_to_dict, force_field_bridge_error_dict,
    force_field_policy_str, snake_case_debug,
};

// ---------------------------------------------------------------------------
// EnsembleV2Config
// ---------------------------------------------------------------------------

/// Configuration for `Mol.conformer_ensemble_v2()`. Unlike
/// `PipelineV2Config`, this constructor is infallible: the Rust
/// `EnsembleV2Config` struct has no invariant of its own to parse or
/// reject at construction time (its one fallible field, `rmsd_threshold`,
/// is validated by `embed_ensemble_v2` itself, at call time -- see this
/// module's doc for why that raises a plain `ValueError` rather than a
/// dedicated exception type).
#[pyclass(name = "EnsembleV2Config", from_py_object)]
#[derive(Clone)]
pub struct PyEnsembleV2Config {
    pub(crate) inner: chematic_3d::EnsembleV2Config,
}

#[pymethods]
impl PyEnsembleV2Config {
    /// `per_conformer`: a `PipelineV2Config`, reused verbatim for every
    /// attempt except its `embed_seed`, which is overridden per-attempt by
    /// a value derived from `base_seed`.
    ///
    /// `count`: number of independent embedding attempts (before RMSD
    /// pruning) -- the kept ensemble may have fewer conformers than this.
    ///
    /// `base_seed`: attempt `i` uses a seed deterministically derived from
    /// this value. The same `base_seed` always reproduces the same
    /// ensemble.
    ///
    /// `rmsd_threshold` (default 0.5 Å): minimum RMSD between kept
    /// conformers within the same force-field group. `0.0` disables
    /// pruning. Must be `0.0` or a positive, finite value -- an invalid
    /// value is accepted here but rejected with `ValueError` when
    /// `conformer_ensemble_v2()` is actually called.
    ///
    /// `use_symmetric_rmsd_pruning` (default `True`): automorphism-aware
    /// RMSD (correct on molecules with interchangeable substituents, e.g.
    /// -CF3, at the cost of enumerating automorphisms per comparison). Set
    /// `False` for cheaper plain Kabsch RMSD.
    ///
    /// `ensemble_timeout_ms` (default `None`): wall-clock budget across
    /// all `count` attempts combined, checked between attempts only.
    /// Distinct from `per_conformer`'s own `total_timeout_ms`, which
    /// budgets a single attempt.
    #[new]
    #[pyo3(signature = (
        per_conformer,
        count,
        base_seed,
        rmsd_threshold = 0.5,
        use_symmetric_rmsd_pruning = true,
        ensemble_timeout_ms = None,
    ))]
    fn new(
        per_conformer: &PyPipelineV2Config,
        count: usize,
        base_seed: u64,
        rmsd_threshold: f64,
        use_symmetric_rmsd_pruning: bool,
        ensemble_timeout_ms: Option<u64>,
    ) -> Self {
        Self {
            inner: chematic_3d::EnsembleV2Config {
                per_conformer: per_conformer.inner.clone(),
                count,
                base_seed,
                rmsd_threshold,
                use_symmetric_rmsd_pruning,
                ensemble_timeout_ms,
            },
        }
    }

    #[getter]
    fn count(&self) -> usize {
        self.inner.count
    }
    #[getter]
    fn base_seed(&self) -> u64 {
        self.inner.base_seed
    }
    #[getter]
    fn rmsd_threshold(&self) -> f64 {
        self.inner.rmsd_threshold
    }
    #[getter]
    fn use_symmetric_rmsd_pruning(&self) -> bool {
        self.inner.use_symmetric_rmsd_pruning
    }
    #[getter]
    fn ensemble_timeout_ms(&self) -> Option<u64> {
        self.inner.ensemble_timeout_ms
    }
    #[getter]
    fn per_conformer(&self) -> PyPipelineV2Config {
        PyPipelineV2Config {
            inner: self.inner.per_conformer.clone(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "EnsembleV2Config(count={}, base_seed={}, rmsd_threshold={}, \
             use_symmetric_rmsd_pruning={}, ensemble_timeout_ms={:?})",
            self.count(),
            self.base_seed(),
            self.rmsd_threshold(),
            self.use_symmetric_rmsd_pruning(),
            self.ensemble_timeout_ms(),
        )
    }
}

// ---------------------------------------------------------------------------
// Result -> PyDict conversion
// ---------------------------------------------------------------------------

fn conformer_disposition_dict<'py>(
    py: Python<'py>,
    d: &ConformerDisposition,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    match d {
        ConformerDisposition::Kept { conformer_index } => {
            dict.set_item("kind", "kept")?;
            dict.set_item("conformer_index", conformer_index)?;
        }
        ConformerDisposition::PrunedAsDuplicate {
            representative_attempt_index,
            rmsd,
            symmetric,
        } => {
            dict.set_item("kind", "pruned_as_duplicate")?;
            dict.set_item("representative_attempt_index", representative_attempt_index)?;
            dict.set_item("rmsd", rmsd)?;
            dict.set_item("symmetric", symmetric)?;
        }
        // `ConformerDisposition` is `#[non_exhaustive]`: a new variant is a
        // real gap in this binding, not something to paper over silently.
        other => {
            return Err(PyValueError::new_err(format!(
                "unhandled ConformerDisposition variant in Python binding: {other:?}"
            )));
        }
    }
    Ok(dict)
}

fn conformer_attempt_dict<'py>(
    py: Python<'py>,
    a: &ConformerAttempt,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("attempt_index", a.attempt_index)?;
    dict.set_item("seed", a.seed)?;
    match &a.outcome {
        Ok(success) => {
            dict.set_item("outcome", "success")?;
            let success_dict = PyDict::new(py);
            success_dict.set_item("energy", success.energy)?;
            success_dict.set_item(
                "actual_force_field_used",
                force_field_policy_str(success.actual_force_field_used),
            )?;
            match &success.fallback_reason {
                Some(e) => success_dict
                    .set_item("fallback_reason", force_field_bridge_error_dict(py, e)?)?,
                None => success_dict.set_item("fallback_reason", py.None())?,
            }
            success_dict.set_item(
                "disposition",
                conformer_disposition_dict(py, &success.disposition)?,
            )?;
            dict.set_item("success", success_dict)?;
            dict.set_item("failure", py.None())?;
        }
        Err(failure) => {
            dict.set_item("outcome", "failure")?;
            dict.set_item("success", py.None())?;
            dict.set_item("failure", failure_to_dict(py, failure)?)?;
        }
    }
    Ok(dict)
}

fn ensemble_v2_result_dict<'py>(
    py: Python<'py>,
    r: &EnsembleV2Result,
) -> PyResult<Bound<'py, PyDict>> {
    let conformer_count = r.ensemble.conformer_count();
    let conformers: Vec<Vec<Vec<f64>>> = (0..conformer_count)
        .map(|i| {
            coords_to_vec(
                r.ensemble
                    .get_conformer(i)
                    .expect("index is within conformer_count()"),
            )
        })
        .collect();

    // Reverse-map Kept{conformer_index} -> the attempt that produced it, so
    // callers never have to manually scan `attempts` to find out what a
    // given `conformers[i]` actually is.
    let mut provenance_by_index: Vec<Option<Bound<'py, PyDict>>> = (0..conformer_count)
        .map(|_| None::<Bound<'py, PyDict>>)
        .collect();
    for attempt in &r.attempts {
        if let Ok(success) = &attempt.outcome
            && let ConformerDisposition::Kept { conformer_index } = &success.disposition
        {
            let entry = PyDict::new(py);
            entry.set_item("attempt_index", attempt.attempt_index)?;
            entry.set_item("seed", attempt.seed)?;
            entry.set_item("energy", success.energy)?;
            entry.set_item(
                "actual_force_field_used",
                force_field_policy_str(success.actual_force_field_used),
            )?;
            provenance_by_index[*conformer_index] = Some(entry);
        }
    }
    let conformer_provenance: Vec<Bound<'py, PyDict>> = provenance_by_index
        .into_iter()
        .enumerate()
        .map(|(i, entry)| {
            entry.unwrap_or_else(|| {
                panic!(
                    "conformer {i} of {conformer_count} has no Kept disposition in `attempts` \
                     -- embed_ensemble_v2's own invariant is broken, not a case to paper over"
                )
            })
        })
        .collect();

    let attempts: Vec<Bound<'py, PyDict>> = r
        .attempts
        .iter()
        .map(|a| conformer_attempt_dict(py, a))
        .collect::<PyResult<Vec<_>>>()?;

    let dict = PyDict::new(py);
    dict.set_item("conformers", conformers)?;
    dict.set_item("conformer_provenance", conformer_provenance)?;
    dict.set_item("attempts", attempts)?;
    dict.set_item("mixed_force_field", r.mixed_force_field)?;
    dict.set_item("termination", snake_case_debug(&r.termination))?;
    dict.set_item("requested_count", r.requested_count)?;
    Ok(dict)
}

// ---------------------------------------------------------------------------
// Entry point called from `Mol::conformer_ensemble_v2` in mol_methods.rs
// ---------------------------------------------------------------------------

pub fn run_embed_ensemble_v2<'py>(
    py: Python<'py>,
    mol: &Molecule,
    config: &PyEnsembleV2Config,
) -> PyResult<Bound<'py, PyDict>> {
    match embed_ensemble_v2(mol, &config.inner) {
        Ok(result) => ensemble_v2_result_dict(py, &result),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}
