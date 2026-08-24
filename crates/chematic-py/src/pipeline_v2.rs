//! Python binding for `chematic_3d::pipeline_v2::embed_pipeline_v2`.
//!
//! Calls `embed_pipeline_v2` directly on the caller's own `Mol.inner` — never
//! canonicalizes/reparses first (see issue #172; `conformer_ensemble`'s fix in
//! `mol_methods.rs` is the precedent this follows).
//!
//! `PipelineV2Config` deliberately has no `Default` on the Rust side (every
//! field, especially `force_field_policy`/`stereo_policy`/`ring_torsion_policy`,
//! must be an explicit judgment call — see that struct's own doc comment). The
//! Python-facing `PipelineV2Config` constructor mirrors that: every field is a
//! required keyword argument. `PipelineV2Config.safe()` is a convenience
//! classmethod, but it still requires `force_field_policy`/`stereo_policy`/
//! `ring_torsion_policy` explicitly rather than defaulting them.

use chematic_3d::distance_geometry_v2::EmbedParameters;
use chematic_3d::etkdg_knowledge::TorsionOptimizationConfig;
use chematic_3d::minimize::ForceFieldPolicy;
use chematic_3d::pipeline_v2::{
    self as pv2, PipelineV2Failure, PipelineV2Result, RingTorsionApplicationPolicy, StereoPolicy,
};
use chematic_core::Molecule;
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

create_exception!(
    chematic,
    PipelineV2Error,
    PyValueError,
    "A failed Mol.embed_pipeline_v2() call. `.diagnostics` carries the same \
     per-stage partial evidence a Rust caller sees on `PipelineV2Failure` -- \
     `diagnostics['last_known_coords']` is diagnostic only, never a usable result."
);

// ---------------------------------------------------------------------------
// String <-> enum conversions (the stable public contract for policy fields)
// ---------------------------------------------------------------------------

fn parse_stereo_policy(s: &str) -> PyResult<StereoPolicy> {
    match s {
        "ignore" => Ok(StereoPolicy::Ignore),
        "verify_only" => Ok(StereoPolicy::VerifyOnly),
        "repair_and_verify" => Ok(StereoPolicy::RepairAndVerify),
        other => Err(PyValueError::new_err(format!(
            "unknown stereo_policy {other:?} -- expected one of: \
             \"ignore\", \"verify_only\", \"repair_and_verify\""
        ))),
    }
}

fn stereo_policy_str(p: StereoPolicy) -> &'static str {
    match p {
        StereoPolicy::Ignore => "ignore",
        StereoPolicy::VerifyOnly => "verify_only",
        StereoPolicy::RepairAndVerify => "repair_and_verify",
    }
}

fn parse_ring_torsion_policy(s: &str) -> PyResult<RingTorsionApplicationPolicy> {
    match s {
        "fail_closed" => Ok(RingTorsionApplicationPolicy::FailClosed),
        "diagnostic_only" => Ok(RingTorsionApplicationPolicy::DiagnosticOnly),
        other => Err(PyValueError::new_err(format!(
            "unknown ring_torsion_policy {other:?} -- expected one of: \
             \"fail_closed\", \"diagnostic_only\""
        ))),
    }
}

fn ring_torsion_policy_str(p: RingTorsionApplicationPolicy) -> &'static str {
    match p {
        RingTorsionApplicationPolicy::FailClosed => "fail_closed",
        RingTorsionApplicationPolicy::DiagnosticOnly => "diagnostic_only",
    }
}

fn parse_force_field_policy(s: &str) -> PyResult<ForceFieldPolicy> {
    match s {
        "mmff94_bond_angle_strict" => Ok(ForceFieldPolicy::Mmff94BondAngleStrict),
        "mmff94_with_uff_fallback" => Ok(ForceFieldPolicy::Mmff94WithUffFallback),
        "uff_only" => Ok(ForceFieldPolicy::UffOnly),
        "dreiding" => Ok(ForceFieldPolicy::Dreiding),
        "none" => Ok(ForceFieldPolicy::None),
        other => Err(PyValueError::new_err(format!(
            "unknown force_field_policy {other:?} -- expected one of: \
             \"mmff94_bond_angle_strict\", \"mmff94_with_uff_fallback\", \"uff_only\", \
             \"dreiding\", \"none\""
        ))),
    }
}

pub(crate) fn force_field_policy_str(p: ForceFieldPolicy) -> &'static str {
    match p {
        ForceFieldPolicy::Mmff94BondAngleStrict => "mmff94_bond_angle_strict",
        ForceFieldPolicy::Mmff94WithUffFallback => "mmff94_with_uff_fallback",
        ForceFieldPolicy::UffOnly => "uff_only",
        ForceFieldPolicy::Dreiding => "dreiding",
        ForceFieldPolicy::None => "none",
    }
}

/// `format!("{value:?}")` -> `snake_case`, for the many small fieldless enums in
/// this pipeline (e.g. `EmbedFailureCause::BoundsSmoothingFailed` ->
/// `"bounds_smoothing_failed"`). Not used for enums with data-carrying variants
/// (those get hand-written conversions above/below so the payload isn't lost).
pub(crate) fn snake_case_debug<T: std::fmt::Debug>(value: &T) -> String {
    let debug = format!("{value:?}");
    let mut out = String::with_capacity(debug.len() + 4);
    for (i, c) in debug.chars().enumerate() {
        if c.is_uppercase() {
            if i != 0 {
                out.push('_');
            }
            out.extend(c.to_lowercase());
        } else {
            out.push(c);
        }
    }
    out
}

// ---------------------------------------------------------------------------
// PipelineV2Config
// ---------------------------------------------------------------------------

/// Configuration for `Mol.embed_pipeline_v2()`. Every field is required
/// (matching the Rust `PipelineV2Config`'s deliberate lack of a `Default` --
/// force-field/stereo/ring-torsion policy are judgment calls, never a hidden
/// default). Use `PipelineV2Config.safe(...)` for a convenience constructor
/// that still requires those three policies explicitly.
#[pyclass(name = "PipelineV2Config", from_py_object)]
#[derive(Clone)]
pub struct PyPipelineV2Config {
    pub(crate) inner: pv2::PipelineV2Config,
}

#[pymethods]
impl PyPipelineV2Config {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        embed_seed,
        max_attempts,
        embed_timeout_ms,
        use_exp_torsions,
        use_small_ring_torsions,
        use_macrocycle_torsions,
        use_macrocycle_14_bounds,
        include_legacy_torsion_heuristic,
        stereo_policy,
        fail_on_unevaluable_stereo,
        force_field_policy,
        force_field_max_iterations,
        gate_mmff94_torsion_oop,
        gate_mmff94_stretch_bend,
        ring_torsion_policy,
        total_timeout_ms,
        enforce_chirality = false,
        expand_implicit_h_through_pipeline = false,
    ))]
    fn new(
        embed_seed: u64,
        max_attempts: usize,
        embed_timeout_ms: Option<u64>,
        use_exp_torsions: bool,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle_14_bounds: bool,
        include_legacy_torsion_heuristic: bool,
        stereo_policy: &str,
        fail_on_unevaluable_stereo: bool,
        force_field_policy: &str,
        force_field_max_iterations: usize,
        gate_mmff94_torsion_oop: bool,
        gate_mmff94_stretch_bend: bool,
        ring_torsion_policy: &str,
        total_timeout_ms: Option<u64>,
        // Trailing, defaulted (`false`, matching `EmbedParameters::default()`) so
        // existing callers' positional/keyword calls keep working unchanged --
        // added after `enforce_chirality` (v0.14.0, issue #285's E/Z bound fix)
        // gained a real production effect. See `distance_geometry_v2.rs`'s module
        // doc for what it does; see `pipeline_v2.rs`'s for why it's compatible
        // with `stereo_policy="ignore"`/`"verify_only"` but not
        // `"repair_and_verify"` (raises `PipelineV2Error` at validate-config
        // otherwise).
        enforce_chirality: bool,
        // Same trailing-defaulted precedent again (issue #291/#383). Requires
        // `enforce_chirality=True` (raises `PipelineV2Error` at validate-config
        // otherwise) -- see `PipelineV2Config::expand_implicit_h_through_pipeline`'s
        // own Rust doc for what it does. Prefer `PipelineV2Config.stereo_safe(...)`
        // over setting this flag alone: it only works correctly combined with
        // `stereo_policy="repair_and_verify"` and `enforce_chirality=True`, and
        // `stereo_safe` sets all three together so a caller can't set one but
        // forget another.
        expand_implicit_h_through_pipeline: bool,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: pv2::PipelineV2Config {
                embed: EmbedParameters {
                    random_seed: embed_seed,
                    max_attempts,
                    timeout_ms: embed_timeout_ms,
                    use_exp_torsions,
                    use_small_ring_torsions,
                    use_macrocycle_torsions,
                    use_macrocycle_14_bounds,
                    enforce_chirality,
                    ..EmbedParameters::default()
                },
                torsion_optimization: TorsionOptimizationConfig::default(),
                include_legacy_torsion_heuristic,
                stereo_policy: parse_stereo_policy(stereo_policy)?,
                fail_on_unevaluable_stereo,
                force_field_policy: parse_force_field_policy(force_field_policy)?,
                force_field_max_iterations,
                gate_mmff94_torsion_oop,
                gate_mmff94_stretch_bend,
                ring_torsion_policy: parse_ring_torsion_policy(ring_torsion_policy)?,
                total_timeout_ms,
                expand_implicit_h_through_pipeline,
            },
        })
    }

    /// Convenience constructor. `force_field` (the force-field policy),
    /// `stereo_policy`, and `ring_torsion_policy` are still required, explicit
    /// arguments -- never hidden defaults -- everything else takes a
    /// conservative default (every torsion-knowledge flag off,
    /// `fail_on_unevaluable_stereo=False`, no timeouts).
    #[staticmethod]
    #[pyo3(signature = (
        force_field,
        stereo_policy,
        ring_torsion_policy,
        fail_on_unevaluable_stereo = false,
        embed_seed = 0xC0FF_EE42_D157_6E02,
        max_attempts = 8,
        embed_timeout_ms = None,
        use_exp_torsions = false,
        use_small_ring_torsions = false,
        use_macrocycle_torsions = false,
        use_macrocycle_14_bounds = false,
        include_legacy_torsion_heuristic = false,
        force_field_max_iterations = 200,
        gate_mmff94_torsion_oop = false,
        gate_mmff94_stretch_bend = false,
        total_timeout_ms = None,
        enforce_chirality = false,
        expand_implicit_h_through_pipeline = false,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn safe(
        force_field: &str,
        stereo_policy: &str,
        ring_torsion_policy: &str,
        fail_on_unevaluable_stereo: bool,
        embed_seed: u64,
        max_attempts: usize,
        embed_timeout_ms: Option<u64>,
        use_exp_torsions: bool,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle_14_bounds: bool,
        include_legacy_torsion_heuristic: bool,
        force_field_max_iterations: usize,
        gate_mmff94_torsion_oop: bool,
        gate_mmff94_stretch_bend: bool,
        total_timeout_ms: Option<u64>,
        enforce_chirality: bool,
        expand_implicit_h_through_pipeline: bool,
    ) -> PyResult<Self> {
        Self::new(
            embed_seed,
            max_attempts,
            embed_timeout_ms,
            use_exp_torsions,
            use_small_ring_torsions,
            use_macrocycle_torsions,
            use_macrocycle_14_bounds,
            include_legacy_torsion_heuristic,
            stereo_policy,
            fail_on_unevaluable_stereo,
            force_field,
            force_field_max_iterations,
            gate_mmff94_torsion_oop,
            gate_mmff94_stretch_bend,
            ring_torsion_policy,
            total_timeout_ms,
            enforce_chirality,
            expand_implicit_h_through_pipeline,
        )
    }

    /// Convenience constructor for the "stereo-safe" configuration (issue
    /// #291/#383): sets `stereo_policy="repair_and_verify"`,
    /// `enforce_chirality=True`, and `expand_implicit_h_through_pipeline=True`
    /// together -- the exact combination measured to correctly handle
    /// ring-fused declared stereocenters (e.g. testosterone, cholesterol) that
    /// `enforce_chirality` alone cannot repair. These three only work
    /// correctly as a set; prefer this over setting them individually via
    /// `safe(...)`/the constructor, where forgetting one silently falls back
    /// to a configuration issue #291 measured as unsound for that molecule
    /// class. `force_field`/`ring_torsion_policy` are still required,
    /// explicit arguments; everything else takes the same conservative
    /// defaults `safe(...)` does.
    #[staticmethod]
    #[pyo3(signature = (
        force_field,
        ring_torsion_policy,
        fail_on_unevaluable_stereo = false,
        embed_seed = 0xC0FF_EE42_D157_6E02,
        max_attempts = 8,
        embed_timeout_ms = None,
        use_exp_torsions = false,
        use_small_ring_torsions = false,
        use_macrocycle_torsions = false,
        use_macrocycle_14_bounds = false,
        include_legacy_torsion_heuristic = false,
        force_field_max_iterations = 200,
        gate_mmff94_torsion_oop = false,
        gate_mmff94_stretch_bend = false,
        total_timeout_ms = None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn stereo_safe(
        force_field: &str,
        ring_torsion_policy: &str,
        fail_on_unevaluable_stereo: bool,
        embed_seed: u64,
        max_attempts: usize,
        embed_timeout_ms: Option<u64>,
        use_exp_torsions: bool,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle_14_bounds: bool,
        include_legacy_torsion_heuristic: bool,
        force_field_max_iterations: usize,
        gate_mmff94_torsion_oop: bool,
        gate_mmff94_stretch_bend: bool,
        total_timeout_ms: Option<u64>,
    ) -> PyResult<Self> {
        let mut inner = pv2::PipelineV2Config::stereo_safe(parse_force_field_policy(force_field)?);
        inner.ring_torsion_policy = parse_ring_torsion_policy(ring_torsion_policy)?;
        inner.fail_on_unevaluable_stereo = fail_on_unevaluable_stereo;
        inner.embed.random_seed = embed_seed;
        inner.embed.max_attempts = max_attempts;
        inner.embed.timeout_ms = embed_timeout_ms;
        inner.embed.use_exp_torsions = use_exp_torsions;
        inner.embed.use_small_ring_torsions = use_small_ring_torsions;
        inner.embed.use_macrocycle_torsions = use_macrocycle_torsions;
        inner.embed.use_macrocycle_14_bounds = use_macrocycle_14_bounds;
        inner.include_legacy_torsion_heuristic = include_legacy_torsion_heuristic;
        inner.force_field_max_iterations = force_field_max_iterations;
        inner.gate_mmff94_torsion_oop = gate_mmff94_torsion_oop;
        inner.gate_mmff94_stretch_bend = gate_mmff94_stretch_bend;
        inner.total_timeout_ms = total_timeout_ms;
        Ok(Self { inner })
    }

    #[getter]
    fn embed_seed(&self) -> u64 {
        self.inner.embed.random_seed
    }
    #[getter]
    fn max_attempts(&self) -> usize {
        self.inner.embed.max_attempts
    }
    #[getter]
    fn embed_timeout_ms(&self) -> Option<u64> {
        self.inner.embed.timeout_ms
    }
    #[getter]
    fn use_exp_torsions(&self) -> bool {
        self.inner.embed.use_exp_torsions
    }
    #[getter]
    fn use_small_ring_torsions(&self) -> bool {
        self.inner.embed.use_small_ring_torsions
    }
    #[getter]
    fn use_macrocycle_torsions(&self) -> bool {
        self.inner.embed.use_macrocycle_torsions
    }
    #[getter]
    fn use_macrocycle_14_bounds(&self) -> bool {
        self.inner.embed.use_macrocycle_14_bounds
    }
    #[getter]
    fn include_legacy_torsion_heuristic(&self) -> bool {
        self.inner.include_legacy_torsion_heuristic
    }
    #[getter]
    fn stereo_policy(&self) -> &'static str {
        stereo_policy_str(self.inner.stereo_policy)
    }
    #[getter]
    fn fail_on_unevaluable_stereo(&self) -> bool {
        self.inner.fail_on_unevaluable_stereo
    }
    #[getter]
    fn force_field_policy(&self) -> &'static str {
        force_field_policy_str(self.inner.force_field_policy)
    }
    #[getter]
    fn force_field_max_iterations(&self) -> usize {
        self.inner.force_field_max_iterations
    }
    #[getter]
    fn gate_mmff94_torsion_oop(&self) -> bool {
        self.inner.gate_mmff94_torsion_oop
    }
    #[getter]
    fn gate_mmff94_stretch_bend(&self) -> bool {
        self.inner.gate_mmff94_stretch_bend
    }
    #[getter]
    fn ring_torsion_policy(&self) -> &'static str {
        ring_torsion_policy_str(self.inner.ring_torsion_policy)
    }
    #[getter]
    fn total_timeout_ms(&self) -> Option<u64> {
        self.inner.total_timeout_ms
    }
    #[getter]
    fn enforce_chirality(&self) -> bool {
        self.inner.embed.enforce_chirality
    }
    #[getter]
    fn expand_implicit_h_through_pipeline(&self) -> bool {
        self.inner.expand_implicit_h_through_pipeline
    }

    fn __repr__(&self) -> String {
        format!(
            "PipelineV2Config(force_field_policy={:?}, stereo_policy={:?}, \
             ring_torsion_policy={:?}, embed_seed={}, max_attempts={})",
            self.force_field_policy(),
            self.stereo_policy(),
            self.ring_torsion_policy(),
            self.embed_seed(),
            self.max_attempts(),
        )
    }
}

// ---------------------------------------------------------------------------
// Result / failure -> PyDict conversion
// ---------------------------------------------------------------------------

pub(crate) fn coords_to_vec(coords: &chematic_3d::coords::Coords3D) -> Vec<Vec<f64>> {
    coords.points.iter().map(|p| vec![p.x, p.y, p.z]).collect()
}

fn embed_stats_dict<'py>(
    py: Python<'py>,
    stats: &chematic_3d::distance_geometry_v2::EmbedStats,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("attempts_used", stats.attempts_used)?;
    let failure_counts = PyDict::new(py);
    for (cause, count) in &stats.failure_counts {
        failure_counts.set_item(snake_case_debug(cause), count)?;
    }
    d.set_item("failure_counts", failure_counts)?;
    d.set_item(
        "negative_eigenvalues_beyond_embedding_dim",
        stats.negative_eigenvalues_beyond_embedding_dim,
    )?;
    d.set_item(
        "max_negative_eigenvalue_magnitude",
        stats.max_negative_eigenvalue_magnitude,
    )?;
    d.set_item(
        "last_smoothing_invariants_ok",
        stats.last_smoothing_invariants_ok,
    )?;
    d.set_item("used_random_coords", stats.used_random_coords)?;
    d.set_item("adjustments_applied", stats.adjustments_applied)?;
    Ok(d)
}

fn bound_adjustment_dict<'py>(
    py: Python<'py>,
    a: &chematic_3d::etkdg_knowledge::PairBoundAdjustment,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("atom1", a.atom_pair.0.0)?;
    d.set_item("atom2", a.atom_pair.1.0)?;
    d.set_item("old_lower", a.old_lower)?;
    d.set_item("new_lower", a.new_lower)?;
    d.set_item("old_upper", a.old_upper)?;
    d.set_item("new_upper", a.new_upper)?;
    d.set_item("rule_id", &a.rule_id)?;
    d.set_item("source", snake_case_debug(&a.source))?;
    d.set_item("ring_size", a.ring_size)?;
    d.set_item("reason", &a.reason)?;
    Ok(d)
}

fn fourier_term_dict<'py>(
    py: Python<'py>,
    t: &chematic_3d::etkdg_knowledge::FourierTorsionTerm,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("periodicity", t.periodicity)?;
    d.set_item("phase_deg", t.phase_deg)?;
    d.set_item("amplitude", t.amplitude)?;
    Ok(d)
}

fn torsion_potential_dict<'py>(
    py: Python<'py>,
    p: &chematic_3d::etkdg_knowledge::TorsionPotential,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("atoms", p.atoms.iter().map(|a| a.0).collect::<Vec<u32>>())?;
    d.set_item("central_bond", (p.central_bond.0.0, p.central_bond.1.0))?;
    d.set_item("source", snake_case_debug(&p.source))?;
    d.set_item("rule_id", &p.rule_id)?;
    let terms = p
        .terms
        .iter()
        .map(|t| fourier_term_dict(py, t))
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("terms", terms)?;
    d.set_item("ring_size", p.ring_size)?;
    Ok(d)
}

fn torsion_diagnostic_dict<'py>(
    py: Python<'py>,
    diag: &chematic_3d::etkdg_knowledge::TorsionKnowledgeDiagnostic,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item(
        "central_bond",
        (diag.central_bond.0.0, diag.central_bond.1.0),
    )?;
    d.set_item("kind", snake_case_debug(&diag.kind))?;
    d.set_item("message", &diag.message)?;
    d.set_item("candidate_rule_ids", &diag.candidate_rule_ids)?;
    Ok(d)
}

fn torsion_knowledge_report_dict<'py>(
    py: Python<'py>,
    r: &chematic_3d::etkdg_knowledge::TorsionKnowledgeReport,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let potentials = r
        .potentials
        .iter()
        .map(|p| torsion_potential_dict(py, p))
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("potentials", potentials)?;
    d.set_item("matched_rule_ids", &r.matched_rule_ids)?;
    d.set_item(
        "unmatched_rotatable_bonds",
        r.unmatched_rotatable_bonds
            .iter()
            .map(|(a, b)| (a.0, b.0))
            .collect::<Vec<(u32, u32)>>(),
    )?;
    let ambiguous = r
        .ambiguous_matches
        .iter()
        .map(|diag| torsion_diagnostic_dict(py, diag))
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("ambiguous_matches", ambiguous)?;
    let skipped = r
        .skipped_bonds
        .iter()
        .map(|diag| torsion_diagnostic_dict(py, diag))
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("skipped_bonds", skipped)?;
    Ok(d)
}

fn ring_torsion_evidence_dict<'py>(
    py: Python<'py>,
    e: &chematic_3d::pipeline_v2::RingTorsionEvidence,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let potentials = e
        .potentials
        .iter()
        .map(|p| {
            let pd = PyDict::new(py);
            pd.set_item("rule_id", &p.rule_id)?;
            pd.set_item("central_bond", (p.central_bond.0.0, p.central_bond.1.0))?;
            pd.set_item("source", snake_case_debug(&p.source))?;
            pd.set_item("applied_to_geometry", p.applied_to_geometry)?;
            Ok(pd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("potentials", potentials)?;
    d.set_item("diagnostic_only", e.diagnostic_only)?;
    d.set_item("n_applied", e.n_applied())?;
    d.set_item("n_scored_only", e.n_scored_only())?;
    Ok(d)
}

fn torsion_optimization_report_dict<'py>(
    py: Python<'py>,
    r: &chematic_3d::etkdg_knowledge::TorsionOptimizationReport,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("energy_before", r.energy_before)?;
    d.set_item("energy_after", r.energy_after)?;
    d.set_item("iterations_used", r.iterations_used)?;
    d.set_item("converged", r.converged)?;
    d.set_item("max_bond_length_delta", r.max_bond_length_delta)?;
    d.set_item("max_ring_closure_delta", r.max_ring_closure_delta)?;
    d.set_item("rotated_bond_count", r.rotated_bond_count)?;
    Ok(d)
}

fn stereo_status_dict<'py>(
    py: Python<'py>,
    status: chematic_3d::stereo_constraints::StereoStatus,
) -> PyResult<Bound<'py, PyDict>> {
    use chematic_3d::stereo_constraints::StereoStatus;
    let d = PyDict::new(py);
    match status {
        StereoStatus::Satisfied => {
            d.set_item("status", "satisfied")?;
            d.set_item("rejection_reason", py.None())?;
        }
        StereoStatus::Violated => {
            d.set_item("status", "violated")?;
            d.set_item("rejection_reason", py.None())?;
        }
        StereoStatus::Unevaluable(reason) => {
            d.set_item("status", "unevaluable")?;
            d.set_item("rejection_reason", snake_case_debug(&reason))?;
        }
    }
    Ok(d)
}

fn stereo_verification_dict<'py>(
    py: Python<'py>,
    v: &chematic_3d::stereo_constraints::StereoVerification,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let tetrahedral = v
        .tetrahedral
        .iter()
        .map(|r| {
            let rd = stereo_status_dict(py, r.status)?;
            rd.set_item("atom", r.atom.0)?;
            Ok(rd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("tetrahedral", tetrahedral)?;
    let double_bond = v
        .double_bond
        .iter()
        .map(|r| {
            let rd = stereo_status_dict(py, r.status)?;
            rd.set_item("bond", r.bond.0)?;
            Ok(rd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("double_bond", double_bond)?;
    d.set_item("n_declared", v.n_declared())?;
    d.set_item("n_satisfied", v.n_satisfied())?;
    d.set_item("n_violations", v.n_violations())?;
    d.set_item("n_unevaluable", v.n_unevaluable())?;
    d.set_item("is_fully_satisfied", v.is_fully_satisfied())?;
    Ok(d)
}

fn repaired_element_dict<'py>(
    py: Python<'py>,
    r: &chematic_3d::stereo_constraints::RepairedElement,
) -> PyResult<Bound<'py, PyDict>> {
    use chematic_3d::stereo_constraints::StereoElement;
    let d = PyDict::new(py);
    match r.element {
        StereoElement::Tetrahedral(atom) => {
            d.set_item("element_kind", "tetrahedral")?;
            d.set_item("atom", atom.0)?;
        }
        StereoElement::DoubleBond(bond) => {
            d.set_item("element_kind", "double_bond")?;
            d.set_item("bond", bond.0)?;
        }
    }
    d.set_item("atoms_moved", r.atoms_moved)?;
    d.set_item("max_displacement", r.max_displacement)?;
    Ok(d)
}

fn stereo_repair_summary_dict<'py>(
    py: Python<'py>,
    s: &chematic_3d::pipeline_v2::StereoRepairSummary,
) -> PyResult<Bound<'py, PyDict>> {
    use chematic_3d::stereo_constraints::StereoElement;
    let d = PyDict::new(py);
    let repaired = s
        .repaired
        .iter()
        .map(|r| repaired_element_dict(py, r))
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("repaired", repaired)?;
    let failures = s
        .failures
        .iter()
        .map(|(elem, reason)| {
            let fd = PyDict::new(py);
            match elem {
                StereoElement::Tetrahedral(atom) => {
                    fd.set_item("element_kind", "tetrahedral")?;
                    fd.set_item("atom", atom.0)?;
                }
                StereoElement::DoubleBond(bond) => {
                    fd.set_item("element_kind", "double_bond")?;
                    fd.set_item("bond", bond.0)?;
                }
            }
            fd.set_item("reason", snake_case_debug(reason))?;
            Ok(fd)
        })
        .collect::<PyResult<Vec<_>>>()?;
    d.set_item("failures", failures)?;
    Ok(d)
}

fn mmff94_missing_term_dict<'py>(
    py: Python<'py>,
    t: &chematic_3d::minimize::Mmff94MissingTerm,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("kind", snake_case_debug(&t.kind))?;
    d.set_item("atoms", t.atoms.iter().map(|a| a.0).collect::<Vec<u32>>())?;
    d.set_item("description", &t.description)?;
    Ok(d)
}

fn mmff94_coverage_dict<'py>(
    py: Python<'py>,
    r: &chematic_3d::minimize::Mmff94CoverageReport,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("bonds_total", r.bonds_total)?;
    d.set_item(
        "bonds_missing",
        r.bonds_missing
            .iter()
            .map(|t| mmff94_missing_term_dict(py, t))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    d.set_item("angles_total", r.angles_total)?;
    d.set_item(
        "angles_missing",
        r.angles_missing
            .iter()
            .map(|t| mmff94_missing_term_dict(py, t))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    d.set_item("torsions_total", r.torsions_total)?;
    d.set_item(
        "torsions_missing",
        r.torsions_missing
            .iter()
            .map(|t| mmff94_missing_term_dict(py, t))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    d.set_item("oop_total", r.oop_total)?;
    d.set_item(
        "oop_missing",
        r.oop_missing
            .iter()
            .map(|t| mmff94_missing_term_dict(py, t))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    d.set_item("stretch_bend_total", r.stretch_bend_total)?;
    d.set_item(
        "stretch_bend_missing",
        r.stretch_bend_missing
            .iter()
            .map(|t| mmff94_missing_term_dict(py, t))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    Ok(d)
}

pub(crate) fn force_field_bridge_error_dict<'py>(
    py: Python<'py>,
    e: &chematic_3d::minimize::ForceFieldBridgeError,
) -> PyResult<Bound<'py, PyDict>> {
    use chematic_3d::minimize::ForceFieldBridgeError;
    let d = PyDict::new(py);
    match e {
        ForceFieldBridgeError::UnsupportedAtomType(msg) => {
            d.set_item("kind", "unsupported_atom_type")?;
            d.set_item("message", msg)?;
        }
        ForceFieldBridgeError::MissingParameters(coverage) => {
            d.set_item("kind", "missing_parameters")?;
            d.set_item("coverage", mmff94_coverage_dict(py, coverage)?)?;
        }
        ForceFieldBridgeError::MinimizationFailed(detail) => {
            d.set_item("kind", "minimization_failed")?;
            d.set_item("policy", force_field_policy_str(detail.policy))?;
            d.set_item("reason", snake_case_debug(&detail.reason))?;
            d.set_item("converged", detail.converged)?;
            d.set_item("iterations", detail.iterations)?;
            d.set_item("max_residual_force", detail.max_residual_force)?;
        }
    }
    Ok(d)
}

fn energy_report_dict<'py>(
    py: Python<'py>,
    e: &chematic_3d::minimize::EnergyReport,
) -> PyResult<Bound<'py, PyDict>> {
    use chematic_3d::minimize::EnergyReport;
    let d = PyDict::new(py);
    match e {
        EnergyReport::Mmff94(b) => {
            d.set_item("kind", "mmff94")?;
            d.set_item("bond", b.bond)?;
            d.set_item("angle", b.angle)?;
            d.set_item("stretch_bend", b.stretch_bend)?;
            d.set_item("torsion", b.torsion)?;
            d.set_item("oop", b.oop)?;
            d.set_item("vdw", b.vdw)?;
            d.set_item("electrostatic", b.electrostatic)?;
            d.set_item("total", b.total)?;
        }
        EnergyReport::Uff { total } => {
            d.set_item("kind", "uff")?;
            d.set_item("total", total)?;
        }
        EnergyReport::Dreiding { total } => {
            d.set_item("kind", "dreiding")?;
            d.set_item("total", total)?;
        }
        EnergyReport::None => {
            d.set_item("kind", "none")?;
            d.set_item("total", 0.0)?;
        }
    }
    Ok(d)
}

fn policy_minimize_result_dict<'py>(
    py: Python<'py>,
    r: &chematic_3d::minimize::PolicyMinimizeResult,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("coords", coords_to_vec(&r.coords))?;
    d.set_item(
        "requested_force_field",
        force_field_policy_str(r.requested_force_field),
    )?;
    d.set_item(
        "actual_force_field_used",
        force_field_policy_str(r.actual_force_field_used),
    )?;
    match &r.fallback_reason {
        Some(e) => d.set_item("fallback_reason", force_field_bridge_error_dict(py, e)?)?,
        None => d.set_item("fallback_reason", py.None())?,
    }
    d.set_item(
        "missing_parameter_classes",
        r.missing_parameter_classes
            .iter()
            .map(|t| mmff94_missing_term_dict(py, t))
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    match &r.coverage {
        Some(c) => d.set_item("coverage", mmff94_coverage_dict(py, c)?)?,
        None => d.set_item("coverage", py.None())?,
    }
    d.set_item("energy_before", energy_report_dict(py, &r.energy_before)?)?;
    d.set_item("energy_after", energy_report_dict(py, &r.energy_after)?)?;
    d.set_item("converged", r.converged)?;
    d.set_item("iterations", r.iterations)?;
    d.set_item("max_residual_force", r.max_residual_force)?;
    d.set_item(
        "starting_geometry",
        r.starting_geometry.map(|g| snake_case_debug(&g)),
    )?;
    Ok(d)
}

fn bounds_conformance_dict<'py>(
    py: Python<'py>,
    b: &chematic_3d::distance_geometry_v2::BoundsConformance,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("n_pairs", b.n_pairs)?;
    d.set_item("n_violations", b.n_violations)?;
    d.set_item("max_rel_violation", b.max_rel_violation)?;
    Ok(d)
}

fn final_validation_dict<'py>(
    py: Python<'py>,
    v: &chematic_3d::pipeline_v2::FinalGeometryValidation,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("all_finite", v.all_finite)?;
    d.set_item("atom_count_unchanged", v.atom_count_unchanged)?;
    d.set_item("worst_bond_length", v.worst_bond_length)?;
    d.set_item("bond_violation_rate_15pct", v.bond_violation_rate_15pct)?;
    d.set_item("bond_violation_rate_50pct", v.bond_violation_rate_50pct)?;
    d.set_item("gross_clash_count", v.gross_clash_count)?;
    d.set_item(
        "bounds_conformance",
        bounds_conformance_dict(py, &v.bounds_conformance)?,
    )?;
    d.set_item("stereo_ok", v.stereo_ok)?;
    d.set_item("torsion_energy_after", v.torsion_energy_after)?;
    d.set_item("ring_closure_delta", v.ring_closure_delta)?;
    d.set_item("sound", v.sound)?;
    Ok(d)
}

fn stage_timings_dict<'py>(
    py: Python<'py>,
    t: &chematic_3d::pipeline_v2::StageTimings,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("torsion_knowledge_ms", t.torsion_knowledge_ms)?;
    d.set_item("bound_adjustment_ms", t.bound_adjustment_ms)?;
    d.set_item("distance_geometry_ms", t.distance_geometry_ms)?;
    d.set_item("torsion_energy_eval_ms", t.torsion_energy_eval_ms)?;
    d.set_item("torsion_optimization_ms", t.torsion_optimization_ms)?;
    d.set_item("stereo_verify_before_ms", t.stereo_verify_before_ms)?;
    d.set_item("stereo_repair_ms", t.stereo_repair_ms)?;
    d.set_item(
        "stereo_verify_after_repair_ms",
        t.stereo_verify_after_repair_ms,
    )?;
    d.set_item("force_field_ms", t.force_field_ms)?;
    d.set_item("final_stereo_verify_ms", t.final_stereo_verify_ms)?;
    d.set_item("final_validation_ms", t.final_validation_ms)?;
    d.set_item("total_ms", t.total_ms)?;
    Ok(d)
}

fn result_to_dict<'py>(py: Python<'py>, r: &PipelineV2Result) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("coords", coords_to_vec(&r.coords))?;
    d.set_item("embed_stats", embed_stats_dict(py, &r.embed_stats)?)?;
    match &r.bound_adjustment_report {
        Some(v) => d.set_item(
            "bound_adjustment_report",
            v.iter()
                .map(|a| bound_adjustment_dict(py, a))
                .collect::<PyResult<Vec<_>>>()?,
        )?,
        None => d.set_item("bound_adjustment_report", py.None())?,
    }
    d.set_item(
        "torsion_knowledge_report",
        torsion_knowledge_report_dict(py, &r.torsion_knowledge_report)?,
    )?;
    d.set_item(
        "ring_torsion_evidence",
        ring_torsion_evidence_dict(py, &r.ring_torsion_evidence)?,
    )?;
    match &r.torsion_optimization_report {
        Some(rep) => d.set_item(
            "torsion_optimization_report",
            torsion_optimization_report_dict(py, rep)?,
        )?,
        None => d.set_item("torsion_optimization_report", py.None())?,
    }
    d.set_item(
        "stereo_before",
        stereo_verification_dict(py, &r.stereo_before)?,
    )?;
    match &r.stereo_repair {
        Some(s) => d.set_item("stereo_repair", stereo_repair_summary_dict(py, s)?)?,
        None => d.set_item("stereo_repair", py.None())?,
    }
    d.set_item(
        "stereo_after_repair",
        stereo_verification_dict(py, &r.stereo_after_repair)?,
    )?;
    d.set_item(
        "force_field",
        policy_minimize_result_dict(py, &r.force_field)?,
    )?;
    d.set_item(
        "final_stereo",
        stereo_verification_dict(py, &r.final_stereo)?,
    )?;
    d.set_item(
        "final_validation",
        final_validation_dict(py, &r.final_validation)?,
    )?;
    d.set_item(
        "elapsed_ms_by_stage",
        stage_timings_dict(py, &r.elapsed_ms_by_stage)?,
    )?;
    Ok(d)
}

fn failure_cause_dict<'py>(
    py: Python<'py>,
    cause: &chematic_3d::pipeline_v2::PipelineV2FailureCause,
) -> PyResult<Bound<'py, PyDict>> {
    use chematic_3d::pipeline_v2::PipelineV2FailureCause;
    let d = PyDict::new(py);
    match cause {
        PipelineV2FailureCause::InvalidConfiguration => {
            d.set_item("kind", "invalid_configuration")?
        }
        PipelineV2FailureCause::BoundAdjustmentFailed => {
            d.set_item("kind", "bound_adjustment_failed")?
        }
        PipelineV2FailureCause::DistanceGeometry(e) => {
            d.set_item("kind", "distance_geometry")?;
            d.set_item("embed_failure_cause", snake_case_debug(e))?;
        }
        PipelineV2FailureCause::TorsionKnowledge(e) => {
            d.set_item("kind", "torsion_knowledge")?;
            d.set_item("torsion_knowledge_error", snake_case_debug(e))?;
        }
        PipelineV2FailureCause::RingTorsionApplicationUnsupported => {
            d.set_item("kind", "ring_torsion_application_unsupported")?
        }
        PipelineV2FailureCause::StereoRepairFailed => d.set_item("kind", "stereo_repair_failed")?,
        PipelineV2FailureCause::StereoUnevaluableUnderStrictPolicy => {
            d.set_item("kind", "stereo_unevaluable_under_strict_policy")?
        }
        PipelineV2FailureCause::ForceField(e) => {
            d.set_item("kind", "force_field")?;
            d.set_item(
                "force_field_bridge_error",
                force_field_bridge_error_dict(py, e)?,
            )?;
        }
        PipelineV2FailureCause::FinalStereoViolation => {
            d.set_item("kind", "final_stereo_violation")?
        }
        PipelineV2FailureCause::FinalGeometryInvalid => {
            d.set_item("kind", "final_geometry_invalid")?
        }
        PipelineV2FailureCause::Timeout => d.set_item("kind", "timeout")?,
    }
    Ok(d)
}

pub(crate) fn failure_to_dict<'py>(
    py: Python<'py>,
    f: &PipelineV2Failure,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("cause", failure_cause_dict(py, &f.cause)?)?;
    d.set_item("stage", snake_case_debug(&f.stage))?;
    // Deliberately named/flagged so it can never be mistaken for a usable
    // result -- see `PipelineV2Failure::last_known_coords`'s own doc comment.
    match &f.last_known_coords {
        Some(c) => d.set_item("last_known_coords", coords_to_vec(c))?,
        None => d.set_item("last_known_coords", py.None())?,
    }
    d.set_item("coords_are_diagnostic_only", true)?;
    match &f.embed_stats {
        Some(s) => d.set_item("embed_stats", embed_stats_dict(py, s)?)?,
        None => d.set_item("embed_stats", py.None())?,
    }
    match &f.bound_adjustment_report {
        Some(v) => d.set_item(
            "bound_adjustment_report",
            v.iter()
                .map(|a| bound_adjustment_dict(py, a))
                .collect::<PyResult<Vec<_>>>()?,
        )?,
        None => d.set_item("bound_adjustment_report", py.None())?,
    }
    match &f.torsion_knowledge_report {
        Some(r) => d.set_item(
            "torsion_knowledge_report",
            torsion_knowledge_report_dict(py, r)?,
        )?,
        None => d.set_item("torsion_knowledge_report", py.None())?,
    }
    match &f.ring_torsion_evidence {
        Some(e) => d.set_item("ring_torsion_evidence", ring_torsion_evidence_dict(py, e)?)?,
        None => d.set_item("ring_torsion_evidence", py.None())?,
    }
    match &f.torsion_optimization_report {
        Some(r) => d.set_item(
            "torsion_optimization_report",
            torsion_optimization_report_dict(py, r)?,
        )?,
        None => d.set_item("torsion_optimization_report", py.None())?,
    }
    match &f.stereo_before {
        Some(s) => d.set_item("stereo_before", stereo_verification_dict(py, s)?)?,
        None => d.set_item("stereo_before", py.None())?,
    }
    match &f.stereo_repair {
        Some(s) => d.set_item("stereo_repair", stereo_repair_summary_dict(py, s)?)?,
        None => d.set_item("stereo_repair", py.None())?,
    }
    match &f.stereo_after_repair {
        Some(s) => d.set_item("stereo_after_repair", stereo_verification_dict(py, s)?)?,
        None => d.set_item("stereo_after_repair", py.None())?,
    }
    match &f.force_field {
        Some(r) => d.set_item("force_field", policy_minimize_result_dict(py, r)?)?,
        None => d.set_item("force_field", py.None())?,
    }
    match &f.final_stereo {
        Some(s) => d.set_item("final_stereo", stereo_verification_dict(py, s)?)?,
        None => d.set_item("final_stereo", py.None())?,
    }
    match &f.final_validation {
        Some(v) => d.set_item("final_validation", final_validation_dict(py, v)?)?,
        None => d.set_item("final_validation", py.None())?,
    }
    d.set_item(
        "elapsed_ms_by_stage",
        stage_timings_dict(py, &f.elapsed_ms_by_stage)?,
    )?;
    Ok(d)
}

// ---------------------------------------------------------------------------
// Entry point called from `Mol::embed_pipeline_v2` in mol_methods.rs
// ---------------------------------------------------------------------------

pub fn run_embed_pipeline_v2<'py>(
    py: Python<'py>,
    mol: &Molecule,
    config: &PyPipelineV2Config,
) -> PyResult<Bound<'py, PyDict>> {
    match pv2::embed_pipeline_v2(mol, &config.inner) {
        Ok(result) => result_to_dict(py, &result),
        Err(failure) => {
            let diagnostics = failure_to_dict(py, &failure)?;
            let message = format!(
                "pipeline_v2 failed at stage {:?}: {:?}",
                failure.stage, failure.cause
            );
            let err = PipelineV2Error::new_err(message);
            err.value(py).setattr("diagnostics", diagnostics)?;
            Err(err)
        }
    }
}
