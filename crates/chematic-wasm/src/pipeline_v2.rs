//! WASM binding for `chematic_3d::pipeline_v2::embed_pipeline_v2`.
//!
//! Mirrors the Python binding (`crates/chematic-py/src/pipeline_v2.rs`, already on
//! `main`) as closely as JS/JSON allow: the same 15 config fields (camelCase JSON
//! keys), the same snake_case policy string values (`"repair_and_verify"`,
//! `"mmff94_with_uff_fallback"`, ...), and the same success/failure evidence fields
//! (camelCase keys, same semantic content). Applies directly to `MolHandle.inner` --
//! never canonicalizes/reparses (the still-open `conformer_ensemble_json` anti-pattern,
//! CHANGELOG-documented as a separate, out-of-scope bug, is deliberately not repeated
//! or touched here).
//!
//! Every JSON output is a tagged union with `schemaVersion: 1`: `{"ok": true,
//! "result": {...}}` or `{"ok": false, "error": {...}}`. `embed_pipeline_v2_json`
//! never throws -- WASM-level safety-limit failures (oversized input, too many
//! atoms, malformed/incomplete config) and real pipeline failures both come back as
//! the same `{"ok": false, ...}` shape, distinguished only by `error.stage` /
//! `error.cause.kind`.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use chematic_3d::distance_geometry_v2::EmbedParameters;
use chematic_3d::etkdg_knowledge::TorsionOptimizationConfig;
use chematic_3d::minimize::ForceFieldPolicy;
use chematic_3d::pipeline_v2 as pv2;

use crate::{MolHandle, WASM_MAX_ATOMS, WASM_MAX_INPUT_BYTES};

const SCHEMA_VERSION: u32 = 1;

/// Pseudo-stage for failures that never reach the real 12-stage pipeline at all
/// (oversized input, too many atoms, malformed/incomplete config JSON) -- kept
/// textually distinct from every real `PipelineStage` variant's own snake_case name.
const STAGE_WASM_INPUT_VALIDATION: &str = "wasm_input_validation";

// ---------------------------------------------------------------------------
// f64 that can never produce invalid JSON
// ---------------------------------------------------------------------------

/// JSON cannot represent NaN/Infinity. Several pipeline fields (eigenvalue
/// magnitudes, energies, residual forces) can go non-finite on degenerate
/// geometry. Finite values serialize at full, round-trippable precision (serde_json's
/// own `f64` writer, never manually rounded); non-finite values serialize as JSON
/// `null` rather than ever emitting a literal `NaN`/`Infinity` token, which would
/// make the output invalid JSON and break `JSON.parse` on the JS side.
#[derive(Debug, Clone, Copy)]
struct FiniteF64(f64);

impl From<f64> for FiniteF64 {
    fn from(v: f64) -> Self {
        FiniteF64(v)
    }
}

impl Serialize for FiniteF64 {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        if self.0.is_finite() {
            serializer.serialize_f64(self.0)
        } else {
            serializer.serialize_none()
        }
    }
}

/// `format!("{value:?}")` -> `snake_case`, for the many small fieldless enums in
/// this pipeline (e.g. `EmbedFailureCause::BoundsSmoothingFailed` ->
/// `"bounds_smoothing_failed"`). Matches the identical helper in the Python binding.
fn snake_case_debug<T: std::fmt::Debug>(value: &T) -> String {
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
// Config: string enum contracts (must match the Python binding's exact values)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
enum StereoPolicyJson {
    #[serde(rename = "ignore")]
    Ignore,
    #[serde(rename = "verify_only")]
    VerifyOnly,
    #[serde(rename = "repair_and_verify")]
    RepairAndVerify,
}

impl From<StereoPolicyJson> for pv2::StereoPolicy {
    fn from(v: StereoPolicyJson) -> Self {
        match v {
            StereoPolicyJson::Ignore => pv2::StereoPolicy::Ignore,
            StereoPolicyJson::VerifyOnly => pv2::StereoPolicy::VerifyOnly,
            StereoPolicyJson::RepairAndVerify => pv2::StereoPolicy::RepairAndVerify,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
enum RingTorsionPolicyJson {
    #[serde(rename = "fail_closed")]
    FailClosed,
    #[serde(rename = "diagnostic_only")]
    DiagnosticOnly,
}

impl From<RingTorsionPolicyJson> for pv2::RingTorsionApplicationPolicy {
    fn from(v: RingTorsionPolicyJson) -> Self {
        match v {
            RingTorsionPolicyJson::FailClosed => pv2::RingTorsionApplicationPolicy::FailClosed,
            RingTorsionPolicyJson::DiagnosticOnly => {
                pv2::RingTorsionApplicationPolicy::DiagnosticOnly
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
enum ForceFieldPolicyJson {
    #[serde(rename = "mmff94_bond_angle_strict")]
    Mmff94BondAngleStrict,
    #[serde(rename = "mmff94_with_uff_fallback")]
    Mmff94WithUffFallback,
    #[serde(rename = "uff_only")]
    UffOnly,
    #[serde(rename = "dreiding")]
    Dreiding,
    #[serde(rename = "none")]
    None,
}

impl From<ForceFieldPolicyJson> for ForceFieldPolicy {
    fn from(v: ForceFieldPolicyJson) -> Self {
        match v {
            ForceFieldPolicyJson::Mmff94BondAngleStrict => ForceFieldPolicy::Mmff94BondAngleStrict,
            ForceFieldPolicyJson::Mmff94WithUffFallback => ForceFieldPolicy::Mmff94WithUffFallback,
            ForceFieldPolicyJson::UffOnly => ForceFieldPolicy::UffOnly,
            ForceFieldPolicyJson::Dreiding => ForceFieldPolicy::Dreiding,
            ForceFieldPolicyJson::None => ForceFieldPolicy::None,
        }
    }
}

fn force_field_policy_str(p: ForceFieldPolicy) -> &'static str {
    match p {
        ForceFieldPolicy::Mmff94BondAngleStrict => "mmff94_bond_angle_strict",
        ForceFieldPolicy::Mmff94WithUffFallback => "mmff94_with_uff_fallback",
        ForceFieldPolicy::UffOnly => "uff_only",
        ForceFieldPolicy::Dreiding => "dreiding",
        ForceFieldPolicy::None => "none",
    }
}

// ---------------------------------------------------------------------------
// Config: input JSON shape
// ---------------------------------------------------------------------------

/// Standard serde idiom for "the key must be present, but its value may be JSON
/// `null`" -- distinct from a bare `Option<T>` field, which serde silently defaults
/// to `None` when the key is *absent* (exactly the silent-default behavior the task
/// requires rejecting for `embedTimeoutMs`/`totalTimeoutMs`, since Python's
/// constructor requires both as explicit, always-present arguments too, even though
/// their value can be `None`).
fn deserialize_present<'de, D, T>(deserializer: D) -> Result<Option<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: serde::Deserialize<'de>,
{
    T::deserialize(deserializer).map(Some)
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
struct PipelineV2ConfigJson {
    embed_seed: u64,
    max_attempts: usize,
    #[serde(deserialize_with = "deserialize_present")]
    embed_timeout_ms: Option<Option<u64>>,
    use_exp_torsions: bool,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
    use_macrocycle_14_bounds: bool,
    include_legacy_torsion_heuristic: bool,
    stereo_policy: StereoPolicyJson,
    fail_on_unevaluable_stereo: bool,
    force_field_policy: ForceFieldPolicyJson,
    force_field_max_iterations: usize,
    gate_mmff94_torsion_oop: bool,
    // #[serde(default)]: added after the WASM binding's 15-field JSON config
    // was already a documented external API (Priority 2, issue #227) --
    // existing callers' configs must keep parsing (as `false`, matching
    // production's own default and every existing arm's unchanged
    // behavior), not start failing deny_unknown_fields' sibling check
    // (a *missing* required field, not an *unknown* one, but the same
    // "never silently break an external caller" principle applies).
    #[serde(default)]
    gate_mmff94_stretch_bend: bool,
    ring_torsion_policy: RingTorsionPolicyJson,
    #[serde(deserialize_with = "deserialize_present")]
    total_timeout_ms: Option<Option<u64>>,
    // #[serde(default)]: same precedent as `gate_mmff94_stretch_bend` above --
    // added after the JSON config was already a documented external API
    // (v0.14.0, issue #285's E/Z bound fix). Existing callers' configs must
    // keep parsing (as `false`, matching `EmbedParameters::default()` and
    // every existing arm's unchanged behavior).
    #[serde(default)]
    enforce_chirality: bool,
    // #[serde(default)]: same precedent again (issue #291/#383). Requires
    // `enforceChirality: true` (raises the same `invalid_configuration` error
    // otherwise) -- see `PipelineV2Config::expand_implicit_h_through_pipeline`'s
    // own Rust doc for what it does. Prefer `pipeline_v2_stereo_safe_config_json`
    // over setting this field alone: it only works correctly combined with
    // `stereoPolicy: "repair_and_verify"` and `enforceChirality: true`, which
    // that helper sets together so a caller can't set one but forget another.
    #[serde(default)]
    expand_implicit_h_through_pipeline: bool,
}

impl PipelineV2ConfigJson {
    fn into_pipeline_config(self) -> Result<pv2::PipelineV2Config, String> {
        let embed_timeout_ms = self.embed_timeout_ms.ok_or_else(|| {
            "missing field `embedTimeoutMs` (must be present; value may be null)".to_string()
        })?;
        let total_timeout_ms = self.total_timeout_ms.ok_or_else(|| {
            "missing field `totalTimeoutMs` (must be present; value may be null)".to_string()
        })?;
        Ok(pv2::PipelineV2Config {
            embed: EmbedParameters {
                random_seed: self.embed_seed,
                max_attempts: self.max_attempts,
                timeout_ms: embed_timeout_ms,
                use_exp_torsions: self.use_exp_torsions,
                use_small_ring_torsions: self.use_small_ring_torsions,
                use_macrocycle_torsions: self.use_macrocycle_torsions,
                use_macrocycle_14_bounds: self.use_macrocycle_14_bounds,
                enforce_chirality: self.enforce_chirality,
                ..EmbedParameters::default()
            },
            torsion_optimization: TorsionOptimizationConfig::default(),
            include_legacy_torsion_heuristic: self.include_legacy_torsion_heuristic,
            stereo_policy: self.stereo_policy.into(),
            fail_on_unevaluable_stereo: self.fail_on_unevaluable_stereo,
            force_field_policy: self.force_field_policy.into(),
            force_field_max_iterations: self.force_field_max_iterations,
            gate_mmff94_torsion_oop: self.gate_mmff94_torsion_oop,
            gate_mmff94_stretch_bend: self.gate_mmff94_stretch_bend,
            ring_torsion_policy: self.ring_torsion_policy.into(),
            total_timeout_ms,
            expand_implicit_h_through_pipeline: self.expand_implicit_h_through_pipeline,
        })
    }
}

/// Parses+validates `config_json` in one step. Never silently defaults an unknown
/// field, a missing required field (including the two nullable timeouts, which must
/// still be explicitly present), an unknown enum string, or an out-of-range/wrong-type
/// integer -- all of these are `serde`'s own `deny_unknown_fields`/required-field/
/// closed-enum/typed-integer behavior, not hand-rolled validation.
fn parse_pipeline_config(config_json: &str) -> Result<pv2::PipelineV2Config, String> {
    let parsed: PipelineV2ConfigJson =
        serde_json::from_str(config_json).map_err(|e| e.to_string())?;
    parsed.into_pipeline_config()
}

// ---------------------------------------------------------------------------
// Result: nested report types (camelCase keys, semantically identical to the
// Python binding's dict conversion in crates/chematic-py/src/pipeline_v2.rs)
// ---------------------------------------------------------------------------

fn coords_to_json(coords: &chematic_3d::coords::Coords3D) -> Vec<[FiniteF64; 3]> {
    coords
        .points
        .iter()
        .map(|p| [p.x.into(), p.y.into(), p.z.into()])
        .collect()
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct EmbedStatsJson {
    attempts_used: usize,
    failure_counts: std::collections::HashMap<String, usize>,
    negative_eigenvalues_beyond_embedding_dim: usize,
    max_negative_eigenvalue_magnitude: FiniteF64,
    last_smoothing_invariants_ok: bool,
    used_random_coords: bool,
    adjustments_applied: usize,
}

fn embed_stats_json(stats: &chematic_3d::distance_geometry_v2::EmbedStats) -> EmbedStatsJson {
    EmbedStatsJson {
        attempts_used: stats.attempts_used,
        failure_counts: stats
            .failure_counts
            .iter()
            .map(|(cause, count)| (snake_case_debug(cause), *count))
            .collect(),
        negative_eigenvalues_beyond_embedding_dim: stats.negative_eigenvalues_beyond_embedding_dim,
        max_negative_eigenvalue_magnitude: stats.max_negative_eigenvalue_magnitude.into(),
        last_smoothing_invariants_ok: stats.last_smoothing_invariants_ok,
        used_random_coords: stats.used_random_coords,
        adjustments_applied: stats.adjustments_applied,
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct BoundAdjustmentJson {
    atom1: u32,
    atom2: u32,
    old_lower: FiniteF64,
    new_lower: FiniteF64,
    old_upper: FiniteF64,
    new_upper: FiniteF64,
    rule_id: String,
    source: String,
    ring_size: usize,
    reason: String,
}

fn bound_adjustment_json(
    a: &chematic_3d::etkdg_knowledge::PairBoundAdjustment,
) -> BoundAdjustmentJson {
    BoundAdjustmentJson {
        atom1: a.atom_pair.0.0,
        atom2: a.atom_pair.1.0,
        old_lower: a.old_lower.into(),
        new_lower: a.new_lower.into(),
        old_upper: a.old_upper.into(),
        new_upper: a.new_upper.into(),
        rule_id: a.rule_id.clone(),
        source: snake_case_debug(&a.source),
        ring_size: a.ring_size,
        reason: a.reason.clone(),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct FourierTermJson {
    periodicity: u8,
    phase_deg: FiniteF64,
    amplitude: FiniteF64,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct TorsionPotentialJson {
    atoms: [u32; 4],
    central_bond: (u32, u32),
    source: String,
    rule_id: String,
    terms: Vec<FourierTermJson>,
    ring_size: Option<usize>,
}

fn torsion_potential_json(
    p: &chematic_3d::etkdg_knowledge::TorsionPotential,
) -> TorsionPotentialJson {
    TorsionPotentialJson {
        atoms: [p.atoms[0].0, p.atoms[1].0, p.atoms[2].0, p.atoms[3].0],
        central_bond: (p.central_bond.0.0, p.central_bond.1.0),
        source: snake_case_debug(&p.source),
        rule_id: p.rule_id.clone(),
        terms: p
            .terms
            .iter()
            .map(|t| FourierTermJson {
                periodicity: t.periodicity,
                phase_deg: t.phase_deg.into(),
                amplitude: t.amplitude.into(),
            })
            .collect(),
        ring_size: p.ring_size,
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct TorsionDiagnosticJson {
    central_bond: (u32, u32),
    kind: String,
    message: String,
    candidate_rule_ids: Vec<String>,
}

fn torsion_diagnostic_json(
    d: &chematic_3d::etkdg_knowledge::TorsionKnowledgeDiagnostic,
) -> TorsionDiagnosticJson {
    TorsionDiagnosticJson {
        central_bond: (d.central_bond.0.0, d.central_bond.1.0),
        kind: snake_case_debug(&d.kind),
        message: d.message.clone(),
        candidate_rule_ids: d.candidate_rule_ids.clone(),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct TorsionKnowledgeReportJson {
    potentials: Vec<TorsionPotentialJson>,
    matched_rule_ids: Vec<String>,
    unmatched_rotatable_bonds: Vec<(u32, u32)>,
    ambiguous_matches: Vec<TorsionDiagnosticJson>,
    skipped_bonds: Vec<TorsionDiagnosticJson>,
}

fn torsion_knowledge_report_json(
    r: &chematic_3d::etkdg_knowledge::TorsionKnowledgeReport,
) -> TorsionKnowledgeReportJson {
    TorsionKnowledgeReportJson {
        potentials: r.potentials.iter().map(torsion_potential_json).collect(),
        matched_rule_ids: r.matched_rule_ids.clone(),
        unmatched_rotatable_bonds: r
            .unmatched_rotatable_bonds
            .iter()
            .map(|(a, b)| (a.0, b.0))
            .collect(),
        ambiguous_matches: r
            .ambiguous_matches
            .iter()
            .map(torsion_diagnostic_json)
            .collect(),
        skipped_bonds: r
            .skipped_bonds
            .iter()
            .map(torsion_diagnostic_json)
            .collect(),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct RingTorsionPotentialEvidenceJson {
    rule_id: String,
    central_bond: (u32, u32),
    source: String,
    applied_to_geometry: bool,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct RingTorsionEvidenceJson {
    potentials: Vec<RingTorsionPotentialEvidenceJson>,
    diagnostic_only: bool,
    n_applied: usize,
    n_scored_only: usize,
}

fn ring_torsion_evidence_json(e: &pv2::RingTorsionEvidence) -> RingTorsionEvidenceJson {
    RingTorsionEvidenceJson {
        potentials: e
            .potentials
            .iter()
            .map(|p| RingTorsionPotentialEvidenceJson {
                rule_id: p.rule_id.clone(),
                central_bond: (p.central_bond.0.0, p.central_bond.1.0),
                source: snake_case_debug(&p.source),
                applied_to_geometry: p.applied_to_geometry,
            })
            .collect(),
        diagnostic_only: e.diagnostic_only,
        n_applied: e.n_applied(),
        n_scored_only: e.n_scored_only(),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct TorsionOptimizationReportJson {
    energy_before: FiniteF64,
    energy_after: FiniteF64,
    iterations_used: usize,
    converged: bool,
    max_bond_length_delta: FiniteF64,
    max_ring_closure_delta: FiniteF64,
    rotated_bond_count: usize,
}

fn torsion_optimization_report_json(
    r: &chematic_3d::etkdg_knowledge::TorsionOptimizationReport,
) -> TorsionOptimizationReportJson {
    TorsionOptimizationReportJson {
        energy_before: r.energy_before.into(),
        energy_after: r.energy_after.into(),
        iterations_used: r.iterations_used,
        converged: r.converged,
        max_bond_length_delta: r.max_bond_length_delta.into(),
        max_ring_closure_delta: r.max_ring_closure_delta.into(),
        rotated_bond_count: r.rotated_bond_count,
    }
}

fn stereo_status_fields(
    status: chematic_3d::stereo_constraints::StereoStatus,
) -> (String, Option<String>) {
    use chematic_3d::stereo_constraints::StereoStatus;
    match status {
        StereoStatus::Satisfied => ("satisfied".to_string(), None),
        StereoStatus::Violated => ("violated".to_string(), None),
        StereoStatus::Unevaluable(reason) => {
            ("unevaluable".to_string(), Some(snake_case_debug(&reason)))
        }
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct TetrahedralReportJson {
    atom: u32,
    status: String,
    rejection_reason: Option<String>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct DoubleBondReportJson {
    bond: u32,
    status: String,
    rejection_reason: Option<String>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct StereoVerificationJson {
    tetrahedral: Vec<TetrahedralReportJson>,
    double_bond: Vec<DoubleBondReportJson>,
    n_declared: usize,
    n_satisfied: usize,
    n_violations: usize,
    n_unevaluable: usize,
    is_fully_satisfied: bool,
}

fn stereo_verification_json(
    v: &chematic_3d::stereo_constraints::StereoVerification,
) -> StereoVerificationJson {
    StereoVerificationJson {
        tetrahedral: v
            .tetrahedral
            .iter()
            .map(|r| {
                let (status, rejection_reason) = stereo_status_fields(r.status);
                TetrahedralReportJson {
                    atom: r.atom.0,
                    status,
                    rejection_reason,
                }
            })
            .collect(),
        double_bond: v
            .double_bond
            .iter()
            .map(|r| {
                let (status, rejection_reason) = stereo_status_fields(r.status);
                DoubleBondReportJson {
                    bond: r.bond.0,
                    status,
                    rejection_reason,
                }
            })
            .collect(),
        n_declared: v.n_declared(),
        n_satisfied: v.n_satisfied(),
        n_violations: v.n_violations(),
        n_unevaluable: v.n_unevaluable(),
        is_fully_satisfied: v.is_fully_satisfied(),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct RepairedElementJson {
    element_kind: String,
    atom: Option<u32>,
    bond: Option<u32>,
    atoms_moved: usize,
    max_displacement: FiniteF64,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct RepairFailureJson {
    element_kind: String,
    atom: Option<u32>,
    bond: Option<u32>,
    reason: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct StereoRepairSummaryJson {
    repaired: Vec<RepairedElementJson>,
    failures: Vec<RepairFailureJson>,
}

fn stereo_element_fields(
    element: chematic_3d::stereo_constraints::StereoElement,
) -> (&'static str, Option<u32>, Option<u32>) {
    use chematic_3d::stereo_constraints::StereoElement;
    match element {
        StereoElement::Tetrahedral(atom) => ("tetrahedral", Some(atom.0), None),
        StereoElement::DoubleBond(bond) => ("double_bond", None, Some(bond.0)),
    }
}

fn stereo_repair_summary_json(s: &pv2::StereoRepairSummary) -> StereoRepairSummaryJson {
    StereoRepairSummaryJson {
        repaired: s
            .repaired
            .iter()
            .map(|r| {
                let (element_kind, atom, bond) = stereo_element_fields(r.element);
                RepairedElementJson {
                    element_kind: element_kind.to_string(),
                    atom,
                    bond,
                    atoms_moved: r.atoms_moved,
                    max_displacement: r.max_displacement.into(),
                }
            })
            .collect(),
        failures: s
            .failures
            .iter()
            .map(|(element, reason)| {
                let (element_kind, atom, bond) = stereo_element_fields(*element);
                RepairFailureJson {
                    element_kind: element_kind.to_string(),
                    atom,
                    bond,
                    reason: snake_case_debug(reason),
                }
            })
            .collect(),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct Mmff94MissingTermJson {
    kind: String,
    atoms: Vec<u32>,
    description: String,
}

fn mmff94_missing_term_json(t: &chematic_3d::minimize::Mmff94MissingTerm) -> Mmff94MissingTermJson {
    Mmff94MissingTermJson {
        kind: snake_case_debug(&t.kind),
        atoms: t.atoms.iter().map(|a| a.0).collect(),
        description: t.description.clone(),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct Mmff94CoverageJson {
    bonds_total: usize,
    bonds_missing: Vec<Mmff94MissingTermJson>,
    angles_total: usize,
    angles_missing: Vec<Mmff94MissingTermJson>,
    torsions_total: usize,
    torsions_missing: Vec<Mmff94MissingTermJson>,
    oop_total: usize,
    oop_missing: Vec<Mmff94MissingTermJson>,
    stretch_bend_total: usize,
    stretch_bend_missing: Vec<Mmff94MissingTermJson>,
}

fn mmff94_coverage_json(r: &chematic_3d::minimize::Mmff94CoverageReport) -> Mmff94CoverageJson {
    Mmff94CoverageJson {
        bonds_total: r.bonds_total,
        bonds_missing: r
            .bonds_missing
            .iter()
            .map(mmff94_missing_term_json)
            .collect(),
        angles_total: r.angles_total,
        angles_missing: r
            .angles_missing
            .iter()
            .map(mmff94_missing_term_json)
            .collect(),
        torsions_total: r.torsions_total,
        torsions_missing: r
            .torsions_missing
            .iter()
            .map(mmff94_missing_term_json)
            .collect(),
        oop_total: r.oop_total,
        oop_missing: r.oop_missing.iter().map(mmff94_missing_term_json).collect(),
        stretch_bend_total: r.stretch_bend_total,
        stretch_bend_missing: r
            .stretch_bend_missing
            .iter()
            .map(mmff94_missing_term_json)
            .collect(),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ForceFieldBridgeErrorJson {
    kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    message: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    coverage: Option<Mmff94CoverageJson>,
    #[serde(skip_serializing_if = "Option::is_none")]
    policy: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    reason: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    converged: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    iterations: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    max_residual_force: Option<FiniteF64>,
}

fn force_field_bridge_error_json(
    e: &chematic_3d::minimize::ForceFieldBridgeError,
) -> ForceFieldBridgeErrorJson {
    use chematic_3d::minimize::ForceFieldBridgeError;
    match e {
        ForceFieldBridgeError::UnsupportedAtomType(msg) => ForceFieldBridgeErrorJson {
            kind: "unsupported_atom_type".to_string(),
            message: Some(msg.clone()),
            coverage: None,
            policy: None,
            reason: None,
            converged: None,
            iterations: None,
            max_residual_force: None,
        },
        ForceFieldBridgeError::MissingParameters(coverage) => ForceFieldBridgeErrorJson {
            kind: "missing_parameters".to_string(),
            message: None,
            coverage: Some(mmff94_coverage_json(coverage)),
            policy: None,
            reason: None,
            converged: None,
            iterations: None,
            max_residual_force: None,
        },
        ForceFieldBridgeError::MinimizationFailed(detail) => ForceFieldBridgeErrorJson {
            kind: "minimization_failed".to_string(),
            message: None,
            coverage: None,
            policy: Some(force_field_policy_str(detail.policy).to_string()),
            reason: Some(snake_case_debug(&detail.reason)),
            converged: Some(detail.converged),
            iterations: Some(detail.iterations),
            max_residual_force: Some(detail.max_residual_force.into()),
        },
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct EnergyReportJson {
    kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    bond: Option<FiniteF64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    angle: Option<FiniteF64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    stretch_bend: Option<FiniteF64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    torsion: Option<FiniteF64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    oop: Option<FiniteF64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    vdw: Option<FiniteF64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    electrostatic: Option<FiniteF64>,
    total: FiniteF64,
}

fn energy_report_json(e: &chematic_3d::minimize::EnergyReport) -> EnergyReportJson {
    use chematic_3d::minimize::EnergyReport;
    match e {
        EnergyReport::Mmff94(b) => EnergyReportJson {
            kind: "mmff94".to_string(),
            bond: Some(b.bond.into()),
            angle: Some(b.angle.into()),
            stretch_bend: Some(b.stretch_bend.into()),
            torsion: Some(b.torsion.into()),
            oop: Some(b.oop.into()),
            vdw: Some(b.vdw.into()),
            electrostatic: Some(b.electrostatic.into()),
            total: b.total.into(),
        },
        EnergyReport::Uff { total } => EnergyReportJson {
            kind: "uff".to_string(),
            bond: None,
            angle: None,
            stretch_bend: None,
            torsion: None,
            oop: None,
            vdw: None,
            electrostatic: None,
            total: (*total).into(),
        },
        EnergyReport::Dreiding { total } => EnergyReportJson {
            kind: "dreiding".to_string(),
            bond: None,
            angle: None,
            stretch_bend: None,
            torsion: None,
            oop: None,
            vdw: None,
            electrostatic: None,
            total: (*total).into(),
        },
        EnergyReport::None => EnergyReportJson {
            kind: "none".to_string(),
            bond: None,
            angle: None,
            stretch_bend: None,
            torsion: None,
            oop: None,
            vdw: None,
            electrostatic: None,
            total: 0.0.into(),
        },
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PolicyMinimizeResultJson {
    coords: Vec<[FiniteF64; 3]>,
    requested_force_field: String,
    actual_force_field_used: String,
    fallback_reason: Option<ForceFieldBridgeErrorJson>,
    missing_parameter_classes: Vec<Mmff94MissingTermJson>,
    coverage: Option<Mmff94CoverageJson>,
    energy_before: EnergyReportJson,
    energy_after: EnergyReportJson,
    converged: bool,
    iterations: usize,
    max_residual_force: FiniteF64,
    starting_geometry: Option<String>,
}

fn policy_minimize_result_json(
    r: &chematic_3d::minimize::PolicyMinimizeResult,
) -> PolicyMinimizeResultJson {
    PolicyMinimizeResultJson {
        coords: coords_to_json(&r.coords),
        requested_force_field: force_field_policy_str(r.requested_force_field).to_string(),
        actual_force_field_used: force_field_policy_str(r.actual_force_field_used).to_string(),
        fallback_reason: r
            .fallback_reason
            .as_ref()
            .map(force_field_bridge_error_json),
        missing_parameter_classes: r
            .missing_parameter_classes
            .iter()
            .map(mmff94_missing_term_json)
            .collect(),
        coverage: r.coverage.as_ref().map(mmff94_coverage_json),
        energy_before: energy_report_json(&r.energy_before),
        energy_after: energy_report_json(&r.energy_after),
        converged: r.converged,
        iterations: r.iterations,
        max_residual_force: r.max_residual_force.into(),
        starting_geometry: r.starting_geometry.map(|g| snake_case_debug(&g)),
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct BoundsConformanceJson {
    n_pairs: usize,
    n_violations: usize,
    max_rel_violation: FiniteF64,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct FinalValidationJson {
    all_finite: bool,
    atom_count_unchanged: bool,
    worst_bond_length: FiniteF64,
    bond_violation_rate_15pct: FiniteF64,
    bond_violation_rate_50pct: FiniteF64,
    gross_clash_count: usize,
    bounds_conformance: BoundsConformanceJson,
    stereo_ok: bool,
    torsion_energy_after: FiniteF64,
    ring_closure_delta: FiniteF64,
    sound: bool,
}

fn final_validation_json(v: &pv2::FinalGeometryValidation) -> FinalValidationJson {
    FinalValidationJson {
        all_finite: v.all_finite,
        atom_count_unchanged: v.atom_count_unchanged,
        worst_bond_length: v.worst_bond_length.into(),
        bond_violation_rate_15pct: v.bond_violation_rate_15pct.into(),
        bond_violation_rate_50pct: v.bond_violation_rate_50pct.into(),
        gross_clash_count: v.gross_clash_count,
        bounds_conformance: BoundsConformanceJson {
            n_pairs: v.bounds_conformance.n_pairs,
            n_violations: v.bounds_conformance.n_violations,
            max_rel_violation: v.bounds_conformance.max_rel_violation.into(),
        },
        stereo_ok: v.stereo_ok,
        torsion_energy_after: v.torsion_energy_after.into(),
        ring_closure_delta: v.ring_closure_delta.into(),
        sound: v.sound,
    }
}

#[derive(Serialize, Default)]
#[serde(rename_all = "camelCase")]
struct StageTimingsJson {
    torsion_knowledge_ms: u64,
    bound_adjustment_ms: u64,
    distance_geometry_ms: u64,
    torsion_energy_eval_ms: u64,
    torsion_optimization_ms: u64,
    stereo_verify_before_ms: u64,
    stereo_repair_ms: u64,
    stereo_verify_after_repair_ms: u64,
    force_field_ms: u64,
    final_stereo_verify_ms: u64,
    final_validation_ms: u64,
    total_ms: u64,
}

fn stage_timings_json(t: &pv2::StageTimings) -> StageTimingsJson {
    StageTimingsJson {
        torsion_knowledge_ms: t.torsion_knowledge_ms,
        bound_adjustment_ms: t.bound_adjustment_ms,
        distance_geometry_ms: t.distance_geometry_ms,
        torsion_energy_eval_ms: t.torsion_energy_eval_ms,
        torsion_optimization_ms: t.torsion_optimization_ms,
        stereo_verify_before_ms: t.stereo_verify_before_ms,
        stereo_repair_ms: t.stereo_repair_ms,
        stereo_verify_after_repair_ms: t.stereo_verify_after_repair_ms,
        force_field_ms: t.force_field_ms,
        final_stereo_verify_ms: t.final_stereo_verify_ms,
        final_validation_ms: t.final_validation_ms,
        total_ms: t.total_ms,
    }
}

// ---------------------------------------------------------------------------
// Top-level success / failure envelopes
// ---------------------------------------------------------------------------

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PipelineV2SuccessJson {
    coords: Vec<[FiniteF64; 3]>,
    embed_stats: EmbedStatsJson,
    bound_adjustment_report: Option<Vec<BoundAdjustmentJson>>,
    torsion_knowledge_report: TorsionKnowledgeReportJson,
    ring_torsion_evidence: RingTorsionEvidenceJson,
    torsion_optimization_report: Option<TorsionOptimizationReportJson>,
    stereo_before: StereoVerificationJson,
    stereo_repair: Option<StereoRepairSummaryJson>,
    stereo_after_repair: StereoVerificationJson,
    force_field: PolicyMinimizeResultJson,
    final_stereo: StereoVerificationJson,
    final_validation: FinalValidationJson,
    elapsed_ms_by_stage: StageTimingsJson,
}

#[derive(Serialize)]
struct SuccessEnvelopeJson {
    #[serde(rename = "schemaVersion")]
    schema_version: u32,
    ok: bool,
    result: PipelineV2SuccessJson,
}

fn result_to_json(r: &pv2::PipelineV2Result) -> String {
    let envelope = SuccessEnvelopeJson {
        schema_version: SCHEMA_VERSION,
        ok: true,
        result: PipelineV2SuccessJson {
            coords: coords_to_json(&r.coords),
            embed_stats: embed_stats_json(&r.embed_stats),
            bound_adjustment_report: r
                .bound_adjustment_report
                .as_ref()
                .map(|v| v.iter().map(bound_adjustment_json).collect()),
            torsion_knowledge_report: torsion_knowledge_report_json(&r.torsion_knowledge_report),
            ring_torsion_evidence: ring_torsion_evidence_json(&r.ring_torsion_evidence),
            torsion_optimization_report: r
                .torsion_optimization_report
                .as_ref()
                .map(torsion_optimization_report_json),
            stereo_before: stereo_verification_json(&r.stereo_before),
            stereo_repair: r.stereo_repair.as_ref().map(stereo_repair_summary_json),
            stereo_after_repair: stereo_verification_json(&r.stereo_after_repair),
            force_field: policy_minimize_result_json(&r.force_field),
            final_stereo: stereo_verification_json(&r.final_stereo),
            final_validation: final_validation_json(&r.final_validation),
            elapsed_ms_by_stage: stage_timings_json(&r.elapsed_ms_by_stage),
        },
    };
    // `SuccessEnvelopeJson` contains only `FiniteF64`/plain-typed fields built from
    // a real, successful pipeline run -- serialization cannot fail.
    serde_json::to_string(&envelope).expect("success envelope must always serialize")
}

#[derive(Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum FailureCauseJson {
    InvalidConfiguration,
    BoundAdjustmentFailed,
    DistanceGeometry {
        #[serde(rename = "embedFailureCause")]
        embed_failure_cause: String,
    },
    TorsionKnowledge {
        #[serde(rename = "torsionKnowledgeError")]
        torsion_knowledge_error: String,
    },
    RingTorsionApplicationUnsupported,
    StereoRepairFailed,
    StereoUnevaluableUnderStrictPolicy,
    ForceField {
        #[serde(rename = "forceFieldBridgeError")]
        force_field_bridge_error: Box<ForceFieldBridgeErrorJson>,
    },
    FinalStereoViolation,
    FinalGeometryInvalid,
    Timeout,
    /// WASM-level, pre-pipeline: the config JSON itself was malformed, missing a
    /// required field, or used an unknown enum value.
    InvalidConfig {
        message: String,
    },
    /// WASM-level: `mol.atomCount() > WASM_MAX_ATOMS`.
    AtomLimitExceeded {
        limit: usize,
        actual: usize,
    },
    /// WASM-level: `config_json.len() > WASM_MAX_INPUT_BYTES`.
    InputTooLarge {
        #[serde(rename = "limitBytes")]
        limit_bytes: usize,
        #[serde(rename = "actualBytes")]
        actual_bytes: usize,
    },
}

fn failure_cause_json(cause: &pv2::PipelineV2FailureCause) -> FailureCauseJson {
    use pv2::PipelineV2FailureCause;
    match cause {
        PipelineV2FailureCause::InvalidConfiguration => FailureCauseJson::InvalidConfiguration,
        PipelineV2FailureCause::BoundAdjustmentFailed => FailureCauseJson::BoundAdjustmentFailed,
        PipelineV2FailureCause::DistanceGeometry(e) => FailureCauseJson::DistanceGeometry {
            embed_failure_cause: snake_case_debug(e),
        },
        PipelineV2FailureCause::TorsionKnowledge(e) => FailureCauseJson::TorsionKnowledge {
            torsion_knowledge_error: snake_case_debug(e),
        },
        PipelineV2FailureCause::RingTorsionApplicationUnsupported => {
            FailureCauseJson::RingTorsionApplicationUnsupported
        }
        PipelineV2FailureCause::StereoRepairFailed => FailureCauseJson::StereoRepairFailed,
        PipelineV2FailureCause::StereoUnevaluableUnderStrictPolicy => {
            FailureCauseJson::StereoUnevaluableUnderStrictPolicy
        }
        PipelineV2FailureCause::ForceField(e) => FailureCauseJson::ForceField {
            force_field_bridge_error: Box::new(force_field_bridge_error_json(e)),
        },
        PipelineV2FailureCause::FinalStereoViolation => FailureCauseJson::FinalStereoViolation,
        PipelineV2FailureCause::FinalGeometryInvalid => FailureCauseJson::FinalGeometryInvalid,
        PipelineV2FailureCause::Timeout => FailureCauseJson::Timeout,
    }
}

#[derive(Serialize, Default)]
#[serde(rename_all = "camelCase")]
struct FailureDiagnosticsJson {
    embed_stats: Option<EmbedStatsJson>,
    bound_adjustment_report: Option<Vec<BoundAdjustmentJson>>,
    torsion_knowledge_report: Option<TorsionKnowledgeReportJson>,
    ring_torsion_evidence: Option<RingTorsionEvidenceJson>,
    torsion_optimization_report: Option<TorsionOptimizationReportJson>,
    stereo_before: Option<StereoVerificationJson>,
    stereo_repair: Option<StereoRepairSummaryJson>,
    stereo_after_repair: Option<StereoVerificationJson>,
    force_field: Option<PolicyMinimizeResultJson>,
    final_stereo: Option<StereoVerificationJson>,
    final_validation: Option<FinalValidationJson>,
    elapsed_ms_by_stage: StageTimingsJson,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ErrorEnvelopeJson {
    stage: String,
    cause: FailureCauseJson,
    last_known_coords: Option<Vec<[FiniteF64; 3]>>,
    /// Always `true`. Named/typed distinctly from a success `coords` field so a
    /// caller can never mistake this for a usable result -- matches
    /// `PipelineV2Failure::last_known_coords`'s own "diagnostic only" contract on
    /// the Rust/Python side.
    coords_are_diagnostic_only: bool,
    diagnostics: FailureDiagnosticsJson,
}

#[derive(Serialize)]
struct FailureEnvelopeJson {
    #[serde(rename = "schemaVersion")]
    schema_version: u32,
    ok: bool,
    error: ErrorEnvelopeJson,
}

fn failure_to_json(f: &pv2::PipelineV2Failure) -> String {
    let envelope = FailureEnvelopeJson {
        schema_version: SCHEMA_VERSION,
        ok: false,
        error: ErrorEnvelopeJson {
            stage: snake_case_debug(&f.stage),
            cause: failure_cause_json(&f.cause),
            last_known_coords: f.last_known_coords.as_ref().map(coords_to_json),
            coords_are_diagnostic_only: true,
            diagnostics: FailureDiagnosticsJson {
                embed_stats: f.embed_stats.as_ref().map(embed_stats_json),
                bound_adjustment_report: f
                    .bound_adjustment_report
                    .as_ref()
                    .map(|v| v.iter().map(bound_adjustment_json).collect()),
                torsion_knowledge_report: f
                    .torsion_knowledge_report
                    .as_ref()
                    .map(torsion_knowledge_report_json),
                ring_torsion_evidence: f
                    .ring_torsion_evidence
                    .as_ref()
                    .map(ring_torsion_evidence_json),
                torsion_optimization_report: f
                    .torsion_optimization_report
                    .as_ref()
                    .map(torsion_optimization_report_json),
                stereo_before: f.stereo_before.as_ref().map(stereo_verification_json),
                stereo_repair: f.stereo_repair.as_ref().map(stereo_repair_summary_json),
                stereo_after_repair: f.stereo_after_repair.as_ref().map(stereo_verification_json),
                force_field: f.force_field.as_ref().map(policy_minimize_result_json),
                final_stereo: f.final_stereo.as_ref().map(stereo_verification_json),
                final_validation: f.final_validation.as_ref().map(final_validation_json),
                elapsed_ms_by_stage: stage_timings_json(&f.elapsed_ms_by_stage),
            },
        },
    };
    serde_json::to_string(&envelope).expect("failure envelope must always serialize")
}

/// A WASM-level failure that never reached the real pipeline at all (oversized
/// input, too many atoms, malformed/incomplete config JSON). Same envelope shape
/// as a real `PipelineV2Failure`, with `diagnostics` entirely empty (nothing was
/// computed) and `elapsedMsByStage` all-zero.
fn wasm_input_error_json(cause: FailureCauseJson) -> String {
    let envelope = FailureEnvelopeJson {
        schema_version: SCHEMA_VERSION,
        ok: false,
        error: ErrorEnvelopeJson {
            stage: STAGE_WASM_INPUT_VALIDATION.to_string(),
            cause,
            last_known_coords: None,
            coords_are_diagnostic_only: true,
            diagnostics: FailureDiagnosticsJson::default(),
        },
    };
    serde_json::to_string(&envelope).expect("failure envelope must always serialize")
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Run the opt-in v2 embedding pipeline, applied directly to `mol`'s own atom
/// order (never canonicalizes/reparses -- see the module doc).
///
/// `config_json` must be an object with the 15 required fields
/// `PipelineV2Config` requires (camelCase keys: `embedSeed`, `maxAttempts`,
/// `embedTimeoutMs`, `useExpTorsions`, `useSmallRingTorsions`,
/// `useMacrocycleTorsions`, `useMacrocycle14Bounds`,
/// `includeLegacyTorsionHeuristic`, `stereoPolicy`, `failOnUnevaluableStereo`,
/// `forceFieldPolicy`, `forceFieldMaxIterations`, `gateMmff94TorsionOop`,
/// `ringTorsionPolicy`, `totalTimeoutMs`), plus one optional field added in
/// Priority 2 (issue #227): `gateMmff94StretchBend` (`#[serde(default)]` ->
/// `false` if omitted, so pre-Priority-2 caller configs keep working
/// unmodified, matching `false`'s meaning of "existing/unchanged behavior"
/// everywhere else in this codebase). An unknown field, a missing *required*
/// field, an unknown `stereoPolicy`/`ringTorsionPolicy`/`forceFieldPolicy`
/// string, or a wrong-typed/out-of-range integer all fail closed rather than
/// silently defaulting.
///
/// Never throws. Always returns a JSON string tagged with `schemaVersion: 1` and
/// `ok: true`/`false` -- see the module doc for both shapes.
#[wasm_bindgen]
pub fn embed_pipeline_v2_json(mol: &MolHandle, config_json: &str) -> String {
    if config_json.len() > WASM_MAX_INPUT_BYTES {
        return wasm_input_error_json(FailureCauseJson::InputTooLarge {
            limit_bytes: WASM_MAX_INPUT_BYTES,
            actual_bytes: config_json.len(),
        });
    }
    if mol.inner.atom_count() > WASM_MAX_ATOMS {
        return wasm_input_error_json(FailureCauseJson::AtomLimitExceeded {
            limit: WASM_MAX_ATOMS,
            actual: mol.inner.atom_count(),
        });
    }

    let config = match parse_pipeline_config(config_json) {
        Ok(c) => c,
        Err(message) => {
            return wasm_input_error_json(FailureCauseJson::InvalidConfig { message });
        }
    };

    match pv2::embed_pipeline_v2(&mol.inner, &config) {
        Ok(result) => result_to_json(&result),
        Err(failure) => failure_to_json(&failure),
    }
}

/// Returns a ready-to-use `embed_pipeline_v2_json` config JSON string for the
/// "stereo-safe" configuration (issue #291/#383): `stereoPolicy:
/// "repair_and_verify"`, `enforceChirality: true`, and
/// `expandImplicitHThroughPipeline: true` together -- the exact combination
/// measured to correctly handle ring-fused declared stereocenters (e.g.
/// testosterone, cholesterol) that `enforceChirality` alone cannot repair.
/// Mirrors `PipelineV2Config::stereo_safe`/the Python binding's
/// `PipelineV2Config.stereo_safe(...)` -- prefer this over setting those three
/// fields individually: they only work correctly as a set, and forgetting one
/// silently falls back to a configuration issue #291 measured as unsound for
/// that molecule class. `forceFieldPolicy`/`ringTorsionPolicy` are still
/// required, explicit arguments; everything else takes the same conservative
/// defaults `embed_pipeline_v2_json`'s own documented examples do. The caller
/// may parse and further override individual fields before passing the result
/// to `embed_pipeline_v2_json` (e.g. a different `embedSeed`).
///
/// Never throws, matching `embed_pipeline_v2_json`'s own convention: an
/// unknown `force_field`/`ring_torsion_policy` string returns the same
/// `{"ok": false, "error": {...}}` shape `embed_pipeline_v2_json` would for an
/// invalid config, tagged `schemaVersion: 1`.
#[wasm_bindgen]
pub fn pipeline_v2_stereo_safe_config_json(force_field: &str, ring_torsion_policy: &str) -> String {
    // Validate against the same closed enums `embed_pipeline_v2_json` itself
    // parses with -- reuses their exact error messages, no hand-rolled match.
    // The caller's own strings are re-serialized verbatim below (already
    // proven to round-trip, since these enums serialize back to exactly the
    // strings they deserialize from).
    if let Err(e) = serde_json::from_value::<ForceFieldPolicyJson>(serde_json::Value::String(
        force_field.to_string(),
    )) {
        return wasm_input_error_json(FailureCauseJson::InvalidConfig {
            message: format!("forceFieldPolicy: {e}"),
        });
    }
    if let Err(e) = serde_json::from_value::<RingTorsionPolicyJson>(serde_json::Value::String(
        ring_torsion_policy.to_string(),
    )) {
        return wasm_input_error_json(FailureCauseJson::InvalidConfig {
            message: format!("ringTorsionPolicy: {e}"),
        });
    }

    serde_json::json!({
        "schemaVersion": SCHEMA_VERSION,
        "ok": true,
        "config": {
            "embedSeed": 0xC0FF_EE42_D157_6E02_u64,
            "maxAttempts": 8,
            "embedTimeoutMs": null,
            "useExpTorsions": false,
            "useSmallRingTorsions": false,
            "useMacrocycleTorsions": false,
            "useMacrocycle14Bounds": false,
            "includeLegacyTorsionHeuristic": false,
            "stereoPolicy": "repair_and_verify",
            "failOnUnevaluableStereo": false,
            "forceFieldPolicy": force_field,
            "forceFieldMaxIterations": 200,
            "gateMmff94TorsionOop": false,
            "gateMmff94StretchBend": false,
            "ringTorsionPolicy": ring_torsion_policy,
            "totalTimeoutMs": null,
            "enforceChirality": true,
            "expandImplicitHThroughPipeline": true,
        }
    })
    .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_depict::parse_smiles;

    fn safe_config_json(
        force_field: &str,
        stereo_policy: &str,
        ring_torsion_policy: &str,
    ) -> String {
        format!(
            r#"{{
                "embedSeed": 7,
                "maxAttempts": 8,
                "embedTimeoutMs": null,
                "useExpTorsions": false,
                "useSmallRingTorsions": false,
                "useMacrocycleTorsions": false,
                "useMacrocycle14Bounds": false,
                "includeLegacyTorsionHeuristic": false,
                "stereoPolicy": "{stereo_policy}",
                "failOnUnevaluableStereo": false,
                "forceFieldPolicy": "{force_field}",
                "forceFieldMaxIterations": 200,
                "gateMmff94TorsionOop": false,
                "gateMmff94StretchBend": false,
                "ringTorsionPolicy": "{ring_torsion_policy}",
                "totalTimeoutMs": null
            }}"#
        )
    }

    #[test]
    fn success_path_has_expected_envelope_shape() {
        let mol = parse_smiles("CCCCCCCCCC").expect("decane"); // decane
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");

        assert_eq!(value["schemaVersion"], 1);
        assert_eq!(value["ok"], true);
        let result = &value["result"];
        for key in [
            "coords",
            "embedStats",
            "boundAdjustmentReport",
            "torsionKnowledgeReport",
            "ringTorsionEvidence",
            "torsionOptimizationReport",
            "stereoBefore",
            "stereoRepair",
            "stereoAfterRepair",
            "forceField",
            "finalStereo",
            "finalValidation",
            "elapsedMsByStage",
        ] {
            assert!(result.get(key).is_some(), "missing result.{key}");
        }
        assert_eq!(result["coords"].as_array().unwrap().len(), 10);
        assert_eq!(result["boundAdjustmentReport"], serde_json::Value::Null);
        assert_eq!(result["stereoRepair"], serde_json::Value::Null);
        assert_eq!(result["torsionOptimizationReport"], serde_json::Value::Null);
        assert_eq!(result["forceField"]["requestedForceField"], "none");
        assert_eq!(result["forceField"]["actualForceFieldUsed"], "none");
        assert_eq!(
            result["forceField"]["fallbackReason"],
            serde_json::Value::Null
        );
    }

    #[test]
    fn same_seed_is_reproducible() {
        // Wall-clock timing (`elapsedMsByStage`) is explicitly excluded from this
        // determinism requirement (it varies with system load, not with the
        // computation itself) -- strip it before comparing, matching the spec's
        // "wall-clock timing numeric equality is not required" instruction.
        fn without_timings(json: &str) -> serde_json::Value {
            let mut value: serde_json::Value = serde_json::from_str(json).unwrap();
            value["result"]["elapsedMsByStage"].take();
            value
        }
        let mol = parse_smiles("CCCCCCCCCC").expect("decane");
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json1 = embed_pipeline_v2_json(&mol, &config);
        let json2 = embed_pipeline_v2_json(&mol, &config);
        assert_eq!(without_timings(&json1), without_timings(&json2));
    }

    #[test]
    fn failure_path_has_expected_envelope_shape() {
        // A saturated small ring (cyclohexane) fused to an acyclic tail: requesting
        // small-ring torsions under the fail-closed policy is a reliably-reachable
        // typed failure (see the Python binding's own test suite for the same
        // fixture/rationale).
        let mol = parse_smiles("C1CCCCC1CCCCCCCCCCCC").expect("cyclohexane+chain");
        let config = r#"{
            "embedSeed": 7,
            "maxAttempts": 8,
            "embedTimeoutMs": null,
            "useExpTorsions": false,
            "useSmallRingTorsions": true,
            "useMacrocycleTorsions": false,
            "useMacrocycle14Bounds": false,
            "includeLegacyTorsionHeuristic": false,
            "stereoPolicy": "ignore",
            "failOnUnevaluableStereo": false,
            "forceFieldPolicy": "dreiding",
            "forceFieldMaxIterations": 200,
            "gateMmff94TorsionOop": false,
            "gateMmff94StretchBend": false,
            "ringTorsionPolicy": "fail_closed",
            "totalTimeoutMs": null
        }"#;
        let json = embed_pipeline_v2_json(&mol, config);
        let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");

        assert_eq!(value["schemaVersion"], 1);
        assert_eq!(value["ok"], false);
        let error = &value["error"];
        assert_eq!(
            error["cause"]["kind"],
            "ring_torsion_application_unsupported"
        );
        assert_eq!(error["stage"], "torsion_optimization");
        assert_eq!(error["coordsAreDiagnosticOnly"], true);
        assert!(error["lastKnownCoords"].as_array().is_some());
        assert!(!value["error"].as_object().unwrap().contains_key("coords"));
        assert!(error["diagnostics"]["embedStats"].is_object());
    }

    #[test]
    fn unknown_config_field_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config["notARealField"] = serde_json::json!(true);
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["stage"], "wasm_input_validation");
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn missing_required_field_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config.as_object_mut().unwrap().remove("maxAttempts");
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn pre_priority2_config_json_without_gate_stretch_bend_still_parses() {
        // Backward-compat regression (Priority 2, issue #227): a caller's
        // pre-existing 15-field config JSON (no `gateMmff94StretchBend` at
        // all) must keep working exactly as before, defaulting to `false`
        // -- not become a "missing required field" error just because a new
        // gate dimension was added. Contrast with
        // `missing_required_field_is_rejected` above: that one IS still a
        // real required field.
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config
            .as_object_mut()
            .unwrap()
            .remove("gateMmff94StretchBend");
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(
            value["ok"], true,
            "old caller config without gateMmff94StretchBend must still succeed: {value:?}"
        );
    }

    #[test]
    fn missing_nullable_timeout_field_is_still_rejected() {
        // embedTimeoutMs/totalTimeoutMs may be `null`, but the key itself must be
        // present -- a caller omitting it entirely must not be silently treated
        // as null (the double-Option serde trick this test pins).
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config.as_object_mut().unwrap().remove("embedTimeoutMs");
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn null_timeout_value_is_accepted() {
        let mol = parse_smiles("CC").expect("ethane");
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true);
    }

    #[test]
    fn present_timeout_value_is_accepted() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config["totalTimeoutMs"] = serde_json::json!(60_000);
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true);
    }

    #[test]
    fn unknown_enum_value_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config["stereoPolicy"] = serde_json::json!("not_a_real_policy");
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn out_of_range_integer_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let mut config: serde_json::Value =
            serde_json::from_str(&safe_config_json("none", "ignore", "fail_closed")).unwrap();
        config["maxAttempts"] = serde_json::json!(-1);
        let json = embed_pipeline_v2_json(&mol, &config.to_string());
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    #[test]
    fn malformed_json_is_rejected() {
        let mol = parse_smiles("CC").expect("ethane");
        let json = embed_pipeline_v2_json(&mol, "{not valid json");
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config");
    }

    /// Builds a `MolHandle` directly from `chematic_core`'s own builder, bypassing
    /// `parse_smiles`'s `enforce_wasm_molecule_size` check entirely -- that check
    /// itself constructs a `JsValue` on the error path, which aborts the process
    /// when run as a native (non-wasm32) test outside a real JS host. This lets
    /// `embed_pipeline_v2_json`'s *own* atom-limit check (which returns a plain
    /// `String`, never a `JsValue`) be exercised in isolation.
    fn oversized_mol_handle(atom_count: usize) -> MolHandle {
        let mut builder = chematic_core::MoleculeBuilder::new();
        for _ in 0..atom_count {
            builder.add_atom(chematic_core::Atom::new(chematic_core::Element::C));
        }
        MolHandle {
            inner: std::rc::Rc::new(builder.build()),
        }
    }

    #[test]
    fn atom_limit_exceeded_is_fail_closed() {
        let mol = oversized_mol_handle(WASM_MAX_ATOMS + 1);
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "atom_limit_exceeded");
        assert_eq!(value["error"]["cause"]["limit"], WASM_MAX_ATOMS);
        assert_eq!(value["error"]["stage"], "wasm_input_validation");
    }

    #[test]
    fn oversized_config_json_is_fail_closed() {
        let mol = parse_smiles("CC").expect("ethane");
        let padding = " ".repeat(WASM_MAX_INPUT_BYTES + 1);
        let json = embed_pipeline_v2_json(&mol, &padding);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["cause"]["kind"], "input_too_large");
    }

    #[test]
    fn finite_f64_serializes_as_number_non_finite_as_null() {
        assert_eq!(serde_json::to_string(&FiniteF64(1.5)).unwrap(), "1.5");
        assert_eq!(serde_json::to_string(&FiniteF64(f64::NAN)).unwrap(), "null");
        assert_eq!(
            serde_json::to_string(&FiniteF64(f64::INFINITY)).unwrap(),
            "null"
        );
        assert_eq!(
            serde_json::to_string(&FiniteF64(f64::NEG_INFINITY)).unwrap(),
            "null"
        );
        // Full round-trippable precision -- not rounded to a fixed number of
        // decimal places, unlike this crate's older 3D JSON helpers.
        let precise = 1.0 / 3.0;
        let rendered = serde_json::to_string(&FiniteF64(precise)).unwrap();
        let parsed: f64 = rendered.parse().unwrap();
        assert_eq!(parsed, precise);
    }

    fn coord_distance(coords: &serde_json::Value, i: usize, j: usize) -> f64 {
        let p = |idx: usize| {
            let c = &coords[idx];
            (
                c[0].as_f64().unwrap(),
                c[1].as_f64().unwrap(),
                c[2].as_f64().unwrap(),
            )
        };
        let (x1, y1, z1) = p(i);
        let (x2, y2, z2) = p(j);
        ((x1 - x2).powi(2) + (y1 - y2).powi(2) + (z1 - z2).powi(2)).sqrt()
    }

    /// Same property `test_conformer_ensemble.py`'s `_assert_consistent_indexing`
    /// checks on the Python side (issue #172): every real bond must have a
    /// plausible length in the returned coords, and no two atoms may coincide --
    /// exactly what breaks if a canonicalize-then-reparse bug (like
    /// `conformer_ensemble_json`'s still-open one) silently desyncs the atom index
    /// space between `mol` and the returned coordinates.
    fn assert_consistent_indexing(mol: &MolHandle, coords: &serde_json::Value) {
        let coords_arr = coords.as_array().expect("coords must be an array");
        assert_eq!(coords_arr.len(), mol.inner.atom_count());
        for (_bond_idx, bond) in mol.inner.bonds() {
            let d = coord_distance(coords, bond.atom1.0 as usize, bond.atom2.0 as usize);
            assert!(
                (0.6..=2.6).contains(&d),
                "bond {}-{} has length {d:.3} -- outside sane range [0.6, 2.6]; \
                 likely atom-index mismatch",
                bond.atom1.0,
                bond.atom2.0
            );
        }
        let n = coords_arr.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = coord_distance(coords, i, j);
                assert!(
                    d >= 0.4,
                    "atoms {i},{j} are {d:.3} apart -- degenerate/clashing geometry"
                );
            }
        }
    }

    #[test]
    fn atom_order_is_preserved_naphthalene() {
        // Naphthalene reorders under canonicalization (same fixture the Python
        // binding's atom-order test suite uses) -- coords must stay indexed to
        // THIS Mol's own atom order, not a reparsed/canonicalized copy.
        let mol = parse_smiles("c1ccc2ccccc2c1").expect("naphthalene");
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true);
        assert_consistent_indexing(&mol, &value["result"]["coords"]);
    }

    #[test]
    fn atom_order_is_preserved_branched_and_aspirin() {
        for smiles in ["CCC(C)C", "CC(=O)Oc1ccccc1C(=O)O"] {
            let mol = parse_smiles(smiles).unwrap_or_else(|_| panic!("{smiles}"));
            let config = safe_config_json("none", "ignore", "fail_closed");
            let json = embed_pipeline_v2_json(&mol, &config);
            let value: serde_json::Value = serde_json::from_str(&json).unwrap();
            assert_eq!(value["ok"], true, "{smiles}");
            assert_consistent_indexing(&mol, &value["result"]["coords"]);
        }
    }

    #[test]
    fn negative_control_reproduces_old_reparse_bug_shape() {
        // Hand-reconstructs the exact anti-pattern `conformer_ensemble_json` still
        // has (canonicalize -> reparse -> compute -> return coords indexed to the
        // REPARSED molecule) and confirms cross-indexing it against the ORIGINAL
        // mol is caught -- proving `assert_consistent_indexing` can actually detect
        // a real atom-order mismatch, not just that it happens to pass.
        let mol = parse_smiles("c1ccc2ccccc2c1").expect("naphthalene");
        let canonical = mol.canonical_smiles();
        assert_ne!(
            canonical, "c1ccc2ccccc2c1",
            "expected canonicalization to reorder naphthalene's atoms"
        );
        let reparsed = parse_smiles(&canonical).expect("reparse canonical form");
        let config = safe_config_json("none", "ignore", "fail_closed");
        let json = embed_pipeline_v2_json(&reparsed, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true);

        let result = std::panic::catch_unwind(|| {
            assert_consistent_indexing(&mol, &value["result"]["coords"]);
        });
        assert!(
            result.is_err(),
            "cross-indexing the reparsed molecule's coords against the original \
             mol should have failed a bond-length sanity check"
        );
    }

    // -----------------------------------------------------------------------
    // enforce_chirality (v0.14.0, issue #285's E/Z bound fix) -- WASM parity
    // -----------------------------------------------------------------------

    fn enforce_chirality_config_json(stereo_policy: &str, enforce_chirality: bool) -> String {
        format!(
            r#"{{
                "embedSeed": 0,
                "maxAttempts": 1,
                "embedTimeoutMs": null,
                "useExpTorsions": false,
                "useSmallRingTorsions": false,
                "useMacrocycleTorsions": false,
                "useMacrocycle14Bounds": false,
                "includeLegacyTorsionHeuristic": false,
                "stereoPolicy": "{stereo_policy}",
                "failOnUnevaluableStereo": false,
                "forceFieldPolicy": "none",
                "forceFieldMaxIterations": 200,
                "gateMmff94TorsionOop": false,
                "gateMmff94StretchBend": false,
                "ringTorsionPolicy": "fail_closed",
                "totalTimeoutMs": null,
                "enforceChirality": {enforce_chirality}
            }}"#
        )
    }

    #[test]
    fn enforce_chirality_true_fixes_but2ene_z_raw_embed() {
        // Direct WASM-level confirmation that `enforceChirality` in the JSON
        // config reaches distance_geometry_v2.rs's apply_declared_ez_bounds
        // (issue #285): but2ene_Z is the exact molecule that fix targets --
        // raw embedding (no force field) must satisfy declared E/Z once
        // enforceChirality is set, across multiple seeds, matching the
        // Rust-level corpus measurement and the Python binding's parity test.
        let mol = parse_smiles(r"C/C=C\C").expect("but2ene_Z");
        for seed in 0..5u64 {
            let config = format!(
                r#"{{
                    "embedSeed": {seed},
                    "maxAttempts": 1,
                    "embedTimeoutMs": null,
                    "useExpTorsions": false,
                    "useSmallRingTorsions": false,
                    "useMacrocycleTorsions": false,
                    "useMacrocycle14Bounds": false,
                    "includeLegacyTorsionHeuristic": false,
                    "stereoPolicy": "ignore",
                    "failOnUnevaluableStereo": false,
                    "forceFieldPolicy": "none",
                    "forceFieldMaxIterations": 200,
                    "gateMmff94TorsionOop": false,
                    "gateMmff94StretchBend": false,
                    "ringTorsionPolicy": "fail_closed",
                    "totalTimeoutMs": null,
                    "enforceChirality": true
                }}"#
            );
            let json = embed_pipeline_v2_json(&mol, &config);
            let value: serde_json::Value = serde_json::from_str(&json).unwrap();
            assert_eq!(value["ok"], true, "seed {seed}: {json}");
            assert_eq!(
                value["result"]["finalStereo"]["isFullySatisfied"], true,
                "seed {seed}: raw embed must already satisfy declared E/Z"
            );
        }
    }

    #[test]
    fn enforce_chirality_defaults_false_missing_field_still_parses() {
        // #[serde(default)] precedent (matches gateMmff94StretchBend): configs
        // written before this field existed must keep parsing, as `false`.
        let mol = parse_smiles("CCCCCCCCCC").expect("decane");
        let config = safe_config_json("none", "ignore", "fail_closed"); // no enforceChirality key
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true, "{json}");
    }

    #[test]
    fn enforce_chirality_with_repair_and_verify_is_allowed() {
        // Revised 2026-08-24 (issue #291 Step A): this combination was
        // previously rejected as invalid_configuration -- now validated (see
        // `chematic_3d::pipeline_v2`'s revised Stage 1 doc entry and
        // `crates/chematic-3d/examples/issue291_repair_policy_measurement.rs`).
        let mol = parse_smiles(r"C/C=C\C").expect("but2ene_Z");
        let config = enforce_chirality_config_json("repair_and_verify", true);
        let json = embed_pipeline_v2_json(&mol, &config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true, "{json}");
        assert_eq!(
            value["result"]["finalStereo"]["isFullySatisfied"], true,
            "{json}"
        );
    }

    // -----------------------------------------------------------------------
    // expand_implicit_h_through_pipeline / pipeline_v2_stereo_safe_config_json
    // (issue #291/#383) -- WASM parity
    // -----------------------------------------------------------------------

    const TESTOSTERONE: &str = "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O";

    #[test]
    fn expand_implicit_h_through_pipeline_requires_enforce_chirality() {
        let mol = parse_smiles(TESTOSTERONE).expect("testosterone");
        let config = r#"{
                "embedSeed": 0,
                "maxAttempts": 8,
                "embedTimeoutMs": null,
                "useExpTorsions": false,
                "useSmallRingTorsions": false,
                "useMacrocycleTorsions": false,
                "useMacrocycle14Bounds": false,
                "includeLegacyTorsionHeuristic": false,
                "stereoPolicy": "repair_and_verify",
                "failOnUnevaluableStereo": false,
                "forceFieldPolicy": "mmff94_with_uff_fallback",
                "forceFieldMaxIterations": 200,
                "gateMmff94TorsionOop": false,
                "gateMmff94StretchBend": false,
                "ringTorsionPolicy": "diagnostic_only",
                "totalTimeoutMs": null,
                "enforceChirality": false,
                "expandImplicitHThroughPipeline": true
            }"#;
        let json = embed_pipeline_v2_json(&mol, config);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false, "{json}");
        assert_eq!(
            value["error"]["cause"]["kind"], "invalid_configuration",
            "{json}"
        );
        assert_eq!(value["error"]["stage"], "validate_config", "{json}");
    }

    #[test]
    fn stereo_safe_config_json_has_expected_shape() {
        let json = pipeline_v2_stereo_safe_config_json("mmff94_with_uff_fallback", "fail_closed");
        let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
        assert_eq!(value["schemaVersion"], 1);
        assert_eq!(value["ok"], true);
        let config = &value["config"];
        assert_eq!(config["stereoPolicy"], "repair_and_verify");
        assert_eq!(config["enforceChirality"], true);
        assert_eq!(config["expandImplicitHThroughPipeline"], true);
        assert_eq!(config["forceFieldPolicy"], "mmff94_with_uff_fallback");
        assert_eq!(config["ringTorsionPolicy"], "fail_closed");
    }

    #[test]
    fn stereo_safe_config_json_rejects_unknown_force_field() {
        let json = pipeline_v2_stereo_safe_config_json("not_a_real_force_field", "fail_closed");
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], false, "{json}");
        assert_eq!(value["error"]["cause"]["kind"], "invalid_config", "{json}");
    }

    #[test]
    fn stereo_safe_config_json_fixes_testosterone_via_wasm_binding() {
        // Same seed/configuration already Rust-level tested and cross-checked
        // against an independent oracle in pipeline_v2.rs's own
        // stereo_safe_matches_the_hand_built_configuration_above test.
        let mol = parse_smiles(TESTOSTERONE).expect("testosterone");
        let generated =
            pipeline_v2_stereo_safe_config_json("mmff94_with_uff_fallback", "diagnostic_only");
        let mut config: serde_json::Value = serde_json::from_str(&generated).unwrap();
        config["config"]["embedSeed"] = serde_json::json!(0);
        let config_json = config["config"].to_string();

        let json = embed_pipeline_v2_json(&mol, &config_json);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["ok"], true, "{json}");
        assert_eq!(
            value["result"]["finalStereo"]["isFullySatisfied"], true,
            "{json}"
        );
        assert_eq!(value["result"]["coords"].as_array().unwrap().len(), 20);
    }

    #[test]
    fn expand_implicit_h_through_pipeline_is_noop_without_declared_stereo() {
        let mol = parse_smiles("CCCCCCCCCC").expect("decane"); // no declared stereo
        let base = embed_pipeline_v2_json(
            &mol,
            &enforce_chirality_config_json("repair_and_verify", true),
        );
        let expanded_config = {
            let mut c: serde_json::Value =
                serde_json::from_str(&enforce_chirality_config_json("repair_and_verify", true))
                    .unwrap();
            c["expandImplicitHThroughPipeline"] = serde_json::json!(true);
            c.to_string()
        };
        let expanded = embed_pipeline_v2_json(&mol, &expanded_config);
        let base_value: serde_json::Value = serde_json::from_str(&base).unwrap();
        let expanded_value: serde_json::Value = serde_json::from_str(&expanded).unwrap();
        assert_eq!(
            base_value["result"]["coords"],
            expanded_value["result"]["coords"]
        );
    }

    // -----------------------------------------------------------------------
    // Cross-binding parity: validation/pipeline_v2_wasm_parity_fixtures.json
    //
    // Generated once by scripts/gen_pipeline_v2_wasm_parity_fixtures.py via the
    // Python binding (crates/chematic-py/src/pipeline_v2.rs). Both the WASM
    // binding under test here and Python ultimately call the exact same
    // `chematic_3d::pipeline_v2::embed_pipeline_v2` -- so this is less "does the
    // algorithm agree" (trivially true) and more "does each binding's own JSON
    // conversion faithfully represent the same underlying evidence, for a fixed
    // molecule/config/seed." Node's pipeline_v2_parity.test.mjs reads the same
    // file and checks the WASM JS entry point the same way.
    // -----------------------------------------------------------------------

    fn parity_fixtures() -> serde_json::Value {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../validation/pipeline_v2_wasm_parity_fixtures.json");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
        serde_json::from_str(&text).expect("valid JSON")
    }

    fn config_json_from_fixture(fixture: &serde_json::Value) -> String {
        fixture["config"].to_string()
    }

    #[test]
    fn wasm_binding_matches_python_reference_fixtures() {
        let fixtures = parity_fixtures();
        for fixture in fixtures["fixtures"].as_array().unwrap() {
            let name = fixture["name"].as_str().unwrap();
            let smiles = fixture["smiles"].as_str().unwrap();
            let config_json = config_json_from_fixture(fixture);
            let mol = parse_smiles(smiles).unwrap_or_else(|_| panic!("{name}: parse {smiles}"));

            assert_eq!(
                mol.inner.atom_count() as i64,
                fixture["atomCount"].as_i64().unwrap(),
                "{name}: atom count"
            );

            let json = embed_pipeline_v2_json(&mol, &config_json);
            let value: serde_json::Value = serde_json::from_str(&json).unwrap();
            let expected_ok = fixture["ok"].as_bool().unwrap();
            assert_eq!(value["ok"], expected_ok, "{name}: ok");

            if expected_ok {
                let result = &value["result"];
                assert_eq!(
                    result["coords"].as_array().unwrap().len() as i64,
                    fixture["coordsLength"].as_i64().unwrap(),
                    "{name}: coords length"
                );
                assert_consistent_indexing(&mol, &result["coords"]);
                assert_eq!(
                    result["stereoBefore"]["nDeclared"], fixture["stereoBeforeDeclared"],
                    "{name}: stereoBefore.nDeclared"
                );
                assert_eq!(
                    result["finalStereo"]["nDeclared"], fixture["stereoAfterDeclared"],
                    "{name}: finalStereo.nDeclared"
                );
                assert_eq!(
                    result["finalStereo"]["nViolations"], fixture["stereoAfterViolations"],
                    "{name}: finalStereo.nViolations"
                );
                assert_eq!(
                    result["forceField"]["requestedForceField"], fixture["forceFieldRequested"],
                    "{name}: forceField.requestedForceField"
                );
                assert_eq!(
                    result["forceField"]["actualForceFieldUsed"], fixture["forceFieldActual"],
                    "{name}: forceField.actualForceFieldUsed"
                );
                assert_eq!(
                    !result["forceField"]["fallbackReason"].is_null(),
                    fixture["hasFallback"].as_bool().unwrap(),
                    "{name}: forceField fallback presence"
                );
                assert_eq!(
                    result["finalValidation"]["sound"], fixture["sound"],
                    "{name}: finalValidation.sound"
                );
            } else {
                assert_eq!(
                    value["error"]["stage"], fixture["stage"],
                    "{name}: error.stage"
                );
                assert_eq!(
                    value["error"]["cause"]["kind"], fixture["causeKind"],
                    "{name}: error.cause.kind"
                );
            }
        }
    }

    /// Same fixtures, but computed via the raw `chematic_3d::pipeline_v2` API
    /// directly (no JSON at all) -- confirms the WASM JSON conversion (tested
    /// above) is not itself hiding a divergence from what the underlying Rust
    /// call actually returns.
    #[test]
    fn raw_rust_api_matches_python_reference_fixtures() {
        let fixtures = parity_fixtures();
        for fixture in fixtures["fixtures"].as_array().unwrap() {
            let name = fixture["name"].as_str().unwrap();
            let smiles = fixture["smiles"].as_str().unwrap();
            let cfg = &fixture["config"];
            let mol = chematic_smiles::parse(smiles).unwrap_or_else(|_| panic!("{name}: parse"));

            let stereo_policy = match cfg["stereoPolicy"].as_str().unwrap() {
                "ignore" => pv2::StereoPolicy::Ignore,
                "verify_only" => pv2::StereoPolicy::VerifyOnly,
                "repair_and_verify" => pv2::StereoPolicy::RepairAndVerify,
                other => panic!("unknown stereoPolicy {other}"),
            };
            let ring_torsion_policy = match cfg["ringTorsionPolicy"].as_str().unwrap() {
                "fail_closed" => pv2::RingTorsionApplicationPolicy::FailClosed,
                "diagnostic_only" => pv2::RingTorsionApplicationPolicy::DiagnosticOnly,
                other => panic!("unknown ringTorsionPolicy {other}"),
            };
            let force_field_policy = match cfg["forceFieldPolicy"].as_str().unwrap() {
                "mmff94_bond_angle_strict" => ForceFieldPolicy::Mmff94BondAngleStrict,
                "mmff94_with_uff_fallback" => ForceFieldPolicy::Mmff94WithUffFallback,
                "uff_only" => ForceFieldPolicy::UffOnly,
                "dreiding" => ForceFieldPolicy::Dreiding,
                "none" => ForceFieldPolicy::None,
                other => panic!("unknown forceFieldPolicy {other}"),
            };
            let config = pv2::PipelineV2Config {
                embed: EmbedParameters {
                    random_seed: cfg["embedSeed"].as_u64().unwrap(),
                    max_attempts: cfg["maxAttempts"].as_u64().unwrap() as usize,
                    timeout_ms: cfg["embedTimeoutMs"].as_u64(),
                    use_exp_torsions: cfg["useExpTorsions"].as_bool().unwrap(),
                    use_small_ring_torsions: cfg["useSmallRingTorsions"].as_bool().unwrap(),
                    use_macrocycle_torsions: cfg["useMacrocycleTorsions"].as_bool().unwrap(),
                    use_macrocycle_14_bounds: cfg["useMacrocycle14Bounds"].as_bool().unwrap(),
                    enforce_chirality: cfg["enforceChirality"].as_bool().unwrap_or(false),
                    ..EmbedParameters::default()
                },
                torsion_optimization: TorsionOptimizationConfig::default(),
                include_legacy_torsion_heuristic: cfg["includeLegacyTorsionHeuristic"]
                    .as_bool()
                    .unwrap(),
                stereo_policy,
                fail_on_unevaluable_stereo: cfg["failOnUnevaluableStereo"].as_bool().unwrap(),
                force_field_policy,
                force_field_max_iterations: cfg["forceFieldMaxIterations"].as_u64().unwrap()
                    as usize,
                gate_mmff94_torsion_oop: cfg["gateMmff94TorsionOop"].as_bool().unwrap(),
                gate_mmff94_stretch_bend: cfg["gateMmff94StretchBend"].as_bool().unwrap_or(false),
                ring_torsion_policy,
                total_timeout_ms: cfg["totalTimeoutMs"].as_u64(),
                expand_implicit_h_through_pipeline: false,
            };

            let expected_ok = fixture["ok"].as_bool().unwrap();
            match pv2::embed_pipeline_v2(&mol, &config) {
                Ok(result) => {
                    assert!(expected_ok, "{name}: expected failure, got success");
                    assert_eq!(
                        result.coords.atom_count() as i64,
                        fixture["coordsLength"].as_i64().unwrap(),
                        "{name}: coords length"
                    );
                    assert_eq!(
                        result.stereo_before.n_declared() as i64,
                        fixture["stereoBeforeDeclared"].as_i64().unwrap(),
                        "{name}: stereo_before.n_declared"
                    );
                    assert_eq!(
                        force_field_policy_str(result.force_field.actual_force_field_used),
                        fixture["forceFieldActual"].as_str().unwrap(),
                        "{name}: force_field.actual_force_field_used"
                    );
                    assert_eq!(
                        result.force_field.fallback_reason.is_some(),
                        fixture["hasFallback"].as_bool().unwrap(),
                        "{name}: force_field fallback presence"
                    );
                    assert_eq!(
                        result.final_validation.sound,
                        fixture["sound"].as_bool().unwrap(),
                        "{name}: final_validation.sound"
                    );
                }
                Err(failure) => {
                    assert!(!expected_ok, "{name}: expected success, got failure");
                    assert_eq!(
                        snake_case_debug(&failure.stage),
                        fixture["stage"].as_str().unwrap(),
                        "{name}: stage"
                    );
                    assert_eq!(
                        failure_cause_kind_str(&failure.cause),
                        fixture["causeKind"].as_str().unwrap(),
                        "{name}: cause kind"
                    );
                }
            }
        }
    }

    /// `FailureCauseJson`'s serde tag value for a cause, without needing to
    /// serialize the whole struct first -- used only by the raw-API parity test.
    fn failure_cause_kind_str(cause: &pv2::PipelineV2FailureCause) -> String {
        let json = serde_json::to_value(failure_cause_json(cause)).unwrap();
        json["kind"].as_str().unwrap().to_string()
    }
}
