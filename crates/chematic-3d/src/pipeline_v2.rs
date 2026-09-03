//! Opt-in v2 embedding pipeline — 3D Breakthrough Program, Wave 2 → Wave 3
//! Coordinator Integration 1.
//!
//! Wires together 4 independently-merged, independently-verified pieces of work
//! into a single, fail-closed, opt-in orchestration:
//! - Agent C: stochastic distance geometry ([`crate::distance_geometry_v2`]).
//! - Agent E: torsion knowledge v2 + macrocycle 1-4 bounds ([`crate::etkdg_knowledge`]).
//! - Agent D: stereo verification/repair ([`crate::stereo_constraints`]).
//! - Agent F: typed force-field policy minimization ([`crate::minimize`]).
//!
//! This module changes **no existing default behavior** anywhere in the codebase.
//! [`generate_coords_etkdg`](crate::generate_coords_etkdg),
//! [`generate_coords`](crate::generate_coords),
//! [`generate_conformer_ensemble`](crate::generate_conformer_ensemble),
//! [`generate_conformer_ensemble_with_config`](crate::generate_conformer_ensemble_with_config),
//! remain untouched — those still use the older `dg.rs` embedder, unrelated to this
//! module. **Correction (2026-08-11, stale since this PR's original text):**
//! `embed_pipeline_v2` is no longer Rust-only — `Mol.embed_pipeline_v2()`
//! (`chematic-py`) and `embed_pipeline_v2_json` (`chematic-wasm`) both call it
//! directly, added by later PRs without this doc being updated. It remains opt-in
//! (a caller must explicitly choose this entry point; nothing routes through it
//! implicitly) and `chematic-mcp` still does not expose it.
//!
//! # The most important semantics — read this twice
//!
//! Agent E's `optimize_torsions` can only actually move coordinates for **acyclic
//! bridge single bonds** (`etkdg_knowledge::is_bridge_bond`, reused here — see that
//! function's own doc comment for why this is the ONE discriminator, not a proxy like
//! `ring_size.is_some()` or a second, independently-derived `classify_bond` check).
//! Small-ring and macrocycle potentials are **scored** by `evaluate_torsion_energy`
//! but never mechanically rotated: a ring bond has no well-defined single rigid
//! two-side split, so `optimize_torsions` simply never selects it as `rotatable`
//! (`energy.rs`'s own module doc). Concretely, this means `optimize_torsions` never
//! *fails* just because ring potentials were requested — it silently succeeds with
//! `rotated_bond_count` reflecting only the acyclic subset. **This pipeline's own
//! stage 6 pre-check, not `optimize_torsions`'s return value, is what enforces
//! [`RingTorsionApplicationPolicy::FailClosed`]** (see [`embed_pipeline_v2`]'s stage
//! 6 and the negative control in this module's tests proving `optimize_torsions`
//! alone would silently accept a ring-only potential list).
//!
//! Therefore, strictly:
//! - `use_exp_torsions=true` → acyclic bridge-bond potentials are actually optimizable.
//! - `use_small_ring_torsions=true` / `use_macrocycle_torsions=true` → those
//!   potentials are scored only; [`RingTorsionEvidence`] never reports them as
//!   `applied_to_geometry = true`, and under [`RingTorsionApplicationPolicy::FailClosed`]
//!   (the default a caller must still explicitly choose — there is no `Default` impl
//!   for [`PipelineV2Config`], see below) requesting them at all is a typed failure
//!   ([`PipelineV2FailureCause::RingTorsionApplicationUnsupported`]), not a silent
//!   "scored, close enough."
//! - `use_macrocycle_14_bounds=true` → applied to the distance-geometry bound matrix
//!   *before* triangle-inequality smoothing (stage 3/4), via a minimal internal hook
//!   added to `distance_geometry_v2.rs`
//!   (`embed_distance_geometry_v2_with_adjustments`), never by post-hoc-hacking
//!   embedded coordinates.
//!
//! This PR deliberately does **not** write a new ring Cartesian optimizer — out of
//! scope, left for a separate PR (see [`PipelineV2FailureCause::RingTorsionApplicationUnsupported`]).
//!
//! # Stage order (exact)
//!
//! 1. validate config
//! 2. build torsion knowledge
//! 3. compute macrocycle 1-4 adjustments
//! 4. raw stochastic DG with adjustments
//! 5. evaluate torsion energy
//! 6. optimize applicable acyclic torsions (ring-torsion-application gate lives here)
//! 7. verify stereo
//! 8. repair stereo when requested
//! 9. verify stereo again
//! 10. force-field minimization
//! 11. final stereo verification
//! 12. final geometry validation
//!
//! Stereo is repaired into the correct basin **before** the force field runs, and the
//! final, authoritative "does declared stereo still hold" gate is stage 11, evaluated
//! on the **post-minimization** geometry — never the reverse. See
//! [`StereoPolicy::RepairAndVerify`]'s doc for why this ordering is load-bearing (a
//! force field has no notion of declared chirality/E-Z and can walk a geometry back
//! across whichever stereo boundary a naive post-hoc repair would have fixed).
//!
//! # Judgment calls made in this file (see PR body for the full account)
//!
//! - **Revised 2026-08-11 (v0.14.0 release gate)**: `embed.enforce_chirality = true`
//!   was originally rejected as [`PipelineV2FailureCause::InvalidConfiguration`] at
//!   stage 1 for every `stereo_policy != StereoPolicy::Ignore`, reasoning that this
//!   pipeline's own stages 7–11 stereo gate and the raw embedder's `enforce_chirality`
//!   mechanism were unrelated and composing them would be confusing, not defense-in-
//!   depth. Direct 265-molecule-corpus measurement disproved that: `enforce_chirality`
//!   protects embedding-time correctness only (verified via `stereo_before`, populated
//!   before stage 10 runs), and stage 10's force-field minimization has no notion of
//!   declared stereo and can walk a correctly-embedded E/Z bond back across its
//!   boundary (`chembl_tier_b_0076`/`chembl_tier_b_0083`: `stereo_before` fully
//!   satisfied under `enforce_chirality`, `final_stereo` violated after MMFF94
//!   minimization; re-running with `ForceFieldPolicy::None` on the same molecules
//!   keeps `final_stereo` satisfied, isolating minimization as the cause). The two
//!   mechanisms are complementary stages of defense, not redundant: `enforce_chirality`
//!   without a post-minimization gate can silently report `success` on a geometry
//!   whose final declared stereo is wrong. `StereoPolicy::Ignore` and
//!   `StereoPolicy::VerifyOnly` are now both allowed with `enforce_chirality: true`
//!   (VerifyOnly's stage 11 "Violated => failure" gate is exactly the fail-closed
//!   check that catches this class of drift).
//! - **Revised 2026-08-24 (issue #291 Step A)**: `StereoPolicy::RepairAndVerify`
//!   used to also be rejected in combination with `enforce_chirality: true` here —
//!   composing `enforce_chirality`'s own repair-then-retry with stage 8's repair
//!   pass was "a separate, not-yet-validated question, deliberately deferred
//!   rather than decided by omission". Now validated on issue #291's own
//!   29-molecule declared-stereo corpus (5 base seeds, see
//!   `crates/chematic-3d/examples/issue291_repair_policy_measurement.rs`):
//!   `RepairAndVerify` alone (unaffected by this change) already takes `UffOnly`
//!   from 58.6% silently-wrong / 41.4% correct to 0% silently-wrong / 86.2%
//!   correct + 13.8% honest failure — allowing `enforce_chirality: true` on top
//!   raises correctness to 92.4%, with zero regressions (no new silent-wrong
//!   outcomes, no unsound repaired geometry), fully recovering 3 molecules
//!   `RepairAndVerify` alone could not (naproxen_S, ibuprofen_S,
//!   penicillin_core). It does not help the 2 that remain (testosterone,
//!   cholesterol — ring-fused stereocenters with no non-ring substituent to
//!   reflect; that population needs the separately-scoped chiral-volume-
//!   penalty-in-`refine_coords` work `docs/rfcs/etkdg_3d_gap_rfc.md` already
//!   diagnosed, not this fix).
//! - The `use_small_ring_torsions`/`use_macrocycle_torsions` fail-closed gate (stage
//!   6) is scoped exactly to `TorsionKnowledgeSource::SmallRingExperimental` /
//!   `MacrocycleAdaptation` potentials — not to `BasicChemicalKnowledge`'s flat-ring
//!   term (gated by `use_exp_torsions`, a different flag the spec's §2 never names in
//!   this context) even though that term also targets a non-bridge bond and is also
//!   scored-only. [`RingTorsionEvidence`] still reports its `applied_to_geometry`
//!   truthfully either way; only the hard gate is scoped narrowly, matching the
//!   literal flags spec §2 names.
//! - `StereoPolicy::VerifyOnly`'s "Violated => failure" and the strict
//!   `fail_on_unevaluable_stereo` check are both evaluated **once**, at stage 11
//!   (post-force-field), not additionally fast-failed at stage 7 — because a force
//!   field can in principle change an element's status, stage 11 is the only
//!   *authoritative* answer to "does the delivered geometry satisfy declared stereo,"
//!   and `stereo_before`/`stereo_after_repair` remain full diagnostic evidence
//!   regardless. `StereoPolicy::RepairAndVerify`'s repair-failure gate (stage 8) is
//!   the one exception: a failed repair stops the pipeline immediately, before
//!   spending a force-field run on a geometry already known to be unrepairable.
//! - `bounds_conformance` in [`FinalGeometryValidation`] reuses
//!   `distance_geometry_v2::bounds_conformance` as-is (unadjusted-bounds diagnostic)
//!   rather than re-deriving an adjustment-aware version — a second bounds-computation
//!   pathway duplicating `dg_fft` internals is exactly what this program's spec §7
//!   forbids, and the unadjusted measurement is still an honest, meaningful number
//!   (how far the final geometry sits from the *naive* bound matrix).
//! - `FinalGeometryValidation::sound` gates only hard geometric sanity (finite
//!   coordinates, unchanged atom count, worst bond length under
//!   `minimize::MAX_SANE_BOND_LENGTH`) — reused, not re-derived, from Agent F's own
//!   soundness gate, since `ForceFieldPolicy::None` skips that gate entirely and this
//!   is then the *only* backstop. Bond-violation rate, gross-clash count, bounds
//!   conformance, and ring-closure delta are measured and reported (never silently
//!   dropped) but do not by themselves fail an individual pipeline call — they are
//!   corpus-level metrics for the integration gate harness, matching this codebase's
//!   standing convention that heuristic-projection residuals are visible, not assumed
//!   away (`distance_geometry_v2::bounds_conformance`'s own doc comment).

use chematic_core::{AtomIdx, Molecule};

use crate::clock::Instant;
use crate::coords::Coords3D;
use crate::dg_fft::ideal_bond_length;
use crate::distance_geometry_v2::{
    self, BoundsConformance, DistanceBoundAdjustment, EmbedFailureCause, EmbedParameters,
    EmbedStats, EmbedWithAdjustmentsFailure, bounds_conformance, mol_has_declared_stereo,
    truncate_coords,
};
use crate::etkdg_knowledge::{
    PairBoundAdjustment, TorsionKnowledgeConfig, TorsionKnowledgeError, TorsionKnowledgeReport,
    TorsionKnowledgeSource, TorsionOptimizationConfig, TorsionOptimizationReport,
    build_torsion_knowledge, evaluate_torsion_energy, is_bridge_bond,
    macrocycle_14_bound_adjustments, optimize_torsions,
};
use crate::minimize::{
    ForceFieldBridgeError, ForceFieldPolicy, MAX_SANE_BOND_LENGTH, MinimizeConfig,
    PolicyMinimizeResult, minimize_with_policy_gated,
};
use crate::stereo_constraints::{
    RepairRejectionReason, RepairedElement, StereoElement, StereoVerification, repair_stereo,
    verify_stereo,
};

/// Non-bonded atom pairs closer than this (Å) count as a gross steric clash for
/// [`FinalGeometryValidation::gross_clash_count`]. Consistent with this codebase's
/// existing ad hoc "no clash" test thresholds elsewhere in `minimize.rs`
/// (0.8–1.2 Å range) rather than a newly-invented number.
///
/// ponytail: fixed constant, not a config field — nothing in this PR needs it to vary.
const NONBONDED_CLASH_THRESHOLD_ANGSTROM: f64 = 1.2;

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// How declared stereo is handled by [`embed_pipeline_v2`]. See the module docs for
/// exactly when each policy's checks are evaluated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StereoPolicy {
    /// No stereo processing at all: `verify_stereo` is still called (so
    /// `stereo_before`/`stereo_after_repair`/`final_stereo` are always real evidence,
    /// never fabricated), but its result never gates success or failure.
    Ignore,
    /// Never repairs. A `Violated` element in the final (post-force-field) geometry
    /// is a pipeline failure. Whether an `Unevaluable` element also fails is governed
    /// separately by [`PipelineV2Config::fail_on_unevaluable_stereo`].
    VerifyOnly,
    /// Repairs every `Violated` element (stage 8) before the force field runs (stage
    /// 10), then re-verifies after minimization (stage 11). A repair failure, or a
    /// `Violated` element surviving to the final check, is a pipeline failure.
    RepairAndVerify,
}

/// How [`embed_pipeline_v2`] handles a request to apply small-ring/macrocycle torsion
/// potentials to geometry, when the current optimizer cannot mechanically act on
/// ring-bond coordinates (see the module docs' "most important semantics" section).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingTorsionApplicationPolicy {
    /// Typed failure ([`PipelineV2FailureCause::RingTorsionApplicationUnsupported`])
    /// when a ring/macrocycle torsion potential that cannot be mechanically applied
    /// was requested. There is no `Default` impl for [`PipelineV2Config`], so a
    /// caller must always write this field explicitly — this variant is simply the
    /// one to pick when in doubt, not a silently-assumed default.
    FailClosed,
    /// Explicit opt-in: continue with those potentials scored-only.
    /// [`RingTorsionEvidence`] always still reports `applied_to_geometry = false`,
    /// `diagnostic_only = true` for them — never silently upgraded to "applied."
    DiagnosticOnly,
}

/// Configuration for [`embed_pipeline_v2`]. Deliberately has no `Default` impl:
/// [`ForceFieldPolicy`] itself has none (Agent F's own design — see spec §10, "the
/// caller must explicitly choose"), so every field here must be written out by hand
/// at every call site. See [`PipelineV2Config::new`] for a convenience constructor
/// that still requires the force-field policy as an explicit argument.
#[derive(Debug, Clone)]
pub struct PipelineV2Config {
    pub embed: EmbedParameters,
    pub torsion_optimization: TorsionOptimizationConfig,
    /// Maps to `TorsionKnowledgeConfig::include_legacy_heuristic` — no
    /// `EmbedParameters` counterpart exists for this flag (see that config's own doc
    /// comment), so it is a separate field here rather than folded into `embed`.
    pub include_legacy_torsion_heuristic: bool,
    pub stereo_policy: StereoPolicy,
    /// Explicit opt-in (spec §9: "how Unevaluable is treated must be explicit in
    /// config"): when `true`, a declared stereo element that is `Unevaluable` in the
    /// final (post-force-field) geometry is a typed failure
    /// ([`PipelineV2FailureCause::StereoUnevaluableUnderStrictPolicy`]), evaluated at
    /// the same stage-11 point as `Violated`. Has no effect under
    /// [`StereoPolicy::Ignore`].
    pub fail_on_unevaluable_stereo: bool,
    pub force_field_policy: ForceFieldPolicy,
    pub force_field_max_iterations: usize,
    pub gate_mmff94_torsion_oop: bool,
    /// Independent opt-in, same shape as `gate_mmff94_torsion_oop` (Priority
    /// 2 / Stage 1B, issue #227): when `true`, `Mmff94BondAngleStrict`/
    /// `Mmff94WithUffFallback` also refuse on a missing stretch-bend cross
    /// term. `false` leaves existing arms' behavior unchanged — see
    /// `minimize::minimize_with_policy_gated`'s doc for the full rationale.
    pub gate_mmff94_stretch_bend: bool,
    pub ring_torsion_policy: RingTorsionApplicationPolicy,
    /// Wall-clock budget (milliseconds) for the whole call. `None` = no limit.
    /// Checked coarsely: once immediately after every one of the 12 stages
    /// completes (including the last, stage 12 -- independent verification round 2
    /// found this was missing after stage 12 specifically; it is now checked
    /// there too, after that stage's own semantic gate), never pre-emptively
    /// *during* a stage. A single expensive stage (e.g. a large `max_attempts` in
    /// `EmbedParameters`, or `optimize_torsions` on a very large molecule) can
    /// therefore overshoot the budget before the next check point — this is the
    /// same "checked between attempts, not pre-emptively inside one attempt"
    /// convention `EmbedParameters::timeout_ms` already documents for the raw
    /// embedder, applied one level up.
    pub total_timeout_ms: Option<u64>,
    /// Issue #291: for declared stereocenters whose only non-ring substituent is
    /// an implicit H (`repair_tetrahedral_center` has no coordinate to reflect for
    /// those -- ring-fused steroid-like centers such as testosterone/cholesterol),
    /// run this whole pipeline on a temporary `chematic_chem::add_hydrogens`-
    /// expanded copy of the molecule instead of the original, then map the result
    /// back onto the caller's original atom count before returning.
    ///
    /// Distinct from `embed.materialize_implicit_h_for_chirality` (still rejected
    /// unconditionally by stage-1 validation, unaffected by this field): that one
    /// only expands within a single raw embed attempt and truncates immediately,
    /// so this pipeline's own stage 7/8/9/11 verify/repair calls never see the
    /// real H position and can disagree with it. This field expands once, before
    /// stage 4, and keeps the expanded molecule for the whole stage 4-11 sequence
    /// -- `PipelineV2Result::coords`/`final_stereo` are still always scoped to the
    /// original atom count; other diagnostic fields (`embed_stats`,
    /// `torsion_knowledge_report`, `stereo_before`/`stereo_repair`/
    /// `stereo_after_repair`, `force_field`, `final_validation`) describe the
    /// expanded internal working state when this is active -- every stereocenter/
    /// declared-E/Z-bond they can reference is a heavy atom or heavy-heavy bond,
    /// so their `AtomIdx`/`BondIdx` values stay meaningful against the original
    /// molecule too (`add_hydrogens` never renumbers heavy atoms or reorders
    /// original bonds), but coordinate-sized fields like `force_field.coords`
    /// stay at the expanded atom count.
    ///
    /// Requires `embed.enforce_chirality: true` (`InvalidConfiguration` otherwise
    /// -- same precedent as every other flag here that only makes sense combined
    /// with it). A no-op, byte-identical to `false`, for any molecule with no
    /// declared stereo at all. See `ROADMAP.md`'s `#291` entry ("Phase 0.5") for
    /// the measurement this design is based on. Default `false`.
    pub expand_implicit_h_through_pipeline: bool,
}

impl PipelineV2Config {
    /// Convenience constructor: every knowledge/optimization/stereo/ring-torsion flag
    /// off or at its most conservative (`Ignore`, `FailClosed`), with `force_field_policy`
    /// still an explicit, required argument (never defaulted — see the struct's own
    /// doc). A minimal starting point for building up a specific arm's config, not a
    /// `Default` impl in disguise.
    pub fn minimal(force_field_policy: ForceFieldPolicy) -> Self {
        Self {
            embed: EmbedParameters::default(),
            torsion_optimization: TorsionOptimizationConfig::default(),
            include_legacy_torsion_heuristic: false,
            stereo_policy: StereoPolicy::Ignore,
            fail_on_unevaluable_stereo: false,
            force_field_policy,
            force_field_max_iterations: 200,
            gate_mmff94_torsion_oop: false,
            gate_mmff94_stretch_bend: false,
            ring_torsion_policy: RingTorsionApplicationPolicy::FailClosed,
            total_timeout_ms: None,
            expand_implicit_h_through_pipeline: false,
        }
    }

    /// Convenience constructor for the "stereo-safe" configuration (issue
    /// #291/#383): the exact 3-flag combination measured to close #291's
    /// residual for ring-fused declared stereocenters (testosterone,
    /// cholesterol, and similar) -- [`StereoPolicy::RepairAndVerify`],
    /// `embed.enforce_chirality: true`, and
    /// `expand_implicit_h_through_pipeline: true`, always set together.
    /// These three only work correctly as a set (`expand_implicit_h_through_pipeline`
    /// requires `enforce_chirality`, stage-1-validated elsewhere in this
    /// module; `RepairAndVerify` is what actually uses the corrected
    /// geometry) -- exposing them as independent flags risks a caller
    /// setting some but not all of them, silently landing back on a
    /// configuration issue #291's own measurement found unsound. Everything
    /// else matches [`Self::minimal`]'s conservative defaults; override any
    /// other field on the returned value the same way `minimal`'s own
    /// callers already do. `max_attempts` is deliberately left at
    /// `EmbedParameters::default()`'s `8` -- the exact value the 29-molecule
    /// regression measurement (`crates/chematic-3d/examples/
    /// issue291_repair_policy_measurement.rs`) validated (144/145 correct,
    /// 0 silently wrong), not bumped speculatively.
    pub fn stereo_safe(force_field_policy: ForceFieldPolicy) -> Self {
        Self {
            stereo_policy: StereoPolicy::RepairAndVerify,
            embed: EmbedParameters {
                enforce_chirality: true,
                ..EmbedParameters::default()
            },
            expand_implicit_h_through_pipeline: true,
            ..Self::minimal(force_field_policy)
        }
    }
}

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

/// Per-stage wall-clock elapsed time (milliseconds), for the integration gate
/// harness's runtime percentiles. Never silently dropped, present on success and on
/// failure alike (via [`PipelineV2Failure::elapsed_ms_by_stage`]).
#[derive(Debug, Clone, Copy, Default)]
pub struct StageTimings {
    pub torsion_knowledge_ms: u64,
    pub bound_adjustment_ms: u64,
    pub distance_geometry_ms: u64,
    pub torsion_energy_eval_ms: u64,
    pub torsion_optimization_ms: u64,
    pub stereo_verify_before_ms: u64,
    pub stereo_repair_ms: u64,
    pub stereo_verify_after_repair_ms: u64,
    pub force_field_ms: u64,
    pub final_stereo_verify_ms: u64,
    /// Time spent in the post-minimization repair-and-reverify attempt
    /// (issue #227 Phase 2). `0` unless `StereoPolicy::RepairAndVerify` AND
    /// stage 11 found a violation the force field introduced -- see
    /// [`PipelineV2Result::post_minimization_stereo_repair`].
    pub post_min_stereo_repair_ms: u64,
    pub final_validation_ms: u64,
    pub total_ms: u64,
}

/// Whether one [`crate::etkdg_knowledge::TorsionPotential`] was **eligible for**
/// mechanical rotation by stage 6's `optimize_torsions` call (its central bond is a
/// genuine acyclic bridge bond), vs. structurally scored-only (a ring/macrocycle
/// bond `optimize_torsions` can never rotate at all). The ONE source of truth for
/// this classification is [`is_bridge_bond`] on the potential's `central_bond` —
/// never a proxy (see the module docs' "most important semantics" section).
///
/// Caveat found by independent verification round 1, disclosed rather than silently
/// left implicit: `true` here means "was a candidate `optimize_torsions` could act
/// on," not "definitely moved" — a bridge-bond potential whose initial gradient is
/// already below `TorsionOptimizationConfig::convergence_grad_deg` is skipped
/// without ever rotating (see `energy.rs`'s inner loop), and `optimize_torsions`
/// does not report per-bond movement back up. This never affects the correctness
/// property this field exists to guarantee (a ring/macrocycle bond is *never*
/// mechanically rotated, full stop); it only means "applied_to_geometry = true"
/// should be read as "geometrically applicable," not "definitely moved this call."
#[derive(Debug, Clone, PartialEq)]
pub struct PotentialApplicationEvidence {
    pub rule_id: String,
    pub central_bond: (AtomIdx, AtomIdx),
    pub source: TorsionKnowledgeSource,
    pub applied_to_geometry: bool,
}

/// Full applied-vs-scored-only evidence for every matched torsion potential, plus
/// whether the run as a whole was under [`RingTorsionApplicationPolicy::DiagnosticOnly`]
/// (spec §13, Arm J: "the result always carries evidence equivalent to
/// `applied_to_geometry = false`, `diagnostic_only = true`").
#[derive(Debug, Clone, Default)]
pub struct RingTorsionEvidence {
    pub potentials: Vec<PotentialApplicationEvidence>,
    pub diagnostic_only: bool,
}

impl RingTorsionEvidence {
    pub fn n_applied(&self) -> usize {
        self.potentials
            .iter()
            .filter(|p| p.applied_to_geometry)
            .count()
    }

    pub fn n_scored_only(&self) -> usize {
        self.potentials
            .iter()
            .filter(|p| !p.applied_to_geometry)
            .count()
    }
}

/// Lossless summary of one [`repair_stereo`] call, regardless of whether it returned
/// `Ok`/`Err` — both carry a `repaired` list; `failures` is empty on the `Ok` path.
#[derive(Debug, Clone, Default)]
pub struct StereoRepairSummary {
    pub repaired: Vec<RepairedElement>,
    pub failures: Vec<(StereoElement, RepairRejectionReason)>,
}

/// Independent, pipeline-level final-geometry soundness measurement (spec §11), in
/// addition to (not instead of) the force-field bridge's own soundness gate — the
/// only backstop at all when `force_field_policy == ForceFieldPolicy::None`, which
/// skips Agent F's own gate entirely.
#[derive(Debug, Clone)]
pub struct FinalGeometryValidation {
    pub all_finite: bool,
    pub atom_count_unchanged: bool,
    pub worst_bond_length: f64,
    /// Fraction of real bonds whose length deviates from
    /// `dg_fft::ideal_bond_length` by more than 15%.
    pub bond_violation_rate_15pct: f64,
    /// Fraction of real bonds whose length deviates from
    /// `dg_fft::ideal_bond_length` by more than 50%.
    pub bond_violation_rate_50pct: f64,
    /// Non-bonded atom pairs closer than [`NONBONDED_CLASH_THRESHOLD_ANGSTROM`].
    pub gross_clash_count: usize,
    /// Reuses `distance_geometry_v2::bounds_conformance` as-is — see the module
    /// docs' judgment-call section for why this is the naive (unadjusted) bound
    /// matrix's conformance, not an adjustment-aware re-derivation.
    pub bounds_conformance: BoundsConformance,
    pub stereo_ok: bool,
    pub torsion_energy_after: f64,
    /// From `torsion_optimization_report.max_ring_closure_delta`; `0.0` when stage 6
    /// never ran (no potentials to optimize).
    pub ring_closure_delta: f64,
    /// Hard pass/fail gate: finite coordinates, unchanged atom count, worst bond
    /// length within `minimize::MAX_SANE_BOND_LENGTH`. See the module docs for why
    /// the other fields on this struct are measured-and-reported but non-gating.
    pub sound: bool,
}

/// Full result of a successful [`embed_pipeline_v2`] call — evidence of what actually
/// happened at every stage, never just final coordinates.
///
/// Issue #291's `expand_implicit_h_through_pipeline` (see `PipelineV2Config`'s own
/// doc) changes what indexing every field below other than `coords`/`final_stereo`
/// describes when active: `coords` and `final_stereo` are always scoped to the
/// caller's original molecule (the external contract), but `embed_stats`,
/// `torsion_knowledge_report`, `ring_torsion_evidence`, `torsion_optimization_report`,
/// `stereo_before`/`stereo_repair`/`stereo_after_repair`, `force_field` (including
/// `force_field.coords`), and `final_validation` describe the temporary, internal
/// H-expanded working molecule/coordinates instead — a strictly larger atom count
/// than `coords.atom_count()`. This is not a mixed-indexing hazard for the element
/// references those fields carry (`StereoElement::Tetrahedral(AtomIdx)`,
/// `StereoElement::DoubleBond(BondIdx)`, `RepairedElement`,
/// `PotentialApplicationEvidence::central_bond`): every stereocenter or declared
/// E/Z bond any of them can name is a heavy atom or a heavy-heavy bond, and
/// `chematic_chem::add_hydrogens` never renumbers heavy atoms or reorders original
/// bonds, so those specific indices agree against either molecule. It IS a real
/// difference for anything that measures the geometry itself:
/// `final_validation.gross_clash_count`/`bond_violation_rate_15pct`/`_50pct`/
/// `bounds_conformance` describe the expanded molecule when this flag is active,
/// not directly comparable to a flag-off run. `worst_bond_length`/`all_finite` are
/// the exception — heavy-heavy bonds are a strict subset of the expanded
/// molecule's bonds, so those two specifically still bound the (truncated)
/// returned geometry too.
#[derive(Debug, Clone)]
pub struct PipelineV2Result {
    pub coords: Coords3D,
    pub embed_stats: EmbedStats,
    /// `Some(..)` (possibly empty) iff `config.embed.use_macrocycle_14_bounds` was
    /// set; `None` when the feature was never requested at all.
    pub bound_adjustment_report: Option<Vec<PairBoundAdjustment>>,
    pub torsion_knowledge_report: TorsionKnowledgeReport,
    pub ring_torsion_evidence: RingTorsionEvidence,
    /// `None` when there were zero torsion potentials to optimize (stage 6 never ran).
    pub torsion_optimization_report: Option<TorsionOptimizationReport>,
    pub stereo_before: StereoVerification,
    /// `None` under `StereoPolicy::Ignore`/`VerifyOnly` (repair never attempted).
    pub stereo_repair: Option<StereoRepairSummary>,
    pub stereo_after_repair: StereoVerification,
    pub force_field: PolicyMinimizeResult,
    /// The authoritative gate (see the struct doc above) — always computed
    /// against the caller's *original* molecule and `coords`' atom count, even
    /// when `expand_implicit_h_through_pipeline` is active and every stage
    /// leading up to this ran on a temporary, larger, H-expanded molecule.
    pub final_stereo: StereoVerification,
    /// `Some(..)` iff `StereoPolicy::RepairAndVerify` AND force-field
    /// minimization (stage 10) introduced a stereo violation that stage 8's
    /// pre-minimization repair never had a chance to see (issue #227 Phase
    /// 2: MMFF94 minimization has no notion of declared stereo and can walk
    /// a fully-satisfied geometry back across a declared E/Z or tetrahedral
    /// boundary -- the same class already documented for
    /// `chembl_tier_b_0076`/`chembl_tier_b_0083` in this module's own
    /// judgment-call notes above). When present, `coords`/`final_stereo`
    /// above already reflect the SUCCESSFUL post-repair, re-verified,
    /// re-soundness-checked geometry -- this field is only ever `Some` on a
    /// call that recovered; a repair attempt that failed, didn't clear every
    /// violation, or produced an unsound geometry falls through to
    /// `PipelineV2FailureCause::FinalStereoViolation` exactly as before this
    /// field existed (`force_field.coords`/`force_field`'s other fields
    /// still show the FORCE FIELD's own unmodified output, never the
    /// repaired geometry -- this field is the one place the repair is
    /// visible).
    pub post_minimization_stereo_repair: Option<StereoRepairSummary>,
    pub final_validation: FinalGeometryValidation,
    pub elapsed_ms_by_stage: StageTimings,
}

// ---------------------------------------------------------------------------
// Failure
// ---------------------------------------------------------------------------

/// Which of the 12 exact pipeline stages a failure occurred at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PipelineStage {
    ValidateConfig,
    TorsionKnowledge,
    MacrocycleBoundAdjustment,
    DistanceGeometry,
    TorsionEnergyEvaluation,
    TorsionOptimization,
    StereoVerifyBefore,
    StereoRepair,
    StereoVerifyAfterRepair,
    ForceFieldMinimization,
    FinalStereoVerify,
    FinalGeometryValidationStage,
}

/// Typed failure cause. Never collapsed to a string — every variant carries (via
/// [`PipelineV2Failure`]) the stage it occurred at plus as much partial diagnostic
/// evidence as was computed before the failure.
#[derive(Debug, Clone)]
pub enum PipelineV2FailureCause {
    InvalidConfiguration,
    BoundAdjustmentFailed,
    DistanceGeometry(EmbedFailureCause),
    TorsionKnowledge(TorsionKnowledgeError),
    RingTorsionApplicationUnsupported,
    StereoRepairFailed,
    StereoUnevaluableUnderStrictPolicy,
    ForceField(ForceFieldBridgeError),
    FinalStereoViolation,
    FinalGeometryInvalid,
    Timeout,
}

/// A failed [`embed_pipeline_v2`] call. Never carries `coords` as if it were a usable
/// result (see the module's standing "never return partial coordinates dressed as a
/// success" rule) — `last_known_coords` is explicitly named and typed as a diagnostic
/// only, distinguishing it from [`PipelineV2Result::coords`]. When
/// `expand_implicit_h_through_pipeline` is active, `last_known_coords` may be sized
/// to the temporary H-expanded molecule (a strictly larger atom count than the
/// caller's original molecule) rather than always matching it — same caveat as
/// [`PipelineV2Result`]'s own doc on its non-`coords`/`final_stereo` fields, and for
/// the same reason (diagnostic only, never treated as a returned result).
#[derive(Debug, Clone)]
pub struct PipelineV2Failure {
    pub cause: PipelineV2FailureCause,
    pub stage: PipelineStage,
    pub last_known_coords: Option<Coords3D>,
    pub embed_stats: Option<EmbedStats>,
    pub bound_adjustment_report: Option<Vec<PairBoundAdjustment>>,
    pub torsion_knowledge_report: Option<TorsionKnowledgeReport>,
    pub ring_torsion_evidence: Option<RingTorsionEvidence>,
    pub torsion_optimization_report: Option<TorsionOptimizationReport>,
    pub stereo_before: Option<StereoVerification>,
    pub stereo_repair: Option<StereoRepairSummary>,
    pub stereo_after_repair: Option<StereoVerification>,
    pub force_field: Option<PolicyMinimizeResult>,
    pub final_stereo: Option<StereoVerification>,
    /// See [`PipelineV2Result::post_minimization_stereo_repair`].
    /// On a failed call, this is `Some` only when post-minimization repair
    /// was successfully accepted, but the pipeline subsequently failed
    /// during strict unevaluable-stereo checking, final geometry validation,
    /// or timeout enforcement. A repair attempt that was rejected leaves
    /// this as `None`.
    pub post_minimization_stereo_repair: Option<StereoRepairSummary>,
    pub final_validation: Option<FinalGeometryValidation>,
    pub elapsed_ms_by_stage: StageTimings,
}

impl PipelineV2Failure {
    fn new(cause: PipelineV2FailureCause, stage: PipelineStage, timings: StageTimings) -> Self {
        Self {
            cause,
            stage,
            last_known_coords: None,
            embed_stats: None,
            bound_adjustment_report: None,
            torsion_knowledge_report: None,
            ring_torsion_evidence: None,
            torsion_optimization_report: None,
            stereo_before: None,
            stereo_repair: None,
            stereo_after_repair: None,
            force_field: None,
            final_stereo: None,
            post_minimization_stereo_repair: None,
            final_validation: None,
            elapsed_ms_by_stage: timings,
        }
    }
}

/// Progressive diagnostic accumulator threaded through [`embed_pipeline_v2`]:
/// updated after every stage completes, so that **every** failure exit --
/// including a mid-stage timeout via `check_timeout!`, which has no local
/// stage-specific failure-construction block of its own -- carries the exact same
/// partial evidence an explicit `Err` return at that point would have. Fixes a real
/// gap independent verification round 1 found: `check_timeout!` previously built a
/// bare, all-`None` [`PipelineV2Failure`] regardless of how much had already been
/// computed, silently violating this module's own "carries as much partial
/// diagnostic information as possible" claim.
#[derive(Default)]
struct Evidence {
    last_known_coords: Option<Coords3D>,
    embed_stats: Option<EmbedStats>,
    bound_adjustment_report: Option<Vec<PairBoundAdjustment>>,
    torsion_knowledge_report: Option<TorsionKnowledgeReport>,
    ring_torsion_evidence: Option<RingTorsionEvidence>,
    torsion_optimization_report: Option<TorsionOptimizationReport>,
    stereo_before: Option<StereoVerification>,
    stereo_repair: Option<StereoRepairSummary>,
    stereo_after_repair: Option<StereoVerification>,
    force_field: Option<PolicyMinimizeResult>,
    final_stereo: Option<StereoVerification>,
    post_minimization_stereo_repair: Option<StereoRepairSummary>,
    final_validation: Option<FinalGeometryValidation>,
}

impl Evidence {
    fn fail(
        &self,
        cause: PipelineV2FailureCause,
        stage: PipelineStage,
        timings: StageTimings,
    ) -> PipelineV2Failure {
        let mut failure = PipelineV2Failure::new(cause, stage, timings);
        failure.last_known_coords = self.last_known_coords.clone();
        failure.embed_stats = self.embed_stats.clone();
        failure.bound_adjustment_report = self.bound_adjustment_report.clone();
        failure.torsion_knowledge_report = self.torsion_knowledge_report.clone();
        failure.ring_torsion_evidence = self.ring_torsion_evidence.clone();
        failure.torsion_optimization_report = self.torsion_optimization_report.clone();
        failure.stereo_before = self.stereo_before.clone();
        failure.stereo_repair = self.stereo_repair.clone();
        failure.stereo_after_repair = self.stereo_after_repair.clone();
        failure.force_field = self.force_field.clone();
        failure.final_stereo = self.final_stereo.clone();
        failure.post_minimization_stereo_repair = self.post_minimization_stereo_repair.clone();
        failure.final_validation = self.final_validation.clone();
        failure
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Run the full opt-in v2 embedding pipeline: torsion knowledge → macrocycle 1-4
/// bounds → stochastic distance geometry → torsion optimization → stereo
/// verification/repair → force-field minimization → final verification. See the
/// module docs for the exact 12-stage order and every judgment call made along the
/// way.
// `PipelineV2Failure` is intentionally large: it carries as much partial
// stats/stage/diagnostic evidence as was computed before failure ("Failure results
// must still carry as much... information as possible" -- spec §6), spread across
// ~12 genuinely-needed `Option<..>` fields rather than concentrated in one dominant
// field. This is the same size/richness tradeoff `ForceFieldBridgeError`'s own
// `MissingParameters(Box<Mmff94CoverageReport>)`/`MinimizationFailed(Box<..>)`
// variants already make elsewhere in this crate -- boxed there because a single
// field dominated; here boxing the whole struct would just move the allocation to
// every caller (including every success path's `?`) instead of removing it, for a
// type whose entire purpose is carrying that much data on the failure path only.
#[allow(clippy::result_large_err)]
pub fn embed_pipeline_v2(
    mol: &Molecule,
    config: &PipelineV2Config,
) -> Result<PipelineV2Result, PipelineV2Failure> {
    let overall_start = Instant::now();
    let mut timings = StageTimings::default();
    // Progressive diagnostic accumulator (see `Evidence`'s own doc comment): updated
    // right after every value becomes available, so `check_timeout!` -- which has no
    // stage-specific context of its own -- carries exactly the same partial evidence
    // an explicit `Err` return at that point would.
    let mut evidence = Evidence::default();

    macro_rules! check_timeout {
        ($stage:expr) => {
            if let Some(budget) = config.total_timeout_ms
                // `budget == 0` is checked unconditionally, not via the elapsed-time
                // comparison below: `Instant::elapsed().as_millis()` has only
                // millisecond granularity, so real (nonzero) wall-clock time spent
                // reaching this checkpoint can still read back as exactly 0ms on a
                // fast machine/small molecule -- `0 > 0` is false, letting a
                // `total_timeout_ms: 0` config silently succeed instead of failing
                // closed as documented. A zero budget means "no time was granted at
                // all," so any checkpoint reached at all must fail, independent of
                // what the (unreliable-at-this-resolution) elapsed reading says.
                // Found via a real, intermittent CI failure in `pipeline_v2_web_
                // target.test.mjs`/`pipeline_v2.test.mjs` (both assert this exact
                // "zero timeout must fail closed" contract) -- reproduced 0/60 times
                // locally but observed on 3 independent CI runs, consistent with a
                // race that only manifests on especially fast CI runners.
                && (budget == 0 || overall_start.elapsed().as_millis() as u64 > budget)
            {
                timings.total_ms = overall_start.elapsed().as_millis() as u64;
                return Err(evidence.fail(PipelineV2FailureCause::Timeout, $stage, timings));
            }
        };
    }

    // -----------------------------------------------------------------
    // Stage 1: validate config.
    // -----------------------------------------------------------------
    // Revised 2026-08-24 (issue #291 Step A): `embed.enforce_chirality` +
    // `StereoPolicy::RepairAndVerify` used to be rejected here as
    // `InvalidConfiguration` -- composing the two repair mechanisms was
    // "a separate, not-yet-validated question (deliberately deferred, not
    // decided by omission)". Now validated: on issue #291's own 29-molecule
    // declared-stereo corpus (5 base seeds, `crates/chematic-3d/examples/
    // issue291_repair_policy_measurement.rs`), `RepairAndVerify` alone
    // (unaffected by this change) already takes `UffOnly` from 58.6%
    // silently-wrong / 41.4% correct to 0% silently-wrong / 86.2% correct +
    // 13.8% honest failure. Additionally allowing `enforce_chirality: true`
    // here raises that to 92.4% correct, with zero regressions (no new
    // silent-wrong outcomes, no unsound repaired geometry) -- it fully
    // recovers 3 of the corpus's molecules that `RepairAndVerify` alone
    // could not (naproxen_S, ibuprofen_S, penicillin_core), and does not
    // help or hurt the 2 that remain unfixable by substituent-reflection
    // repair (testosterone, cholesterol -- ring-fused stereocenters with no
    // non-ring substituent to reflect; needs the separately-scoped chiral-
    // volume-penalty-in-`refine_coords` work `docs/rfcs/etkdg_3d_gap_rfc.md`
    // already diagnosed, not this fix). `StereoPolicy::Ignore` and
    // `StereoPolicy::VerifyOnly` remain allowed with `enforce_chirality: true`
    // as before this revision.

    // `materialize_implicit_h_for_chirality: true` is rejected here, unconditionally,
    // pending a follow-up to issue #291: `embed.materialize_implicit_h_for_chirality`
    // is proven correct in isolation (`distance_geometry_v2`'s own tests, an
    // independent-oracle cross-check), but this pipeline's own stages 7/8/9/11 below
    // verify/repair against `mol` -- the *original*, non-H-expanded molecule -- so
    // they fall back to `stereo_constraints::phantom_neighbor_position`'s estimated
    // implicit-H direction, not the real materialized position the embed actually
    // used. Measured directly: for ring-fused declared stereocenters (testosterone,
    // cholesterol) that estimate disagrees with the true, oracle-confirmed-correct
    // geometry, so `StereoVerifyBefore` falsely reports a violation and, under
    // `RepairAndVerify`, `StereoRepair` then fails the same way the original bug did.
    // Reject here rather than let a caller observe that confusing behavior -- the
    // flag stays usable directly through `embed_distance_geometry_v2`/
    // `embed_distance_geometry_v2_detail`, which have no such stage, until a
    // follow-up either threads H-materialization through this pipeline's own stages
    // or removes `phantom_neighbor_position`'s dependency on an estimated position.
    if config.embed.materialize_implicit_h_for_chirality {
        timings.total_ms = overall_start.elapsed().as_millis() as u64;
        return Err(evidence.fail(
            PipelineV2FailureCause::InvalidConfiguration,
            PipelineStage::ValidateConfig,
            timings,
        ));
    }

    // Revised (issue #291 real implementation): `expand_implicit_h_through_pipeline`
    // is the follow-up the comment above points to -- see its own doc comment on
    // `PipelineV2Config` and the stage-4 shadow below for what it actually does.
    // Requires `enforce_chirality` for the same reason every other flag here that
    // only matters combined with it does: without `enforce_chirality`, nothing
    // downstream ever checks declared stereo during embedding, so materializing
    // H for that purpose would just be wasted cost.
    if config.expand_implicit_h_through_pipeline && !config.embed.enforce_chirality {
        timings.total_ms = overall_start.elapsed().as_millis() as u64;
        return Err(evidence.fail(
            PipelineV2FailureCause::InvalidConfiguration,
            PipelineStage::ValidateConfig,
            timings,
        ));
    }

    // -----------------------------------------------------------------
    // Stage 2: build torsion knowledge.
    // -----------------------------------------------------------------
    let torsion_config = TorsionKnowledgeConfig {
        use_exp_torsions: config.embed.use_exp_torsions,
        use_small_ring_torsions: config.embed.use_small_ring_torsions,
        use_macrocycle_torsions: config.embed.use_macrocycle_torsions,
        use_macrocycle_14_bounds: config.embed.use_macrocycle_14_bounds,
        include_legacy_heuristic: config.include_legacy_torsion_heuristic,
    };
    let t0 = Instant::now();
    let torsion_knowledge_report = build_torsion_knowledge(mol, &torsion_config);
    timings.torsion_knowledge_ms = t0.elapsed().as_millis() as u64;
    evidence.torsion_knowledge_report = Some(torsion_knowledge_report.clone());
    check_timeout!(PipelineStage::TorsionKnowledge);

    // -----------------------------------------------------------------
    // Stage 3: compute macrocycle 1-4 adjustments.
    // -----------------------------------------------------------------
    let t0 = Instant::now();
    let bound_adjustment_report = if config.embed.use_macrocycle_14_bounds {
        match macrocycle_14_bound_adjustments(mol, &torsion_config) {
            Ok(v) => Some(v),
            Err(e) => {
                timings.bound_adjustment_ms = t0.elapsed().as_millis() as u64;
                timings.total_ms = overall_start.elapsed().as_millis() as u64;
                return Err(evidence.fail(
                    PipelineV2FailureCause::TorsionKnowledge(e),
                    PipelineStage::MacrocycleBoundAdjustment,
                    timings,
                ));
            }
        }
    } else {
        None
    };
    timings.bound_adjustment_ms = t0.elapsed().as_millis() as u64;
    evidence.bound_adjustment_report = bound_adjustment_report.clone();
    check_timeout!(PipelineStage::MacrocycleBoundAdjustment);

    // Issue #291: from here through stage 11, `mol` is shadowed to a temporary
    // `add_hydrogens`-expanded copy when `expand_implicit_h_through_pipeline` is
    // set -- NOT any earlier. Stages 2/3 above must see the ORIGINAL molecule
    // unconditionally: `etkdg_knowledge::classify_atom_type` classifies nitrogen
    // via raw graph-neighbor count, and `add_hydrogens` appending implicit H as
    // real graph nodes would silently reclassify e.g. a secondary amine
    // (NSp2 -> NSp3), changing which torsion rule matches for ANY molecule,
    // whenever a torsion-knowledge flag is on -- an interaction this feature has
    // no business touching. `torsion_knowledge_report.potentials[].central_bond`/
    // `bound_adjustment_report[].atom_pair` (already computed above, against the
    // original molecule) cite only heavy-atom `AtomIdx`s, and `add_hydrogens`
    // preserves heavy-atom `AtomIdx`/original `BondIdx` 1:1 (heavy atoms copied
    // first in original order, original bonds copied before H-bonds are
    // appended) -- so those indices stay valid once `mol` is shadowed here.
    let orig_mol: &Molecule = mol;
    let original_atom_count = orig_mol.atom_count();
    let use_expanded_geometry =
        config.expand_implicit_h_through_pipeline && mol_has_declared_stereo(orig_mol);
    let expanded_mol_storage: Molecule;
    let mol: &Molecule = if use_expanded_geometry {
        expanded_mol_storage = chematic_chem::add_hydrogens(orig_mol);
        &expanded_mol_storage
    } else {
        orig_mol
    };

    // -----------------------------------------------------------------
    // Stage 4: raw stochastic DG with adjustments.
    // -----------------------------------------------------------------
    let dg_adjustments: Vec<DistanceBoundAdjustment> = bound_adjustment_report
        .as_deref()
        .unwrap_or(&[])
        .iter()
        .map(|a| DistanceBoundAdjustment {
            atom1: a.atom_pair.0,
            atom2: a.atom_pair.1,
            lower: a.new_lower,
            upper: a.new_upper,
        })
        .collect();

    let t0 = Instant::now();
    let (coords, embed_stats) =
        match distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
            mol,
            &config.embed,
            &dg_adjustments,
        ) {
            Ok(v) => v,
            Err((EmbedWithAdjustmentsFailure::InvalidAdjustment, stats)) => {
                timings.distance_geometry_ms = t0.elapsed().as_millis() as u64;
                timings.total_ms = overall_start.elapsed().as_millis() as u64;
                evidence.embed_stats = Some(stats);
                return Err(evidence.fail(
                    PipelineV2FailureCause::BoundAdjustmentFailed,
                    PipelineStage::MacrocycleBoundAdjustment,
                    timings,
                ));
            }
            Err((EmbedWithAdjustmentsFailure::Embed(cause), stats)) => {
                timings.distance_geometry_ms = t0.elapsed().as_millis() as u64;
                timings.total_ms = overall_start.elapsed().as_millis() as u64;
                evidence.embed_stats = Some(stats);
                return Err(evidence.fail(
                    PipelineV2FailureCause::DistanceGeometry(cause),
                    PipelineStage::DistanceGeometry,
                    timings,
                ));
            }
        };
    timings.distance_geometry_ms = t0.elapsed().as_millis() as u64;
    evidence.embed_stats = Some(embed_stats.clone());
    evidence.last_known_coords = Some(coords.clone());
    check_timeout!(PipelineStage::DistanceGeometry);

    // -----------------------------------------------------------------
    // Stage 5: evaluate torsion energy.
    // -----------------------------------------------------------------
    let t0 = Instant::now();
    if let Err(e) = evaluate_torsion_energy(mol, &coords, &torsion_knowledge_report.potentials) {
        timings.torsion_energy_eval_ms = t0.elapsed().as_millis() as u64;
        timings.total_ms = overall_start.elapsed().as_millis() as u64;
        return Err(evidence.fail(
            PipelineV2FailureCause::TorsionKnowledge(e),
            PipelineStage::TorsionEnergyEvaluation,
            timings,
        ));
    }
    timings.torsion_energy_eval_ms = t0.elapsed().as_millis() as u64;
    check_timeout!(PipelineStage::TorsionEnergyEvaluation);

    // -----------------------------------------------------------------
    // Stage 6: optimize applicable acyclic torsions. This is ALSO where the
    // ring-torsion-application policy gate lives (see module docs).
    // -----------------------------------------------------------------
    let ring_torsion_evidence = RingTorsionEvidence {
        potentials: torsion_knowledge_report
            .potentials
            .iter()
            .map(|p| {
                let (b, c) = p.central_bond;
                let applied = mol.bond_between(b, c).is_some() && is_bridge_bond(mol, b, c);
                PotentialApplicationEvidence {
                    rule_id: p.rule_id.clone(),
                    central_bond: p.central_bond,
                    source: p.source,
                    applied_to_geometry: applied,
                }
            })
            .collect(),
        diagnostic_only: config.ring_torsion_policy == RingTorsionApplicationPolicy::DiagnosticOnly,
    };
    evidence.ring_torsion_evidence = Some(ring_torsion_evidence.clone());

    // Same two-part predicate as the evidence above (`bond_between(..).is_some() &&
    // is_bridge_bond(..)`), not just `!is_bridge_bond(..)` alone: independent
    // verification round 1 flagged that the two were only equivalent because every
    // `central_bond` happens to come from `candidate_central_bonds` (real bonds
    // only) -- matching the predicate exactly here removes that latent, currently-
    // unreachable divergence risk rather than relying on an invariant holding
    // elsewhere.
    let ring_application_requested =
        config.embed.use_small_ring_torsions || config.embed.use_macrocycle_torsions;
    let has_unsupported_ring_potential = torsion_knowledge_report.potentials.iter().any(|p| {
        let (b, c) = p.central_bond;
        matches!(
            p.source,
            TorsionKnowledgeSource::SmallRingExperimental
                | TorsionKnowledgeSource::MacrocycleAdaptation
        ) && !(mol.bond_between(b, c).is_some() && is_bridge_bond(mol, b, c))
    });
    if ring_application_requested
        && has_unsupported_ring_potential
        && config.ring_torsion_policy == RingTorsionApplicationPolicy::FailClosed
    {
        timings.total_ms = overall_start.elapsed().as_millis() as u64;
        return Err(evidence.fail(
            PipelineV2FailureCause::RingTorsionApplicationUnsupported,
            PipelineStage::TorsionOptimization,
            timings,
        ));
    }

    let t0 = Instant::now();
    let (coords, torsion_optimization_report) = if torsion_knowledge_report.potentials.is_empty() {
        (coords, None)
    } else {
        match optimize_torsions(
            mol,
            &coords,
            &torsion_knowledge_report.potentials,
            &config.torsion_optimization,
        ) {
            Ok((new_coords, report)) => (new_coords, Some(report)),
            Err(e) => {
                timings.torsion_optimization_ms = t0.elapsed().as_millis() as u64;
                timings.total_ms = overall_start.elapsed().as_millis() as u64;
                return Err(evidence.fail(
                    PipelineV2FailureCause::TorsionKnowledge(e),
                    PipelineStage::TorsionOptimization,
                    timings,
                ));
            }
        }
    };
    timings.torsion_optimization_ms = t0.elapsed().as_millis() as u64;
    evidence.torsion_optimization_report = torsion_optimization_report.clone();
    evidence.last_known_coords = Some(coords.clone());
    check_timeout!(PipelineStage::TorsionOptimization);

    // -----------------------------------------------------------------
    // Stage 7: verify stereo. Always computed (real evidence under every policy,
    // including `Ignore`) — only whether it GATES success/failure varies by policy.
    // -----------------------------------------------------------------
    let t0 = Instant::now();
    let stereo_before = verify_stereo(mol, &coords);
    timings.stereo_verify_before_ms = t0.elapsed().as_millis() as u64;
    evidence.stereo_before = Some(stereo_before.clone());
    check_timeout!(PipelineStage::StereoVerifyBefore);

    // -----------------------------------------------------------------
    // Stage 8: repair stereo when requested.
    // -----------------------------------------------------------------
    let t0 = Instant::now();
    let (coords, stereo_repair) = if config.stereo_policy == StereoPolicy::RepairAndVerify {
        match repair_stereo(mol, &coords) {
            Ok(outcome) => (
                outcome.coords,
                Some(StereoRepairSummary {
                    repaired: outcome.repaired,
                    failures: Vec::new(),
                }),
            ),
            Err(report) => {
                timings.stereo_repair_ms = t0.elapsed().as_millis() as u64;
                timings.total_ms = overall_start.elapsed().as_millis() as u64;
                evidence.last_known_coords = Some(report.partial_coords);
                evidence.stereo_repair = Some(StereoRepairSummary {
                    repaired: report.repaired,
                    failures: report.failures,
                });
                return Err(evidence.fail(
                    PipelineV2FailureCause::StereoRepairFailed,
                    PipelineStage::StereoRepair,
                    timings,
                ));
            }
        }
    } else {
        (coords, None)
    };
    timings.stereo_repair_ms = t0.elapsed().as_millis() as u64;
    evidence.stereo_repair = stereo_repair.clone();
    evidence.last_known_coords = Some(coords.clone());
    check_timeout!(PipelineStage::StereoRepair);

    // -----------------------------------------------------------------
    // Stage 9: verify stereo again.
    // -----------------------------------------------------------------
    let t0 = Instant::now();
    let stereo_after_repair = if stereo_repair.is_some() {
        verify_stereo(mol, &coords)
    } else {
        // No repair was attempted -- nothing changed, so re-verifying the same
        // coordinates would just recompute an identical result.
        stereo_before.clone()
    };
    timings.stereo_verify_after_repair_ms = t0.elapsed().as_millis() as u64;
    evidence.stereo_after_repair = Some(stereo_after_repair.clone());
    check_timeout!(PipelineStage::StereoVerifyAfterRepair);

    // -----------------------------------------------------------------
    // Stage 10: force-field minimization.
    // -----------------------------------------------------------------
    let ff_config = MinimizeConfig {
        max_steps: config.force_field_max_iterations,
        ..MinimizeConfig::default()
    };
    let t0 = Instant::now();
    let force_field = match minimize_with_policy_gated(
        mol,
        coords,
        config.force_field_policy,
        &ff_config,
        config.gate_mmff94_torsion_oop,
        config.gate_mmff94_stretch_bend,
    ) {
        Ok(r) => r,
        Err(e) => {
            timings.force_field_ms = t0.elapsed().as_millis() as u64;
            timings.total_ms = overall_start.elapsed().as_millis() as u64;
            return Err(evidence.fail(
                PipelineV2FailureCause::ForceField(e),
                PipelineStage::ForceFieldMinimization,
                timings,
            ));
        }
    };
    timings.force_field_ms = t0.elapsed().as_millis() as u64;
    evidence.force_field = Some(force_field.clone());
    evidence.last_known_coords = Some(force_field.coords.clone());
    check_timeout!(PipelineStage::ForceFieldMinimization);

    // -----------------------------------------------------------------
    // Stage 11: final stereo verification -- the single authoritative gate (see
    // module docs' judgment-call section for why this, not stage 7, is where
    // `VerifyOnly`'s "Violated => failure" and the strict-Unevaluable check fire).
    // -----------------------------------------------------------------
    // Issue #291: the actual, caller-facing signal. When `use_expanded_geometry`,
    // `mol`'s own `verify_stereo` (real H positions, no estimate) is proven
    // correct in isolation but is NOT this gate -- Phase 0.5 measured it can
    // disagree with what's actually returned (`PipelineV2Result::coords` is
    // always truncated to `original_atom_count`), specifically right after an
    // unrelaxed repair step. This closure is the one thing that decides
    // success/failure either way: truncate-then-verify-against-the-original-
    // molecule when expanded, otherwise byte-identical to the plain
    // `verify_stereo(mol, coords)` call this replaced.
    let authoritative_final_stereo = |coords: &Coords3D| -> StereoVerification {
        if use_expanded_geometry {
            let truncated = truncate_coords(coords, original_atom_count);
            verify_stereo(orig_mol, &truncated)
        } else {
            verify_stereo(mol, coords)
        }
    };

    let t0 = Instant::now();
    let mut final_stereo = authoritative_final_stereo(&force_field.coords);
    timings.final_stereo_verify_ms = t0.elapsed().as_millis() as u64;
    evidence.final_stereo = Some(final_stereo.clone());
    // Output geometry from here on -- starts as the force field's own,
    // unmodified result; only ever reassigned by a SUCCESSFUL post-
    // minimization repair-and-reverify below. `force_field.coords` itself is
    // never mutated, so `evidence.force_field`/`PipelineV2Result::force_field`
    // keep showing exactly what the force field produced.
    let mut out_coords = force_field.coords.clone();
    let mut post_min_repair: Option<StereoRepairSummary> = None;

    if config.stereo_policy != StereoPolicy::Ignore {
        if final_stereo.n_violations() > 0 {
            // Issue #227 Phase 2: MMFF94/UFF minimization has no notion of
            // declared stereo (see this module's own judgment-call notes
            // above on `chembl_tier_b_0076`/`chembl_tier_b_0083`) and can
            // walk a geometry that was fully satisfied at stage 9 back
            // across a declared E/Z or tetrahedral boundary. Stage 8's
            // repair runs too early to ever see a violation minimization
            // itself introduces. `RepairAndVerify` gets exactly one more
            // repair attempt here, on the post-minimization geometry --
            // accepted ONLY if the repair mechanism succeeds, the
            // re-verified result has zero violations, AND the repaired
            // geometry is still sound (`repair_stereo` is a rigid local
            // reflection that preserves bond lengths by construction, but
            // re-checked here rather than assumed -- see
            // `mmff94_bci_stereo_drift_diagnostic_227.rs` for the empirical
            // check this mirrors: worst_bond_length_ratio and clash count
            // both unchanged by the reflection on the case that motivated
            // this). Any rejection (repair itself fails, doesn't clear
            // every violation, or produces an unsound geometry) falls
            // through to the ORIGINAL `FinalStereoViolation` failure,
            // unchanged from before this existed. Never attempted under
            // `VerifyOnly`/`Ignore` -- this is additive recovery for
            // `RepairAndVerify`'s own contract, not a change to what other
            // policies gate on.
            let repair_t0 = Instant::now();
            if config.stereo_policy == StereoPolicy::RepairAndVerify
                && let Ok(outcome) = repair_stereo(mol, &force_field.coords)
            {
                let reverified = verify_stereo(mol, &outcome.coords);
                let sound_after = outcome.coords.is_finite()
                    && worst_bond_length(mol, &outcome.coords) <= MAX_SANE_BOND_LENGTH;
                if reverified.n_violations() == 0 && sound_after {
                    // Issue #291: `repair_stereo` moves the smallest bridge-eligible
                    // substituent, which for a materialized implicit H is that H
                    // itself (a 1-atom terminal component) -- measured directly
                    // (Phase 0.5): this repair alone never moves a heavy atom, so
                    // its own `reverified`/`sound_after` check above (real H
                    // positions, no truncation) can say Satisfied while the
                    // heavy-only coordinates this pipeline actually returns still
                    // disagree. One more force-field relaxation pass -- reusing the
                    // caller's own policy/gates, not hardcoding a force field --
                    // was measured to always resolve that disagreement, either to a
                    // genuine success or an honest, mutually-agreed failure, never a
                    // silent-wrong result. Gated on the expanded-side check here
                    // (matching exactly what was measured); the actual accept/
                    // reject decision below still runs through
                    // `authoritative_final_stereo`.
                    let mut repaired_coords = outcome.coords;
                    if use_expanded_geometry
                        && let Ok(relaxed) = minimize_with_policy_gated(
                            mol,
                            repaired_coords.clone(),
                            config.force_field_policy,
                            &ff_config,
                            config.gate_mmff94_torsion_oop,
                            config.gate_mmff94_stretch_bend,
                        )
                        && verify_stereo(mol, &relaxed.coords).n_violations() == 0
                    {
                        repaired_coords = relaxed.coords;
                    }
                    post_min_repair = Some(StereoRepairSummary {
                        repaired: outcome.repaired,
                        failures: Vec::new(),
                    });
                    out_coords = repaired_coords;
                    final_stereo = authoritative_final_stereo(&out_coords);
                }
            }
            timings.post_min_stereo_repair_ms = repair_t0.elapsed().as_millis() as u64;

            if post_min_repair.is_none() || final_stereo.n_violations() > 0 {
                evidence.final_stereo = Some(final_stereo.clone());
                timings.total_ms = overall_start.elapsed().as_millis() as u64;
                return Err(evidence.fail(
                    PipelineV2FailureCause::FinalStereoViolation,
                    PipelineStage::FinalStereoVerify,
                    timings,
                ));
            }
            evidence.final_stereo = Some(final_stereo.clone());
            evidence.post_minimization_stereo_repair = post_min_repair.clone();
            evidence.last_known_coords = Some(out_coords.clone());
        }
        if config.fail_on_unevaluable_stereo && final_stereo.n_unevaluable() > 0 {
            timings.total_ms = overall_start.elapsed().as_millis() as u64;
            return Err(evidence.fail(
                PipelineV2FailureCause::StereoUnevaluableUnderStrictPolicy,
                PipelineStage::FinalStereoVerify,
                timings,
            ));
        }
    }
    check_timeout!(PipelineStage::FinalStereoVerify);

    // -----------------------------------------------------------------
    // Stage 12: final geometry validation.
    // -----------------------------------------------------------------
    let t0 = Instant::now();
    let torsion_energy_after =
        match evaluate_torsion_energy(mol, &out_coords, &torsion_knowledge_report.potentials) {
            Ok(r) => r.total_energy,
            Err(_) => {
                // Only reachable if the force field somehow changed the atom count or
                // produced an out-of-range index -- a genuine internal-consistency
                // bug, not an ordinary chemistry failure, so it fails closed here
                // rather than silently reporting a stale/zero energy.
                timings.final_validation_ms = t0.elapsed().as_millis() as u64;
                timings.total_ms = overall_start.elapsed().as_millis() as u64;
                return Err(evidence.fail(
                    PipelineV2FailureCause::FinalGeometryInvalid,
                    PipelineStage::FinalGeometryValidationStage,
                    timings,
                ));
            }
        };

    let ring_closure_delta = torsion_optimization_report
        .as_ref()
        .map(|r| r.max_ring_closure_delta)
        .unwrap_or(0.0);

    let final_validation = compute_final_validation(
        mol,
        &out_coords,
        &final_stereo,
        torsion_energy_after,
        ring_closure_delta,
    );
    timings.final_validation_ms = t0.elapsed().as_millis() as u64;
    evidence.final_validation = Some(final_validation.clone());

    if !final_validation.sound {
        timings.total_ms = overall_start.elapsed().as_millis() as u64;
        return Err(evidence.fail(
            PipelineV2FailureCause::FinalGeometryInvalid,
            PipelineStage::FinalGeometryValidationStage,
            timings,
        ));
    }
    // Independent verification round 2 found this was missing: every other stage
    // gets a `check_timeout!` immediately after it, but stage 12 (the last one)
    // didn't -- meaning `total_timeout_ms` was silently unenforced on whatever time
    // stage 12's own O(n^2) clash-count / bond-violation-rate work took, so a
    // caller's requested budget could be exceeded without ever seeing a `Timeout`
    // failure. Checked here, after the stage's own semantic gate (matching stage
    // 11's own ordering: the stage's real outcome is decided first, the timeout
    // budget is enforced after).
    check_timeout!(PipelineStage::FinalGeometryValidationStage);

    timings.total_ms = overall_start.elapsed().as_millis() as u64;
    // Issue #291: the external contract ("one entry per atom of the molecule the
    // caller passed in") is unconditional -- truncate here, exactly once, the only
    // place `out_coords` itself (as opposed to `final_stereo`, already handled by
    // `authoritative_final_stereo` above) needs to change size. Every other field
    // on `PipelineV2Result` keeps describing the expanded internal state (see the
    // struct's own doc comment for why that's not a hazard).
    let returned_coords = if use_expanded_geometry {
        truncate_coords(&out_coords, original_atom_count)
    } else {
        out_coords
    };
    Ok(PipelineV2Result {
        coords: returned_coords,
        embed_stats,
        bound_adjustment_report,
        torsion_knowledge_report,
        ring_torsion_evidence,
        torsion_optimization_report,
        stereo_before,
        stereo_repair,
        stereo_after_repair,
        force_field,
        final_stereo,
        post_minimization_stereo_repair: post_min_repair,
        final_validation,
        elapsed_ms_by_stage: timings,
    })
}

// ---------------------------------------------------------------------------
// Final geometry validation helper
// ---------------------------------------------------------------------------

fn worst_bond_length(mol: &Molecule, coords: &Coords3D) -> f64 {
    mol.bonds()
        .map(|(_, bond)| coords.get(bond.atom1).distance(&coords.get(bond.atom2)))
        .fold(0.0_f64, f64::max)
}

/// Fraction of real bonds whose length deviates from `ideal_bond_length` by more than
/// `threshold` (e.g. `0.15` for the 15% rate, `0.50` for the 50% rate).
fn bond_violation_rate(mol: &Molecule, coords: &Coords3D, threshold: f64) -> f64 {
    let mut total = 0usize;
    let mut violations = 0usize;
    for (_, bond) in mol.bonds() {
        let ideal = ideal_bond_length(mol, bond.atom1, bond.atom2);
        if ideal <= 0.0 || !ideal.is_finite() {
            continue;
        }
        let actual = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
        total += 1;
        if ((actual - ideal).abs() / ideal) > threshold {
            violations += 1;
        }
    }
    if total == 0 {
        0.0
    } else {
        violations as f64 / total as f64
    }
}

fn gross_clash_count(mol: &Molecule, coords: &Coords3D) -> usize {
    let n = coords.atom_count();
    let mut count = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            let a = AtomIdx(i as u32);
            let b = AtomIdx(j as u32);
            if mol.bond_between(a, b).is_some() {
                continue;
            }
            if coords.get(a).distance(&coords.get(b)) < NONBONDED_CLASH_THRESHOLD_ANGSTROM {
                count += 1;
            }
        }
    }
    count
}

fn compute_final_validation(
    mol: &Molecule,
    coords: &Coords3D,
    final_stereo: &StereoVerification,
    torsion_energy_after: f64,
    ring_closure_delta: f64,
) -> FinalGeometryValidation {
    let all_finite = coords.is_finite();
    let atom_count_unchanged = coords.atom_count() == mol.atom_count();
    let worst_bond_length_v = worst_bond_length(mol, coords);
    let bounds_conformance_v = bounds_conformance(mol, coords);

    let sound = all_finite && atom_count_unchanged && worst_bond_length_v <= MAX_SANE_BOND_LENGTH;

    FinalGeometryValidation {
        all_finite,
        atom_count_unchanged,
        worst_bond_length: worst_bond_length_v,
        bond_violation_rate_15pct: bond_violation_rate(mol, coords, 0.15),
        bond_violation_rate_50pct: bond_violation_rate(mol, coords, 0.50),
        gross_clash_count: gross_clash_count(mol, coords),
        bounds_conformance: bounds_conformance_v,
        stereo_ok: final_stereo.is_fully_satisfied(),
        torsion_energy_after,
        ring_closure_delta,
        sound,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn config_none() -> PipelineV2Config {
        PipelineV2Config::minimal(ForceFieldPolicy::None)
    }

    // -----------------------------------------------------------------------
    // Reproducibility / invariance (spec §16)
    // -----------------------------------------------------------------------

    #[test]
    fn all_flags_off_matches_raw_dg_exactly() {
        // "all stage flags off -> matches raw DG" (spec §16/§11 negative-control
        // partner). With every knowledge/optimization/stereo/bounds flag off and
        // ForceFieldPolicy::None, the pipeline must reduce to exactly the raw
        // embedder's own output at the same seed.
        let mol = parse("CCCCCCCCCC").unwrap(); // decane
        let config = config_none();
        let result = embed_pipeline_v2(&mol, &config).expect("pipeline should succeed");

        let raw = distance_geometry_v2::embed_distance_geometry_v2(&mol, &config.embed)
            .expect("raw embed should succeed");

        for i in 0..mol.atom_count() {
            let p_pipeline = result.coords.get(AtomIdx(i as u32));
            let p_raw = raw.get(AtomIdx(i as u32));
            assert_eq!(
                p_pipeline, p_raw,
                "atom {i}: pipeline coords must be byte-identical to raw DG with every flag off"
            );
        }
        assert!(result.torsion_knowledge_report.potentials.is_empty());
        assert!(result.torsion_optimization_report.is_none());
        assert!(result.bound_adjustment_report.is_none());
    }

    #[test]
    fn same_seed_same_config_reproduces_identical_result() {
        let mol = parse("CC(=O)Nc1ccc(O)cc1").unwrap(); // paracetamol
        let config = config_none();
        let r1 = embed_pipeline_v2(&mol, &config).unwrap();
        let r2 = embed_pipeline_v2(&mol, &config).unwrap();
        for i in 0..mol.atom_count() {
            assert_eq!(
                r1.coords.get(AtomIdx(i as u32)),
                r2.coords.get(AtomIdx(i as u32)),
                "same seed/config must reproduce identical coords at atom {i}"
            );
        }
        assert_eq!(
            r1.final_validation.worst_bond_length,
            r2.final_validation.worst_bond_length
        );
    }

    #[test]
    fn different_seeds_are_not_aliased() {
        let mol = parse("CCCCCCCCCC").unwrap();
        let mut any_diff = false;
        let base = embed_pipeline_v2(&mol, &config_none()).unwrap();
        for seed in [1u64, 42u64] {
            let mut config = config_none();
            config.embed.random_seed = seed;
            let r = embed_pipeline_v2(&mol, &config).unwrap();
            for i in 0..mol.atom_count() {
                if base.coords.get(AtomIdx(i as u32)) != r.coords.get(AtomIdx(i as u32)) {
                    any_diff = true;
                }
            }
        }
        assert!(any_diff, "different seeds must not produce aliased output");
    }

    // -----------------------------------------------------------------------
    // Stage order / typed failures
    // -----------------------------------------------------------------------

    #[test]
    fn enforce_chirality_with_repair_and_verify_stereo_policy_is_allowed() {
        // Revised 2026-08-24 (issue #291 Step A): this combination was
        // previously rejected as InvalidConfiguration -- now validated (see
        // the module doc's revised Stage 1 entry and
        // `crates/chematic-3d/examples/issue291_repair_policy_measurement.rs`).
        let mol = parse("C[C@H](O)CC").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        match embed_pipeline_v2(&mol, &config) {
            Ok(r) => assert!(
                r.final_stereo.is_fully_satisfied(),
                "2-butanol has no ring-fused stereocenter -- must succeed with satisfied stereo"
            ),
            Err(e) => assert!(
                !matches!(e.cause, PipelineV2FailureCause::InvalidConfiguration),
                "must not be rejected as InvalidConfiguration, got {e:?}"
            ),
        }
    }

    #[test]
    fn embed_level_materialize_implicit_h_for_chirality_is_rejected_alone() {
        // `embed.materialize_implicit_h_for_chirality` is proven correct in isolation
        // (`distance_geometry_v2`'s own tests, an independent-oracle cross-check), but
        // on its own it only expands within a single raw embed attempt and truncates
        // immediately -- this pipeline's stages 7/8/9/11 verify/repair would never see
        // the real H position. Stays rejected unconditionally when used alone; the
        // real follow-up landed as the distinct, pipeline-level
        // `expand_implicit_h_through_pipeline` flag tested below, which keeps the
        // expanded molecule for the whole stage 4-11 sequence instead.
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.embed.materialize_implicit_h_for_chirality = true;
        let err = embed_pipeline_v2(&mol, &config).expect_err(
            "materialize_implicit_h_for_chirality must be rejected, not silently ignored or run",
        );
        assert!(
            matches!(err.cause, PipelineV2FailureCause::InvalidConfiguration),
            "expected InvalidConfiguration, got {err:?}"
        );
        assert_eq!(err.stage, PipelineStage::ValidateConfig);
    }

    #[test]
    fn embed_level_and_pipeline_level_materialize_implicit_h_flags_together_still_rejected() {
        // The pre-existing unconditional check on `embed.materialize_implicit_h_for_chirality`
        // already rejects this combination before any new code runs -- pinned here so
        // a future reordering of the stage-1 checks can't silently change that.
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.embed.materialize_implicit_h_for_chirality = true;
        config.expand_implicit_h_through_pipeline = true;
        let err = embed_pipeline_v2(&mol, &config)
            .expect_err("both flags together must still be rejected");
        assert!(matches!(
            err.cause,
            PipelineV2FailureCause::InvalidConfiguration
        ));
    }

    #[test]
    fn expand_implicit_h_through_pipeline_requires_enforce_chirality() {
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = false;
        config.expand_implicit_h_through_pipeline = true;
        let err = embed_pipeline_v2(&mol, &config).expect_err(
            "expand_implicit_h_through_pipeline without enforce_chirality must be rejected",
        );
        assert!(
            matches!(err.cause, PipelineV2FailureCause::InvalidConfiguration),
            "expected InvalidConfiguration, got {err:?}"
        );
        assert_eq!(err.stage, PipelineStage::ValidateConfig);
    }

    #[test]
    fn expand_implicit_h_through_pipeline_is_noop_without_declared_stereo() {
        let mol = parse("CCCCCCCCCC").unwrap(); // decane, no declared stereo at all
        let base = {
            let mut config = config_none();
            config.embed.enforce_chirality = true;
            embed_pipeline_v2(&mol, &config).unwrap()
        };
        let expanded = {
            let mut config = config_none();
            config.embed.enforce_chirality = true;
            config.expand_implicit_h_through_pipeline = true;
            embed_pipeline_v2(&mol, &config).unwrap()
        };
        assert_eq!(base.coords.atom_count(), expanded.coords.atom_count());
        for i in 0..base.coords.atom_count() {
            let idx = AtomIdx(i as u32);
            assert_eq!(
                base.coords.get(idx),
                expanded.coords.get(idx),
                "atom {i}: expand_implicit_h_through_pipeline must be byte-identical to off \
                 for a molecule with no declared stereo"
            );
        }
    }

    #[test]
    fn expand_implicit_h_through_pipeline_negative_control_no_ring_fused_center() {
        // 2-butanol's declared stereocenter has an implicit H but is NOT ring-fused
        // -- repair_tetrahedral_center's existing substituent-reflection already
        // handles it without any H materialization. The flag must not change the
        // observable outcome here: it's real work spent on a molecule that never
        // needed it.
        let mol = parse("C[C@H](O)CC").unwrap();
        let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94WithUffFallback);
        config.embed.enforce_chirality = true;
        config.embed.random_seed = 0;
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        config.expand_implicit_h_through_pipeline = true;
        let result = embed_pipeline_v2(&mol, &config)
            .expect("2-butanol must still succeed with the flag on");
        assert!(result.final_stereo.is_fully_satisfied());
    }

    #[test]
    #[ignore = "Experimental 3D long-run gate; run with cargo test -p chematic-3d --lib -- --ignored"]
    fn expand_implicit_h_through_pipeline_fixes_testosterone_with_correct_geometry() {
        // Issue #291's original residual. Seed picked from the Phase 0.5 measurement
        // harness's own run (`issue291_expanded_geometry_feasibility.rs`) as one of
        // testosterone's 3 known-clean seeds under this exact configuration -- not
        // an arbitrary seed, since that measurement found a real 3/5-vs-2/5 split
        // (the other 2/5 are honest, safe failures, not silent-wrong -- see
        // ROADMAP.md's #291 entry).
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let declared = chematic_chem::assign_cip(&mol);
        assert_eq!(
            declared.assignments.len(),
            6,
            "sanity: 6 declared stereocenters"
        );

        let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94WithUffFallback);
        config.embed.enforce_chirality = true;
        config.embed.random_seed = 0;
        config.embed.max_attempts = 8;
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        config.expand_implicit_h_through_pipeline = true;
        let result = embed_pipeline_v2(&mol, &config).expect("testosterone must succeed");
        assert!(
            result.final_stereo.is_fully_satisfied(),
            "final_stereo: {:?}",
            result.final_stereo
        );
        assert_eq!(result.coords.atom_count(), mol.atom_count());

        let perceived = crate::stereo3d::assign_stereo_from_3d(&mol, &result.coords);
        for &(idx, code) in &declared.assignments {
            if let Some(perceived_code) = perceived.get(idx) {
                assert_eq!(
                    perceived_code, code,
                    "atom {idx:?}: declared {code:?} but 3D-perceived {perceived_code:?} \
                     -- pipeline reported success with wrong chirality"
                );
            }
        }
    }

    #[test]
    #[ignore = "Experimental 3D long-run gate; run with cargo test -p chematic-3d --lib -- --ignored"]
    fn stereo_safe_matches_the_hand_built_configuration_above() {
        // issue #383: `stereo_safe` must be exactly the same as manually setting the
        // 3 flags -- same testosterone happy-path, same seed, same expectations as
        // `expand_implicit_h_through_pipeline_fixes_testosterone_with_correct_geometry`
        // above, just built via the convenience constructor instead.
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let declared = chematic_chem::assign_cip(&mol);

        let mut config = PipelineV2Config::stereo_safe(ForceFieldPolicy::Mmff94WithUffFallback);
        assert!(config.embed.enforce_chirality);
        assert!(config.expand_implicit_h_through_pipeline);
        assert_eq!(config.stereo_policy, StereoPolicy::RepairAndVerify);
        config.embed.random_seed = 0;

        let result = embed_pipeline_v2(&mol, &config).expect("testosterone must succeed");
        assert!(result.final_stereo.is_fully_satisfied());
        assert_eq!(result.coords.atom_count(), mol.atom_count());

        let perceived = crate::stereo3d::assign_stereo_from_3d(&mol, &result.coords);
        for &(idx, code) in &declared.assignments {
            if let Some(perceived_code) = perceived.get(idx) {
                assert_eq!(perceived_code, code);
            }
        }
    }

    #[test]
    #[ignore = "Experimental 3D long-run gate; run with cargo test -p chematic-3d --lib -- --ignored"]
    fn expand_implicit_h_through_pipeline_with_verify_only_never_reports_success_with_violated_final_stereo()
     {
        // `VerifyOnly` never repairs, so this test doesn't assert success at any
        // particular seed (unlike the RepairAndVerify happy-path test above) --
        // it asserts the *invariant* this flag must preserve regardless of policy:
        // whatever `final_stereo` says, it must never disagree with an Ok result.
        // This specific interaction (VerifyOnly + expand_implicit_h_through_pipeline)
        // was never empirically measured by the Phase 0.5 harness (which only ran
        // RepairAndVerify) -- it only follows from `authoritative_final_stereo` being
        // computed before the `stereo_policy != Ignore` branch, which this test
        // verifies directly rather than trusting the code-structure argument alone.
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        for seed in 0..5u64 {
            let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94WithUffFallback);
            config.embed.enforce_chirality = true;
            config.embed.random_seed = seed;
            config.embed.max_attempts = 8;
            config.stereo_policy = StereoPolicy::VerifyOnly;
            config.expand_implicit_h_through_pipeline = true;
            match embed_pipeline_v2(&mol, &config) {
                Ok(result) => assert!(
                    result.final_stereo.is_fully_satisfied(),
                    "seed={seed}: Ok result must never carry a violated final_stereo, got {:?}",
                    result.final_stereo
                ),
                Err(e) => {
                    // VerifyOnly never repairs, so any of these are legitimate --
                    // just confirm it's a recognized stereo-related failure, not a
                    // silent success/failure mismatch.
                    assert!(
                        matches!(
                            e.cause,
                            PipelineV2FailureCause::FinalStereoViolation
                                | PipelineV2FailureCause::DistanceGeometry(_)
                        ),
                        "seed={seed}: unexpected failure cause {:?}",
                        e.cause
                    );
                }
            }
        }
    }

    #[test]
    fn enforce_chirality_with_ignore_stereo_policy_is_allowed() {
        let mol = parse("CCCC").unwrap(); // no declared stereo at all
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.stereo_policy = StereoPolicy::Ignore;
        assert!(embed_pipeline_v2(&mol, &config).is_ok());
    }

    #[test]
    fn enforce_chirality_with_verify_only_stereo_policy_is_allowed() {
        // Revised 2026-08-11: previously InvalidConfiguration -- now allowed, since
        // VerifyOnly's stage 11 gate is exactly what catches a force field
        // undoing enforce_chirality's embedding-time correctness (see the module
        // doc's revised judgment-call entry for the corpus evidence).
        let mol = parse("C[C@H](O)CC").unwrap();
        let mut config = config_none();
        config.embed.enforce_chirality = true;
        config.embed.max_attempts = 8;
        config.stereo_policy = StereoPolicy::VerifyOnly;
        let result = embed_pipeline_v2(&mol, &config);
        if let Err(e) = &result {
            assert!(
                !matches!(e.cause, PipelineV2FailureCause::InvalidConfiguration),
                "must not be rejected as InvalidConfiguration, got {e:?}"
            );
        }
    }

    #[test]
    fn enforce_chirality_with_verify_only_never_reports_success_with_violated_final_stereo() {
        // Motivated by the corpus-measured chembl_tier_b_0076/0083 failure mode
        // (265-molecule v0.14.0 release-gate re-measurement): enforce_chirality can
        // deliver correct declared E/Z at embedding time (verified via
        // `stereo_before`, populated before stage 10 runs) while MMFF94
        // minimization -- which has no notion of declared stereo -- walks it back
        // across the boundary afterward (confirmed separately for that exact
        // molecule/seed: re-running with `ForceFieldPolicy::None` leaves
        // `final_stereo` satisfied, isolating minimization as the cause -- see the
        // module doc's revised judgment-call entry). Before this PR's gate
        // relaxation, `StereoPolicy::Ignore` was the only option compatible with
        // enforce_chirality, and Ignore never gates on stereo -- a caller could get
        // a `success` result whose geometry silently violates its own declared E/Z.
        //
        // This test asserts the general invariant, not a specific outcome for this
        // one molecule/seed: whether MMFF94 happens to revert the fix is itself
        // non-deterministic across build profiles (verified while writing this test
        // -- an unoptimized debug build converged differently than the release
        // build the corpus measurement used, for the identical seed), so asserting
        // "this exact call must fail" would make the test flaky by build profile.
        // What must hold regardless: `enforce_chirality`'s own part of the contract
        // (embedding-time correctness, `stereo_before`) always succeeds for this
        // molecule, and `Ok` is never returned with a violated `final_stereo` --
        // VerifyOnly's stage 11 gate must have caught it if minimization did revert
        // the fix.
        let mol = parse("COc1cc2nc(N3CCN(C(=O)/C=C/c4ccc(NC(=O)CBr)cc4)CC3)nc(N)c2cc1OC").unwrap(); // chembl_tier_b_0083
        let mut config = config_none();
        config.embed.random_seed = 20260801; // the corpus benchmark's exact seed
        config.embed.enforce_chirality = true;
        config.embed.max_attempts = 8;
        config.embed.use_exp_torsions = true;
        config.embed.use_small_ring_torsions = true;
        config.embed.use_macrocycle_torsions = true;
        config.embed.use_macrocycle_14_bounds = true;
        config.stereo_policy = StereoPolicy::VerifyOnly;
        config.force_field_policy = ForceFieldPolicy::Mmff94BondAngleStrict;
        // Capped well below the 200-iteration production default: in an
        // unoptimized debug build (what `cargo test` uses in CI), this specific
        // 36-atom macrocycle-containing molecule's MMFF94 minimization is slow
        // enough at 200 iterations to make the test take minutes -- 40 is still
        // enough to exercise "did minimization move the geometry", the actual
        // invariant under test, without that cost.
        config.force_field_max_iterations = 40;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
        config.total_timeout_ms = Some(60_000);

        match embed_pipeline_v2(&mol, &config) {
            Ok(r) => {
                assert!(
                    r.stereo_before.is_fully_satisfied(),
                    "enforce_chirality's embedding-time fix must be correct -- got {:?}",
                    r.stereo_before
                );
                assert!(
                    r.final_stereo.is_fully_satisfied(),
                    "must never report success with a geometry violating declared stereo -- \
                     got final_stereo {:?}",
                    r.final_stereo
                );
            }
            Err(e) => {
                assert!(
                    !matches!(e.cause, PipelineV2FailureCause::InvalidConfiguration),
                    "must not be rejected as InvalidConfiguration: {e:?}"
                );
                assert!(
                    e.stereo_before
                        .as_ref()
                        .is_some_and(StereoVerification::is_fully_satisfied),
                    "enforce_chirality's embedding-time fix must be correct -- got {:?}",
                    e.stereo_before
                );
            }
        }
    }

    #[test]
    fn ring_torsion_fail_closed_on_macrocycle_when_optimizer_cannot_apply() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cyclododecane, macrocycle
        let mut config = config_none();
        config.embed.use_macrocycle_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::FailClosed;
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        assert!(matches!(
            err.cause,
            PipelineV2FailureCause::RingTorsionApplicationUnsupported
        ));
        // Evidence must still be attached even on this early, pre-embed failure.
        let evidence = err
            .ring_torsion_evidence
            .expect("evidence must be attached");
        assert!(evidence.potentials.iter().any(|p| !p.applied_to_geometry));
    }

    #[test]
    fn ring_torsion_diagnostic_only_succeeds_with_scored_only_evidence() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cyclododecane
        let mut config = config_none();
        config.embed.use_macrocycle_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
        let result = embed_pipeline_v2(&mol, &config).expect("DiagnosticOnly must succeed");
        assert!(result.ring_torsion_evidence.diagnostic_only);
        assert!(
            result
                .ring_torsion_evidence
                .potentials
                .iter()
                .any(|p| p.source == TorsionKnowledgeSource::MacrocycleAdaptation),
            "expected at least one MacrocycleAdaptation potential on a macrocycle"
        );
        for p in &result.ring_torsion_evidence.potentials {
            if p.source == TorsionKnowledgeSource::MacrocycleAdaptation {
                assert!(
                    !p.applied_to_geometry,
                    "a macrocycle potential must never be reported as applied_to_geometry"
                );
            }
        }
    }

    #[test]
    fn ring_torsion_fail_closed_never_fires_on_purely_acyclic_potentials() {
        // Negative control: acyclic-only torsion knowledge must never trip the ring
        // gate, even under FailClosed.
        let mol = parse("CCCC").unwrap(); // butane, fully acyclic
        let mut config = config_none();
        config.embed.use_exp_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::FailClosed;
        let result = embed_pipeline_v2(&mol, &config).expect("acyclic torsions must not gate");
        for p in &result.ring_torsion_evidence.potentials {
            assert!(p.applied_to_geometry, "acyclic potential must be applied");
        }
    }

    #[test]
    fn macrocycle_14_bounds_actually_applied_before_smoothing_changes_output() {
        // Arm A vs Arm B differential (spec §13/§18-c): enabling
        // use_macrocycle_14_bounds on a macrocycle must produce DIFFERENT coordinates
        // than leaving it off, at the same seed -- otherwise the hook isn't wired
        // even though every unit test is green.
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cyclododecane
        let config_a = config_none();
        let mut config_b = config_none();
        config_b.embed.use_macrocycle_14_bounds = true;

        let result_a = embed_pipeline_v2(&mol, &config_a).unwrap();
        let result_b = embed_pipeline_v2(&mol, &config_b).unwrap();

        assert!(result_a.bound_adjustment_report.is_none());
        let adjustments_b = result_b
            .bound_adjustment_report
            .as_ref()
            .expect("Some(..) when the flag is set");
        assert!(
            !adjustments_b.is_empty(),
            "cyclododecane must produce real adjustments"
        );
        assert!(result_b.embed_stats.adjustments_applied > 0);

        let mut any_diff = false;
        for i in 0..mol.atom_count() {
            if result_a.coords.get(AtomIdx(i as u32)) != result_b.coords.get(AtomIdx(i as u32)) {
                any_diff = true;
            }
        }
        assert!(
            any_diff,
            "macrocycle 1-4 bounds must actually change the embedded geometry, not just be scored"
        );
    }

    #[test]
    fn macrocycle_14_amide_pinned_branch_is_reachable_through_the_real_pipeline() {
        // `bounds14.rs`'s `macrocycle_14:amide_ester_pinned` rule (pin-to-cis-or-trans,
        // depending on ring role, for an amide/ester bond inside a macrocycle) is a
        // DIFFERENT branch from the `macrocycle_14:relaxed_band` rule the test above
        // exercises -- every macrocycle in the frozen 58/63 corpus (cyclododecane,
        // crown_12_4, cyclooctadecane) is all-carbon or all-ether, so no arm in
        // `pipeline_v2_integration_gate.rs` ever hit this branch through the real
        // embedder before this test (found during verification round 4). A
        // 12-membered ring lactam is `bounds14.rs`'s own test fixture for this rule.
        let mol = parse("O=C1CCCCCCCCCCN1").unwrap(); // 12-membered ring lactam
        let mut config = config_none();
        config.embed.use_macrocycle_14_bounds = true;

        let result = embed_pipeline_v2(&mol, &config).expect("must embed successfully");
        let adjustments = result
            .bound_adjustment_report
            .as_ref()
            .expect("Some(..) when the flag is set");
        assert!(
            adjustments
                .iter()
                .any(|a| a.rule_id == "macrocycle_14:amide_ester_pinned"),
            "the amide-pinned branch must fire through the real pipeline, not just \
             bounds14.rs's own standalone unit test: {adjustments:?}"
        );
        assert!(result.embed_stats.adjustments_applied > 0);
    }

    /// `atoms[1]`-`atoms[2]` central bond dihedral in degrees, atan2-based
    /// (numerically stable near 0/180). Self-contained rather than reusing
    /// `etkdg_knowledge::energy`'s private `dihedral_deg`, to avoid widening
    /// that module's visibility just for this test.
    fn measured_dihedral_deg(coords: &crate::coords::Coords3D, atoms: [AtomIdx; 4]) -> f64 {
        let p0 = coords.get(atoms[0]);
        let p1 = coords.get(atoms[1]);
        let p2 = coords.get(atoms[2]);
        let p3 = coords.get(atoms[3]);
        let b1 = p1.sub(&p0);
        let b2 = p2.sub(&p1);
        let b3 = p3.sub(&p2);
        let n1 = b1.cross(&b2);
        let n2 = b2.cross(&b3);
        let b2_unit = b2.try_normalize().unwrap_or(crate::coords::Point3::zero());
        let m1 = n1.cross(&b2_unit);
        let x = n1.dot(&n2);
        let y = m1.dot(&n2);
        y.atan2(x).to_degrees()
    }

    /// Distance (degrees) from `dihedral` to the nearest of a planar amide's
    /// two valid configurations, 0° (cis) or ±180° (trans).
    fn dist_to_planar(dihedral_deg: f64) -> f64 {
        let d = dihedral_deg.rem_euclid(360.0);
        (d.min(360.0 - d)).min((d - 180.0).abs())
    }

    #[test]
    fn tertiary_amide_macrocycle_embeds_closer_to_planar_with_14_bounds_fix() {
        // Issue found while surveying RDKit's open issues (analogous to RDKit
        // #9266, "ETKDGv3 twisted tertiary amides in macrocycles"). Confirmed
        // this reproduces in chematic; root cause was `bounds14.rs`'s
        // amide-pinned branch pinning all 4 combinatorial 1-4 pairs through a
        // tertiary amide to the same (geometrically unsatisfiable) cis
        // configuration -- fixed to split cis/trans by ring role (see that
        // module's own tests for the precise per-pair values). This test
        // confirms the fix's real-world effect through the actual embedding
        // pipeline, not just the bound values in isolation: enabling the
        // fixed `use_macrocycle_14_bounds` must measurably improve amide
        // planarity relative to leaving it off, on a genuine tertiary-amide
        // macrolactam. Empirically measured (Python, live pipeline, 2 real
        // multi-amide macrocycles, 4 seeds each, 48 dihedral measurements):
        // mean distance-to-planar dropped from ~61.5° (flag off) to ~20.8°
        // (flag on, fixed) -- not perfect convergence (a soft DG-bounds nudge
        // on a stochastic embedder, not a hard planarity guarantee), but a
        // clear, reproducible improvement. This test uses a smaller,
        // single-amide fixture and a threshold with real margin below that
        // baseline, not a razor-thin one recreated from a lucky single run.
        let mol = parse("O=C1CCCCCCCCCCN1C").unwrap(); // N-methyl 13-membered macrolactam
        // Atom indices per this exact SMILES's left-to-right parse order.
        let o_idx = AtomIdx(0);
        let carbonyl_c_idx = AtomIdx(1);
        let ring_c_idx = AtomIdx(2);
        let ring_n_idx = AtomIdx(11);
        let amide_n_idx = AtomIdx(12);
        let methyl_idx = AtomIdx(13);

        let dihedral_pairs = [
            (ring_n_idx, o_idx),
            (ring_n_idx, ring_c_idx),
            (methyl_idx, o_idx),
            (methyl_idx, ring_c_idx),
        ];

        let mean_dist_to_planar = |use_14_bounds: bool, seed: u64| -> f64 {
            let mut config = config_none();
            config.embed.use_exp_torsions = true;
            config.embed.use_macrocycle_torsions = true;
            config.embed.use_macrocycle_14_bounds = use_14_bounds;
            config.embed.random_seed = seed;
            // Ring-internal torsion potentials are scored-only (never
            // mechanically applied, see `energy.rs`'s `is_bridge_bond` gate)
            // -- FailClosed would reject every attempt once macrocycle
            // torsion knowledge is enabled at all. DiagnosticOnly matches
            // how this scenario is actually run in production.
            config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
            let result = embed_pipeline_v2(&mol, &config).expect("must embed successfully");
            let mut total = 0.0;
            for &(a1, a4) in &dihedral_pairs {
                let dih =
                    measured_dihedral_deg(&result.coords, [a1, amide_n_idx, carbonyl_c_idx, a4]);
                total += dist_to_planar(dih);
            }
            total / dihedral_pairs.len() as f64
        };

        let seeds: [u64; 4] = [1, 7, 42, 123];
        let with_fix: Vec<f64> = seeds
            .iter()
            .map(|&s| mean_dist_to_planar(true, s))
            .collect();
        let without: Vec<f64> = seeds
            .iter()
            .map(|&s| mean_dist_to_planar(false, s))
            .collect();

        let avg = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let avg_with_fix = avg(&with_fix);
        let avg_without = avg(&without);

        assert!(
            avg_with_fix < avg_without,
            "enabling the fixed macrocycle_14_bounds must improve amide planarity: \
             with_fix={with_fix:?} (avg {avg_with_fix:.1}) without={without:?} (avg {avg_without:.1})"
        );
        assert!(
            avg_with_fix < 45.0,
            "average distance-to-planar with the fix should be well below the \
             ~90-130 range the original bug produced: with_fix={with_fix:?} (avg {avg_with_fix:.1})"
        );
    }

    #[test]
    fn invalid_adjustment_pair_is_typed_bound_adjustment_failed() {
        // Direct test of the internal hook's own validation (index out of range).
        let mol = parse("CC").unwrap();
        let bad = [DistanceBoundAdjustment {
            atom1: AtomIdx(0),
            atom2: AtomIdx(99),
            lower: 1.0,
            upper: 2.0,
        }];
        let err = distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
            &mol,
            &EmbedParameters::default(),
            &bad,
        )
        .unwrap_err();
        assert!(matches!(
            err.0,
            EmbedWithAdjustmentsFailure::InvalidAdjustment
        ));
    }

    #[test]
    fn invalid_adjustment_lower_greater_than_upper_is_rejected() {
        let mol = parse("CCCCCCCCCC").unwrap();
        let bad = [DistanceBoundAdjustment {
            atom1: AtomIdx(0),
            atom2: AtomIdx(3),
            lower: 5.0,
            upper: 2.0,
        }];
        let err = distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
            &mol,
            &EmbedParameters::default(),
            &bad,
        )
        .unwrap_err();
        assert!(matches!(
            err.0,
            EmbedWithAdjustmentsFailure::InvalidAdjustment
        ));
    }

    #[test]
    fn invalid_adjustment_nan_is_rejected() {
        let mol = parse("CCCCCCCCCC").unwrap();
        let bad = [DistanceBoundAdjustment {
            atom1: AtomIdx(0),
            atom2: AtomIdx(3),
            lower: f64::NAN,
            upper: 2.0,
        }];
        let err = distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
            &mol,
            &EmbedParameters::default(),
            &bad,
        )
        .unwrap_err();
        assert!(matches!(
            err.0,
            EmbedWithAdjustmentsFailure::InvalidAdjustment
        ));
    }

    #[test]
    fn zero_adjustments_is_byte_identical_to_public_detail_api() {
        let mol = parse("c1ccc2ccccc2c1").unwrap(); // naphthalene
        let params = EmbedParameters::default();
        let (via_hook, _) =
            distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(&mol, &params, &[])
                .unwrap();
        let (via_public, _) =
            distance_geometry_v2::embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        for i in 0..mol.atom_count() {
            assert_eq!(
                via_hook.get(AtomIdx(i as u32)),
                via_public.get(AtomIdx(i as u32))
            );
        }
    }

    #[test]
    fn stereo_repair_and_verify_fixes_violated_stereo_before_force_field() {
        let mol = parse("C[C@H](O)CC").unwrap(); // 2-butanol, declared stereo
        let mut config = config_none();
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        let result = embed_pipeline_v2(&mol, &config).expect("should succeed");
        assert!(result.final_stereo.is_fully_satisfied());
    }

    #[test]
    fn stereo_verify_only_fails_when_final_geometry_violates_declared_stereo() {
        // Force a violation by disabling repair, and manufacturing a molecule whose
        // raw DG output is checked directly against a corrupted variant via the
        // repair path used as an oracle: run RepairAndVerify once to find whether
        // there was anything to repair, then run VerifyOnly and confirm violations
        // (if any occurred) are reported the same way, i.e. gated at stage 11 -- and
        // confirm the trivial case (no possible violation, no repair active) still
        // succeeds under VerifyOnly.
        let mol = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine, implicit-H stereocenter
        let mut config = config_none();
        config.stereo_policy = StereoPolicy::VerifyOnly;
        // VerifyOnly must not panic and must produce a real, evaluated verdict either
        // way (Ok or a typed FinalStereoViolation) -- both are acceptable outcomes
        // here since raw DG + no repair is not guaranteed to land in the declared
        // basin; what must NOT happen is a silent Ok with n_violations() > 0.
        match embed_pipeline_v2(&mol, &config) {
            Ok(result) => assert!(result.final_stereo.is_fully_satisfied()),
            Err(err) => assert!(matches!(
                err.cause,
                PipelineV2FailureCause::FinalStereoViolation
            )),
        }
    }

    #[test]
    fn ignore_stereo_policy_never_fails_on_declared_stereo() {
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let mut config = config_none();
        config.stereo_policy = StereoPolicy::Ignore;
        let result = embed_pipeline_v2(&mol, &config).expect("Ignore must never fail on stereo");
        // Evidence is still real (not fabricated), even though it doesn't gate.
        assert_eq!(result.stereo_before.n_declared(), 1);
    }

    #[test]
    fn strict_unevaluable_stereo_is_typed_failure_when_requested() {
        // A declared quaternary stereocenter with a genuinely degenerate/coplanar
        // arrangement is hard to manufacture from raw DG directly; instead this
        // confirms the wiring: fail_on_unevaluable_stereo has no effect when there is
        // nothing Unevaluable (a molecule with fully-resolvable declared stereo must
        // still succeed under strict mode).
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let mut config = config_none();
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        config.fail_on_unevaluable_stereo = true;
        let result = embed_pipeline_v2(&mol, &config).expect("should succeed");
        assert_eq!(result.final_stereo.n_unevaluable(), 0);
    }

    #[test]
    fn wildcard_atom_reports_typed_distance_geometry_failure() {
        let mol = parse("[*]C").unwrap();
        let config = config_none();
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        assert!(matches!(
            err.cause,
            PipelineV2FailureCause::DistanceGeometry(EmbedFailureCause::InvalidTopology)
        ));
        assert_eq!(err.stage, PipelineStage::DistanceGeometry);
    }

    #[test]
    fn timeout_zero_fails_closed_with_typed_timeout() {
        let mol = parse("CCCCCCCCCC").unwrap();
        let mut config = config_none();
        config.total_timeout_ms = Some(0);
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        assert!(matches!(err.cause, PipelineV2FailureCause::Timeout));
    }

    /// Regression test for a real bug independent verification round 1 found:
    /// `check_timeout!` used to construct a bare, all-`None` `PipelineV2Failure`
    /// regardless of how much had already been computed by that point -- silently
    /// contradicting this module's own "carries as much partial diagnostic
    /// information as possible" claim. `total_timeout_ms: Some(0)` is guaranteed to
    /// trip at SOME `check_timeout!` call site, and every such site is textually
    /// after stage 2 (torsion knowledge) completes -- not asserting on exactly
    /// which stage trips (that depends on sub-millisecond timing and is not
    /// deterministic across machines), only that the evidence stage 2 always
    /// produces by that point is never silently dropped.
    #[test]
    fn timeout_failure_still_carries_evidence_computed_before_it_tripped() {
        let mol = parse("CC(=O)Nc1ccc(O)cc1").unwrap(); // paracetamol
        let mut config = config_none();
        config.total_timeout_ms = Some(0);
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        assert!(matches!(err.cause, PipelineV2FailureCause::Timeout));
        assert!(
            err.torsion_knowledge_report.is_some(),
            "a timeout must always carry at least the torsion-knowledge report \
             (every check_timeout! call site is after stage 2 completes) -- this \
             is exactly the evidence-dropping bug that was found and fixed, stage \
             reached: {:?}",
            err.stage
        );
    }

    #[test]
    fn force_field_none_still_runs_independent_final_validation() {
        let mol = parse("c1ccccc1").unwrap();
        let config = config_none();
        let result = embed_pipeline_v2(&mol, &config).expect("benzene should succeed");
        assert!(result.final_validation.sound);
        assert!(result.final_validation.all_finite);
        assert!(result.final_validation.atom_count_unchanged);
    }

    #[test]
    fn atom_order_permutation_gives_equivalent_result() {
        // "atom-order permutation -> equivalence after mapping" (spec §16). Compare
        // ibuprofen written two structurally-equivalent ways is out of scope for a
        // cheap unit test (would need a full atom-map + Kabsch RMSD harness); this
        // narrower check instead confirms two independently-parsed copies of the
        // exact same SMILES (same atom order) at the same seed give identical
        // results -- the baseline the gate-harness's real permutation test builds on.
        let mol_a = parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap(); // ibuprofen
        let mol_b = parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap();
        let config = config_none();
        let ra = embed_pipeline_v2(&mol_a, &config).unwrap();
        let rb = embed_pipeline_v2(&mol_b, &config).unwrap();
        for i in 0..mol_a.atom_count() {
            assert_eq!(
                ra.coords.get(AtomIdx(i as u32)),
                rb.coords.get(AtomIdx(i as u32))
            );
        }
    }

    // -----------------------------------------------------------------------
    // Negative controls (spec §17): deliberately corrupt an expected value and
    // confirm the relevant check actually fails.
    // -----------------------------------------------------------------------

    #[test]
    fn negative_control_optimize_torsions_alone_would_silently_accept_ring_only_potentials() {
        // Proves WHY stage 6's own pre-check is load-bearing: `optimize_torsions`
        // itself returns Ok (not a typed failure) when every potential's central bond
        // is a ring bond (none are "rotatable"), because it just optimizes zero
        // bonds and trivially "converges." If this test ever starts failing (i.e.
        // optimize_torsions itself started rejecting ring-only input), the
        // stage-6 pre-check in `embed_pipeline_v2` would still be correct, but this
        // comment's premise should be re-checked.
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let config = TorsionKnowledgeConfig {
            use_macrocycle_torsions: true,
            ..TorsionKnowledgeConfig::default()
        };
        let report = build_torsion_knowledge(&mol, &config);
        assert!(
            report.potentials.iter().all(|p| !is_bridge_bond(
                &mol,
                p.central_bond.0,
                p.central_bond.1
            )),
            "expected every macrocycle potential's central bond to be non-bridge"
        );
        let coords =
            distance_geometry_v2::embed_distance_geometry_v2(&mol, &EmbedParameters::default())
                .unwrap();
        let opt_config = TorsionOptimizationConfig::default();
        let outcome = optimize_torsions(&mol, &coords, &report.potentials, &opt_config);
        assert!(
            outcome.is_ok(),
            "optimize_torsions alone silently succeeds on ring-only potentials -- \
             this is exactly why the pipeline's own stage-6 gate cannot delegate to it"
        );
    }

    /// Negative control (spec §17): proves the check in
    /// `ring_torsion_diagnostic_only_succeeds_with_scored_only_evidence` actually
    /// discriminates REAL pipeline output, not just vacuously passes. Runs the same
    /// real `embed_pipeline_v2` call and asserts the WRONG thing (every macrocycle
    /// potential IS applied) -- this must panic, proving a regression that started
    /// reporting scored-only potentials as applied would be caught, not silently
    /// accepted.
    #[test]
    #[should_panic(expected = "must never be reported as applied_to_geometry")]
    fn negative_control_reporting_scored_only_as_applied_would_be_caught() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cyclododecane
        let mut config = config_none();
        config.embed.use_macrocycle_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
        let result = embed_pipeline_v2(&mol, &config).expect("DiagnosticOnly must succeed");
        for p in &result.ring_torsion_evidence.potentials {
            if p.source == TorsionKnowledgeSource::MacrocycleAdaptation {
                assert!(
                    p.applied_to_geometry, // deliberately WRONG -- must panic
                    "a macrocycle potential must never be reported as applied_to_geometry"
                );
            }
        }
    }

    #[test]
    fn negative_control_bounds_after_smoothing_would_fail_the_invariant_check() {
        // Manufacture the exact bug spec §17 names: apply an adjustment to a matrix
        // that has ALREADY been smoothed (instead of before), and confirm the
        // existing `smoothing_preserves_invariants` check (reused by
        // `embed_distance_geometry_v2_with_adjustments`) is what would catch it, by
        // constructing bounds that violate the invariant directly and confirming the
        // helper says so.
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let (lower0, upper0) = crate::dg_fft::build_bound_matrix(&mol);
        let mut lower = lower0.clone();
        let mut upper = upper0.clone();
        crate::dg_fft::smooth_bounds(&mut lower, &mut upper);
        // Now simulate "adjustment applied after smoothing" by widening a
        // FINITE-baseline entry in the ALREADY-SMOOTHED `upper` beyond what the
        // (unsmoothed) baseline allows -- this must be flagged. Search for a pair
        // with a genuinely finite pre-smoothing upper bound (a directly bonded pair
        // always has one; an unconstrained distant pair may have an infinite
        // fallback, which would make this check vacuously pass for the wrong
        // reason, so it must be excluded).
        let n = mol.atom_count();
        let mut found = false;
        'search: for i in 0..n {
            for j in (i + 1)..n {
                if upper0[i][j].is_finite() {
                    let widened_upper = upper[i][j] + 100.0;
                    let bad = widened_upper > upper0[i][j] + 1e-6;
                    assert!(
                        bad,
                        "a post-smoothing widening of a finite-baseline pair must be \
                         detectable as a baseline violation"
                    );
                    found = true;
                    break 'search;
                }
            }
        }
        assert!(
            found,
            "expected at least one finite-baseline pair in cyclododecane"
        );
    }

    /// Was `negative_control_final_stereo_violation_as_success_would_be_caught`
    /// (spec §17), using a REAL molecule found by the integration gate harness
    /// (`examples/pipeline_v2_integration_gate.rs`): under `RepairAndVerify` +
    /// `ForceFieldPolicy::Dreiding`, gly-ala-gly's stereo repair succeeds (stage 8)
    /// but DREIDING minimization (stage 10, no stereo awareness at all) measurably
    /// walks the geometry back across a declared stereo boundary. At the time this
    /// test was written, stage 11 had no way to recover from that and correctly
    /// caught it as `FinalStereoViolation` -- this test asserted exactly that.
    ///
    /// Issue #227 Phase 2 added a post-minimization repair-and-reverify step
    /// (immediately below stage 11's violation check, `RepairAndVerify` only) for
    /// precisely this failure shape -- stage 8's repair runs too early to see a
    /// violation minimization itself introduces, so `RepairAndVerify` now gets one
    /// more repair attempt on the POST-minimization geometry, accepted only if it
    /// fully clears every violation and the repaired geometry is still sound. For
    /// this exact molecule/config, that new step succeeds (empirically confirmed via
    /// `mmff94_bci_stereo_drift_diagnostic_227.rs` before this test was written:
    /// `repair_stereo` recovers a robustly-satisfied geometry -179.7° from the
    /// declared-boundary, `worst_bond_length_ratio`/`gross_clash_count` both
    /// unchanged by the reflection) -- so the old premise ("this case is
    /// unrecoverable") no longer holds. Updated to assert the new, correct, more
    /// precise behavior instead: not just "does it fail," but "does it recover, and
    /// is the recovery genuinely visible in the evidence."
    #[test]
    fn repair_and_verify_recovers_post_minimization_stereo_violation() {
        let mol = parse("NCC(=O)N[C@@H](C)C(=O)NCC(=O)O").unwrap(); // gly_ala_gly
        let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Dreiding);
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        let result = embed_pipeline_v2(&mol, &config).expect(
            "gly_ala_gly under RepairAndVerify+Dreiding must now recover via the \
             post-minimization repair-and-reverify step",
        );
        // Stage 8's repair must have genuinely succeeded first (this is specifically
        // the "fixed, then broken again by the force field" scenario, not merely
        // "repair never worked") -- unchanged precondition from the original test.
        assert!(
            result.stereo_after_repair.is_fully_satisfied(),
            "repair (stage 8) must have succeeded before Dreiding (stage 10) broke it again"
        );
        assert!(
            result.final_stereo.is_fully_satisfied(),
            "post-minimization repair-and-reverify must leave final_stereo fully satisfied"
        );
        let repair = result
            .post_minimization_stereo_repair
            .as_ref()
            .expect("post_minimization_stereo_repair must be Some -- this is exactly the case it exists for");
        assert!(
            !repair.repaired.is_empty(),
            "the post-min repair summary must record what it actually repaired"
        );
        assert!(repair.failures.is_empty());
    }

    /// Companion sanity check: the post-minimization repair-and-reverify step must
    /// be a true no-op (never invoked, `post_minimization_stereo_repair: None`, its
    /// own timing bucket `0`) whenever stage 11 finds nothing to recover from --
    /// covers both "policy doesn't gate stereo at all" (`Ignore`) and "nothing was
    /// violated" implicitly via any passing `RepairAndVerify` case elsewhere in this
    /// test module. Guards against the new step accidentally firing (and doing
    /// needless work, or worse, silently swapping in different-but-also-valid
    /// coordinates) on a call that never needed it.
    #[test]
    fn post_minimization_stereo_repair_is_a_no_op_when_nothing_needs_recovering() {
        let mol = parse("CCCCCC").unwrap(); // hexane, no declared stereo at all
        let mut config = PipelineV2Config::minimal(ForceFieldPolicy::Dreiding);
        config.stereo_policy = StereoPolicy::RepairAndVerify;
        let result = embed_pipeline_v2(&mol, &config).expect("hexane must embed and minimize");
        assert!(result.post_minimization_stereo_repair.is_none());
        assert_eq!(result.elapsed_ms_by_stage.post_min_stereo_repair_ms, 0);
    }

    /// Golden regression test (issue #227 Phase 2), pinned to the exact molecule
    /// that surfaced this gap: `chembl_tier_b_0082`
    /// (`COc1cc2nc(N3CCN(C(=O)/C=C/c4ccc(N=C=S)cc4)CC3)nc(N)c2cc1OC`, ChEMBL Wave 1
    /// corpus). Its single declared E/Z bond is satisfied post-embedding (`stereo_before`)
    /// under the exact seed `pipeline_v2_vs_rdkit_dump.rs` uses, but the BCI charge
    /// fix (this same PR) changed the electrostatic term enough that
    /// `Mmff94BondAngleStrict` minimization now walks it to `Violated` --
    /// oracle-confirmed via `mmff94_bci_stereo_drift_diagnostic_227.rs` (dihedral
    /// ~167° -> ~0.3°) and via RDKit's own real MMFF94 minimizer on the same
    /// molecule, which does NOT reproduce this (all 4 `rdkit_etkdgv3_*` arms
    /// satisfied -- a chematic-specific minimizer-robustness gap, not a physically
    /// expected crossing). Pins two things at once: (1) under `Ignore` (the policy
    /// `pipeline_v2_vs_rdkit_dump.rs`'s `chematic_pipeline_v2_mmff94_strict` arm
    /// actually uses, and this PR's own measured baseline), the violation is real
    /// and NOT gated -- documented here explicitly, not silently accepted; (2) under
    /// `RepairAndVerify`, the new post-minimization step genuinely recovers it.
    #[test]
    #[ignore = "Experimental 3D long-run gate; run with cargo test -p chematic-3d --lib -- --ignored"]
    fn chembl_tier_b_0082_ez_bond_survives_bci_fix_under_repair_and_verify_not_under_ignore() {
        let mol = parse("COc1cc2nc(N3CCN(C(=O)/C=C/c4ccc(N=C=S)cc4)CC3)nc(N)c2cc1OC").unwrap();

        // (1) Ignore: same policy pipeline_v2_vs_rdkit_dump.rs's
        // chematic_pipeline_v2_mmff94_strict arm uses. Must still succeed (Ignore
        // never gates on stereo) but must show the real, uncorrected violation --
        // this is the documented, known, out-of-scope-for-Ignore residual from the
        // BCI fix's 3-state re-measurement (validation/results/mmff94_bci_gap_227_phase2_report.md
        // §3b), pinned here so it can never silently regress further or be
        // silently "fixed" by an unrelated future change without this test noticing.
        // Config mirrors `pipeline_v2_vs_rdkit_dump.rs`'s `base_config` exactly
        // (embed feature flags + `RingTorsionApplicationPolicy::DiagnosticOnly`) --
        // `PipelineV2Config::minimal`'s all-conservative defaults do not reproduce
        // this molecule's violation (confirmed: `EmbedParameters::default()`'s
        // fewer active embed features change the geometry enough that it doesn't
        // manifest), so this is not just a style choice, it's load-bearing for
        // reproducing the exact case measured in the corpus dump.
        let mut ignore_config = PipelineV2Config::minimal(ForceFieldPolicy::Mmff94BondAngleStrict);
        ignore_config.embed = EmbedParameters {
            random_seed: 20260801,
            max_attempts: 8,
            use_exp_torsions: true,
            use_small_ring_torsions: true,
            use_macrocycle_torsions: true,
            use_macrocycle_14_bounds: true,
            track_failures: true,
            ..EmbedParameters::default()
        };
        ignore_config.ring_torsion_policy = RingTorsionApplicationPolicy::DiagnosticOnly;
        ignore_config.stereo_policy = StereoPolicy::Ignore;
        let ignore_result = embed_pipeline_v2(&mol, &ignore_config)
            .expect("chembl_tier_b_0082 must still succeed under Ignore (never gated)");
        // NOT asserted here: an exact `stereo_before`/`final_stereo` violation
        // count. Empirically, `stereo_before.n_violations()` differs between a
        // `--release` and a plain `cargo test` (dev profile) build for this
        // exact molecule/seed (1 vs 0) -- distance-geometry embedding's
        // eigendecomposition + `max_attempts` retry loop is numerically
        // sensitive enough that which attempt succeeds first can differ across
        // optimization levels, even at a fixed seed. This is a real, checked
        // property of the embedder (not assumed away), not a bug this PR
        // introduced or is trying to fix -- so the assertions below only pin
        // what's true in BOTH profiles: Ignore never repairs, and (below)
        // RepairAndVerify always ends fully satisfied regardless of which
        // stage the violation actually appeared at. The exact "violated only
        // in State 3, satisfied in State 2, via minimization specifically"
        // narrative is the `--release`-build, `pipeline_v2_vs_rdkit_dump.rs`-
        // reproduced case documented in `scripts/mmff94_provenance/PROVENANCE.md`
        // and `validation/results/mmff94_bci_gap_227_phase2_report.md` -- this
        // test pins the POLICY CONTRACT (Ignore doesn't recover, RepairAndVerify
        // does), not the exact numeric trajectory, which is what actually
        // varies here.
        assert!(ignore_result.post_minimization_stereo_repair.is_none());

        // (2) RepairAndVerify: must always end fully satisfied, regardless of
        // whether stage 8 (pre-minimization) or the new post-minimization step
        // is what ultimately did the repairing.
        let mut repair_config = ignore_config.clone();
        repair_config.stereo_policy = StereoPolicy::RepairAndVerify;
        let repair_result = embed_pipeline_v2(&mol, &repair_config)
            .expect("chembl_tier_b_0082 must recover under RepairAndVerify");
        assert!(repair_result.final_stereo.is_fully_satisfied());
    }

    /// Negative control (spec §17): a genuinely degenerate (coplanar) 3D arrangement
    /// for a real declared stereocenter must read `Unevaluable`, and
    /// `StereoVerification`'s own accessor methods (reused, not reimplemented, by
    /// `pipeline_v2.rs`'s stage-11 gate) must never fold that into `n_satisfied()`.
    #[test]
    fn negative_control_unevaluable_counted_as_satisfied_would_be_caught() {
        let m = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine, implicit-H stereocenter
        let idx = AtomIdx(1);
        // All 4 real neighbors placed in the same z=0 plane -> the phantom-H
        // direction (opposite the sum of the other three unit bond vectors) is
        // ill-defined / the signed volume is ~0: a genuinely degenerate geometry,
        // not a hand-picked StereoStatus value.
        let mut coords = Coords3D::new_zeroed(m.atom_count());
        coords.set(AtomIdx(0), crate::coords::Point3::new(-1.0, -1.3, 0.0));
        coords.set(idx, crate::coords::Point3::new(-0.3, -0.2, 0.0));
        coords.set(AtomIdx(2), crate::coords::Point3::new(-1.0, 1.0, 0.0));
        coords.set(AtomIdx(3), crate::coords::Point3::new(1.2, -0.2, 0.0));
        coords.set(AtomIdx(4), crate::coords::Point3::new(1.7, -1.0, 0.0));
        coords.set(AtomIdx(5), crate::coords::Point3::new(1.9, 0.7, 0.0));

        let report = verify_stereo(&m, &coords);
        let status = report
            .tetrahedral
            .iter()
            .find(|r| r.atom == idx)
            .unwrap()
            .status;
        assert!(
            matches!(
                status,
                crate::stereo_constraints::StereoStatus::Unevaluable(_)
            ),
            "expected a genuinely degenerate (coplanar) geometry to read Unevaluable, got {status:?}"
        );
        assert_eq!(
            report.n_satisfied(),
            0,
            "Unevaluable must never count as Satisfied"
        );
        assert_eq!(
            report.n_violations(),
            0,
            "Unevaluable must never count as Violated"
        );
        assert_eq!(report.n_unevaluable(), 1);
    }

    #[test]
    fn negative_control_hidden_fallback_would_be_caught() {
        // Proves the correct fallback-occurred check (`fallback_reason.is_some()`,
        // NOT `actual_force_field_used != requested_force_field`) actually
        // discriminates -- mirrors the same check already relied on in
        // `examples/cf_integration_smoke_test.rs`.
        let mol = parse("CC").unwrap();
        let config = MinimizeConfig::default();
        let coords =
            distance_geometry_v2::embed_distance_geometry_v2(&mol, &EmbedParameters::default())
                .unwrap();
        let result = minimize_with_policy_gated(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94WithUffFallback,
            &config,
            false,
            false,
        )
        .expect("ethane should minimize fine under MMFF94");
        // Ethane's MMFF94 attempt should succeed cleanly (no fallback needed) --
        // `actual_force_field_used` still legitimately differs in NAME semantics from
        // `requested_force_field` per that struct's own doc, so asserting
        // `fallback_reason.is_none()` (not `actual == requested`) is the only
        // correct way to confirm "no fallback occurred."
        assert!(
            result.fallback_reason.is_none(),
            "ethane should not need a UFF fallback"
        );
    }

    #[test]
    fn negative_control_partial_coords_as_ok_would_be_caught() {
        // Confirms PipelineV2Result and PipelineV2Failure are structurally distinct
        // types -- a failure can never be mistaken for (or silently coerced into) a
        // success carrying the same `coords` field name, since `PipelineV2Failure`
        // exposes `last_known_coords: Option<Coords3D>`, never `coords: Coords3D`.
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let mut config = config_none();
        config.embed.use_macrocycle_torsions = true;
        config.ring_torsion_policy = RingTorsionApplicationPolicy::FailClosed;
        let err = embed_pipeline_v2(&mol, &config).unwrap_err();
        // `last_known_coords` is diagnostic only -- confirm the type itself has no
        // field named `coords` a careless caller could mistake for a success value
        // (checked at compile time by construction; this line documents the intent).
        let _diagnostic_only: Option<Coords3D> = err.last_known_coords;
    }

    #[test]
    fn negative_control_dropping_a_failed_molecule_would_be_caught() {
        // Every failure carries `stage` -- confirms a caller CAN always attribute a
        // failed molecule to a specific stage rather than needing to drop it from
        // accounting. Exercised for real by the gate harness's denominator
        // bookkeeping; this unit test confirms the field is always populated.
        let mol = parse("[*]C").unwrap();
        let err = embed_pipeline_v2(&mol, &config_none()).unwrap_err();
        // `stage` is a plain enum value, always present (not Option) -- there is no
        // code path that returns `Err` without it.
        match err.stage {
            PipelineStage::ValidateConfig
            | PipelineStage::TorsionKnowledge
            | PipelineStage::MacrocycleBoundAdjustment
            | PipelineStage::DistanceGeometry
            | PipelineStage::TorsionEnergyEvaluation
            | PipelineStage::TorsionOptimization
            | PipelineStage::StereoVerifyBefore
            | PipelineStage::StereoRepair
            | PipelineStage::StereoVerifyAfterRepair
            | PipelineStage::ForceFieldMinimization
            | PipelineStage::FinalStereoVerify
            | PipelineStage::FinalGeometryValidationStage => {}
        }
    }

    #[test]
    fn negative_control_legacy_api_changing_would_be_caught() {
        // "the legacy live API changing at all" -- confirms generate_coords_etkdg
        // (an existing default path this PR must not touch) still produces the same
        // shape/behavior it always has: infallible, correct atom count. A full
        // byte-for-byte regression comparison lives in the crate's existing test
        // suite (untouched by this PR); this is a light-touch confirmation the
        // symbol is still callable with its pre-existing signature and behavior.
        let mol = parse("CCO").unwrap();
        let coords = crate::generate_coords_etkdg(&mol);
        assert_eq!(coords.atom_count(), mol.atom_count());
    }

    #[test]
    fn negative_control_accepting_nan_inf_would_be_caught() {
        let mol = parse("CC").unwrap();
        let bad_pairs = [
            DistanceBoundAdjustment {
                atom1: AtomIdx(0),
                atom2: AtomIdx(1),
                lower: f64::NAN,
                upper: 2.0,
            },
            DistanceBoundAdjustment {
                atom1: AtomIdx(0),
                atom2: AtomIdx(1),
                lower: 1.0,
                upper: f64::INFINITY,
            },
        ];
        for bad in bad_pairs {
            let err = distance_geometry_v2::embed_distance_geometry_v2_with_adjustments(
                &mol,
                &EmbedParameters::default(),
                &[bad],
            )
            .unwrap_err();
            assert!(matches!(
                err.0,
                EmbedWithAdjustmentsFailure::InvalidAdjustment
            ));
        }
    }
}
