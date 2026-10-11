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
    /// with it). Expansion is activated only when a declared tetrahedral center
    /// has an implicit H; E/Z-only molecules and tetrahedral centers without an
    /// implicit H stay on the original heavy-atom path. See
    /// `ROADMAP.md`'s `#291` entry ("Phase 0.5") for the measurement this design is
    /// based on. Default `false`.
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
            // Batch 14: the A6 molecules stop at the iteration limit, not at a
            // stationary-point problem (265 rows: 300 -> 199 converged, 1000 ->
            // 263, 23% more time); RDKit's own default is 200.
            force_field_max_iterations: 1000,
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

