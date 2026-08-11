//! Real (stochastic) distance-geometry conformer embedding — 3D Breakthrough Program,
//! Wave 1, Agent C. See `docs/3d_breakthrough_master_plan.md` §1b/§4 and
//! `docs/etkdg_3d_gap_rfc.md` for the diagnosis this module answers.
//!
//! # What this is
//!
//! `embed_distance_geometry_v2` is a from-scratch orchestration of the bounds/
//! smoothing/Gram/eigendecomposition machinery that already existed, unused, in
//! [`crate::dg_fft`] (`build_bound_matrix`, `smooth_bounds`, `distance_to_gram_matrix`,
//! `jacobi_eigendecompose`, `refine_coords`), plus one piece that was genuinely
//! missing: **stochastic metrization**. `dg_fft::generate_coords_dg` (kept unchanged,
//! still has its own tests) fixes every pairwise target distance to the *midpoint* of
//! its `[lower, upper]` bound — deterministic, and not what "distance geometry" means
//! in the ETKDG/RDKit sense. This module instead *samples* each target distance
//! uniformly at random from within its smoothed bounds, using a caller-controlled seed
//! (`EmbedParameters.random_seed`), so that different seeds produce genuinely different
//! initial geometries (verified empirically in this module's tests) — the actual
//! "stochastic" half of stochastic distance geometry.
//!
//! # What this is NOT (forward-compat placeholders, explicit — see PR body)
//!
//! `EmbedParameters` includes several fields that belong conceptually to later Wave 2
//! work (torsion knowledge = Agent E, force-field bridge = Agent F, stereo constraints
//! = Agent D) so the public API doesn't need a breaking change when those land:
//! `use_exp_torsions`, `use_small_ring_torsions`, `use_macrocycle_torsions`,
//! `use_macrocycle_14_bounds`, `prune_rms_threshold`, `num_threads` are accepted but
//! **currently no-op** — this module never reads them to change behavior. Likewise
//! `EmbedFailureCause::{ConstraintOptimizationFailed, MissingForceFieldParameters,
//! MinimizationFailed}` are reserved variants this module never constructs (no force
//! field or constraint-optimization step runs here at all — the acceptance gate for
//! this PR is measured *before* any such step, see the master plan §4).
//! `enforce_chirality` is NOT one of these no-op placeholders — see the
//! "`enforce_chirality`" section below for what it actually does.
//!
//! # Not wired into the live pipeline
//!
//! Nothing in `etkdg.rs` or `Mol.conformer_ensemble()` calls this module. Per the
//! master plan §1b, that integration is an explicit Wave 2 Coordinator step performed
//! only after this PR, Agent F's force-field bridge, and Agent E's torsion knowledge
//! are all separately merged.
//!
//! # `enforce_chirality` (added 2026-08-11, issue #291/#293)
//!
//! **Important, not obvious from the name**: a pairwise distance matrix is
//! reflection-invariant (a molecule and its mirror image have identical pairwise
//! distances — see `docs/etkdg_3d_gap_rfc.md`'s Phase 3 "Correction" section), so
//! nothing in `build_bound_matrix`/`smooth_bounds`/the MDS embedding step above can
//! ever encode which chirality to prefer. `enforce_chirality` does NOT inject a
//! constraint into that machinery. Instead, per attempt, after the raw MDS/random
//! placement (which operates on real coordinates, unlike the bound matrix, and so
//! *can* carry a chirality sign): apply [`stereo_constraints::repair_stereo`] to
//! fix whatever it can (bridge-eligible substituents only — ring-fused
//! stereocenters cannot be fixed this way, see that module's docs) BEFORE
//! `refine_coords` runs, so the bounds-driven relaxation absorbs whatever bond-
//! length distortion the repair's rigid-subtree translation introduced. Then,
//! after `refine_coords`, independently re-verify the *actual returned* geometry
//! (never trust that refinement preserved the repair) — if any declared element is
//! still `Violated`, this attempt fails with [`EmbedFailureCause::
//! StereoConstraintFailed`] and the existing per-attempt retry loop tries a new
//! seed. Molecules with no declared stereo are unaffected (`repair_stereo`/
//! `verify_stereo` are no-ops on them); `enforce_chirality: false` (the default)
//! is byte-identical to before this change — opt-in only, see `ROADMAP.md`'s
//! v0.14.0 S-tier item 1 for why the wider default-path decision (issue #291's
//! `Ignore`-policy population) was deliberately deferred, not folded in here.
//!
//! For ring-fused stereocenters where `repair_stereo` cannot help, retrying with a
//! new stochastic seed is the only mechanism this provides, and its odds are poor
//! for molecules with several declared centers (each is close to an independent
//! coin flip under a chirality-blind embedder, so naively the per-attempt success
//! rate is on the order of 2^-k for k declared centers) — `max_attempts`
//! exhaustion returning `StereoConstraintFailed` for such molecules is expected,
//! correct, fail-closed behavior, not a bug to work around by raising the attempt
//! count unboundedly.
//!
//! **Measured (2026-08-11, `max_attempts: 8`, the 29 stereo-bearing molecules of
//! `scripts/etkdg_vs_rdkit_gap.py::CORPUS`, `random_seed` swept over 0..5, probed
//! and then removed per this project's own measure-before-claiming convention — see
//! the PR body for the raw per-seed runs). This measures the embedder alone (no UFF
//! minimization, unlike issue #291's `embed_pipeline_v2`-level 18/29 figure — the
//! two numbers are not directly comparable, see the PR body):** 25-27/29
//! (86.2%-93.1%) succeed across the 5 base seeds tested — one draw is not a stable
//! percentage, report the range, not a single decimal. Two failures recur at every
//! base seed tested (5/5):
//! - `testosterone`, `cholesterol` — expected: ring-fused declared stereocenters,
//!   the documented `NoBridgeEligibleSubstituent` case above.
//! - `but2ene_Z` (`C/C=C\C`, a plain acyclic Z-alkene, no ring at all, fails 4/5
//!   base seeds tested) — NOT expected from the design above, and NOT a
//!   retry-odds problem: at a fixed base seed, the *raw* (pre-repair,
//!   `enforce_chirality: false`) embedding is already `Violated` for all of 10
//!   tested derived seeds, deterministically, not ~50/50 as the "coin flip per
//!   center" framing predicts. Isolated further (2026-08-11, follow-up
//!   diagnosis): the raw coordinates *do* vary seed to seed (in-plane spread
//!   ~3.1-3.5 Å / ~1.3-1.5 Å, out-of-plane spread ~0.4-0.9 Å -- real 3D content,
//!   not a planar collapse), yet the sign of the C0-C1=C2-C3 dihedral was
//!   identical across all 10 draws. Contrasted directly against a larger,
//!   more conformationally flexible molecule with the same simple 1-declared-
//!   E/Z-bond shape (`cinnamic_acid_E`, `OC(=O)/C=C/c1ccccc1`, same probe, same
//!   10 seeds): its raw dihedral sign split 5-positive/5-negative -- genuine
//!   per-seed variability, exactly what the retry loop is designed to exploit.
//!   **The measured fact is this contrast** (0/10 sign flips for the small/rigid
//!   molecule vs. 5/10 for the larger/flexible one) -- the *mechanism* behind it
//!   (why the small molecule's sampled distance matrices consistently reconstruct
//!   to the same handedness) is not yet isolated; a plausible but unverified
//!   candidate is that `jacobi_eigendecompose`'s eigenvector-sign outcome is a
//!   deterministic, presumably-continuous function of the input Gram matrix, and
//!   a tightly-bounded (low-conformational-freedom) molecule's sampled distance
//!   matrices across different seeds stay too numerically close to one another
//!   to cross whatever boundary would flip it -- not confirmed, do not treat as
//!   established. **Operational consequence, which *is* established regardless
//!   of mechanism**: for at least this molecule class, `max_attempts` retry is
//!   structurally unable to help, because there is no real per-seed variability
//!   to exploit -- this is a distinct, disclosed gap from both the ring-fused
//!   case and the general reflection-invariance point above; tracked as a
//!   follow-up (issue #285 comment), not fixed here. The remaining
//!   single-base-seed-only
//!   failures (`atorvastatin_fragment` at one seed, `cinnamic_acid_E` at another)
//!   match the expected `max_attempts`-retry-loop variance already documented above
//!   (occasional exhaustion of 8 attempts for a molecule the mechanism *can* fix),
//!   not a new distinct cause.

use std::collections::HashMap;

use crate::clock::Instant;

use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule};

use crate::coords::{Coords3D, Point3};
use crate::dg_fft::{
    DG_MAX_ATOMS, build_bound_matrix, center_coordinates, distance_to_gram_matrix,
    jacobi_eigendecompose, refine_coords, smooth_bounds,
};
use crate::prng::Prng;
use crate::stereo_constraints::{repair_stereo, verify_stereo};

/// Number of bounds-driven SHAKE-like refinement passes after the initial MDS/random
/// placement. Higher than `dg_fft::generate_coords_dg`'s 300 because a stochastic
/// (not midpoint) initial distance draw starts further from a self-consistent geometry
/// on average and needs more passes to converge onto the frozen-58 gate.
///
/// ponytail: fixed constant, not exposed as a parameter — tune here if a larger/harder
/// corpus needs it; not worth a config field until something actually needs to vary it.
const REFINE_ITERS: usize = 800;

/// Any Gram-matrix eigenvalue below this (in a matrix whose scale is Å², i.e. tens to
/// low hundreds) is treated as genuinely negative rather than floating-point noise
/// around zero.
const NEGATIVE_EIGENVALUE_EPS: f64 = 1e-6;

/// Tolerance for the smoothing invariant check (Å) — accounts for floating-point
/// accumulation across an O(n) chain of additions/subtractions in Floyd-Warshall, not
/// a relaxation of the invariant itself.
const INVARIANT_EPS: f64 = 1e-6;

/// A bonded pair placed closer than this (Å) after refinement is a degenerate/
/// coincident placement, not a valid geometry, regardless of what the bounds said.
const MIN_ACCEPTABLE_BOND_DIST: f64 = 0.1;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parameters controlling [`embed_distance_geometry_v2`].
///
/// See the module docs for which fields this implementation actually consumes vs.
/// accepts-but-does-not-yet-act-on (forward compatibility for later waves).
#[derive(Debug, Clone)]
pub struct EmbedParameters {
    /// Seed controlling every random draw this call makes (stochastic metrization,
    /// and the random-coordinate fallback if `use_random_coords` is set). The same
    /// seed on the same target/thread-count reproduces the same output — this
    /// codebase's existing reproducibility convention (see `prng.rs`), not a
    /// cross-platform bit-exactness claim. **Consumed.**
    pub random_seed: u64,
    /// Maximum number of independent embedding attempts (each with its own derived
    /// seed) before giving up. A bad stochastic draw is retried, not fatal by itself.
    /// **Consumed.**
    pub max_attempts: usize,
    /// Wall-clock budget in milliseconds across all attempts combined. Checked
    /// between attempts (not preemptively inside one attempt). `None` = no limit.
    /// **Consumed.**
    pub timeout_ms: Option<u64>,
    /// When true, skip bounds-sampled MDS entirely and start from uniform-random
    /// coordinates in a box sized to the molecule, then run the same bounds-driven
    /// refinement pass. Mirrors ETKDG's own `useRandomCoords` fallback semantics.
    /// **Consumed.**
    pub use_random_coords: bool,
    /// When true, every embedding attempt is checked against the molecule's
    /// declared tetrahedral (`@`/`@@`) and double-bond (`/`/`\`) stereo after
    /// refinement; violations are repaired where possible (bridge-eligible
    /// substituents) before the check, and an attempt that still violates any
    /// declared element after that is treated as failed, so the retry loop tries a
    /// new seed. Exhausting `max_attempts` without a fully-satisfying geometry
    /// returns `EmbedFailureCause::StereoConstraintFailed` (fail-closed — never
    /// silently returns a geometry that violates declared stereo). See the module
    /// doc's "`enforce_chirality`" section for why this can't be a bound-matrix
    /// constraint and what it can/cannot fix. Molecules with **no** declared stereo
    /// are unaffected (zero extra cost). `false` (default) is byte-identical to
    /// this field's behavior before 2026-08-11. **Consumed.**
    pub enforce_chirality: bool,
    /// Reserved for Agent E's experimental-torsion-preference integration (Wave 2).
    /// **Not consumed** — accepted for forward API compatibility only.
    pub use_exp_torsions: bool,
    /// Reserved for Agent E's small-ring torsion handling (Wave 2). **Not consumed.**
    pub use_small_ring_torsions: bool,
    /// Reserved for Agent E's macrocycle torsion handling (Wave 2). **Not consumed.**
    pub use_macrocycle_torsions: bool,
    /// Reserved for macrocycle-specific 1-4 bound relaxation (Wave 2). **Not consumed.**
    pub use_macrocycle_14_bounds: bool,
    /// Reserved for ensemble-level RMSD pruning once this feeds `conformer.rs`'s
    /// multi-conformer path (Wave 2 integration). A single-embed function has nothing
    /// to prune against. **Not consumed.**
    pub prune_rms_threshold: Option<f64>,
    /// Reserved for parallelizing independent attempts across threads. This
    /// implementation runs attempts sequentially regardless of this value.
    /// **Not consumed.**
    pub num_threads: usize,
    /// When true, [`EmbedStats::failure_counts`] accumulates real per-cause counts
    /// across every failed attempt. When false, failures still count toward
    /// `attempts_used` but the per-cause breakdown is left empty (cheaper bookkeeping,
    /// not a lossy summary of data that was collected anyway). **Consumed.**
    pub track_failures: bool,
}

impl Default for EmbedParameters {
    fn default() -> Self {
        Self {
            random_seed: 0xC0FF_EE42_D157_6E02,
            max_attempts: 8,
            timeout_ms: None,
            use_random_coords: false,
            enforce_chirality: false,
            use_exp_torsions: false,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: false,
            use_macrocycle_14_bounds: false,
            prune_rms_threshold: None,
            num_threads: 1,
            track_failures: false,
        }
    }
}

/// Why a single embedding attempt (or the whole call, after retries) failed.
///
/// See the module docs for which variants this implementation can actually return
/// vs. reserved-for-later-waves placeholders.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EmbedFailureCause {
    /// The molecule itself can't be embedded: contains a wildcard (`*`) atom with no
    /// real element to derive bond lengths from, or has zero atoms where a non-empty
    /// embedding was expected.
    InvalidTopology,
    /// A structurally inconsistent bound was produced for an actually-bonded pair
    /// before smoothing even ran (e.g. `lower > upper`) — a bug in bounds
    /// construction, not a smoothing or sampling issue.
    BoundsConstructionFailed,
    /// Triangle-inequality smoothing produced bounds that violate the invariants this
    /// module checks (loosened rather than tightened a bound, or produced
    /// `lower > upper`). See [`EmbedStats::last_smoothing_invariants_ok`].
    BoundsSmoothingFailed,
    /// Reserved: the sampled/averaged distance matrix was detected as unrecoverably
    /// non-Euclidean (this implementation currently only surfaces this as a NaN/Inf
    /// coordinate after refinement, not as an upfront distance-matrix check).
    NonEuclideanDistanceMatrix,
    /// Classical MDS produced no usable positive eigenvalue (or refinement produced
    /// non-finite or degenerate/coincident coordinates for a bonded pair).
    EigenEmbeddingFailed,
    /// Reserved for Wave 2 force-field/constraint integration. Not constructed here.
    ConstraintOptimizationFailed,
    /// `enforce_chirality` was requested and, after genuinely trying every attempt
    /// (repair where possible, retry with a new seed otherwise), no attempt
    /// produced a geometry satisfying every declared tetrahedral/E-Z stereo
    /// element. Expected, correct outcome for molecules whose declared stereo
    /// includes at least one ring-fused center `repair_stereo` cannot fix — not
    /// itself evidence of a bug.
    StereoConstraintFailed,
    /// Reserved for Wave 2 force-field integration (Agent F). Not constructed here —
    /// this module never runs a force field.
    MissingForceFieldParameters,
    /// Reserved for Wave 2 force-field integration (Agent F). Not constructed here —
    /// this module never runs a minimizer (the acceptance gate is measured on raw,
    /// pre-minimization output by design).
    MinimizationFailed,
    /// The wall-clock budget (`EmbedParameters.timeout_ms`) was exceeded before an
    /// attempt succeeded.
    Timeout,
    /// `mol.atom_count()` exceeds [`crate::dg_fft::DG_MAX_ATOMS`] (currently 500):
    /// O(n²) bound-matrix memory and O(n³) Floyd-Warshall smoothing make larger
    /// molecules prohibitive for this implementation.
    AtomLimitExceeded,
    /// The molecule contains a wildcard atom (see `InvalidTopology`) — kept as a
    /// distinct variant name per the required API shape; this implementation returns
    /// `InvalidTopology` for wildcards today (both describe the same root cause here;
    /// kept separate in the enum for a future, more granular per-element check).
    UnsupportedElement,
}

/// Diagnostics accumulated across every attempt of one `embed_distance_geometry_v2*`
/// call. Always available — on success **and** on final failure (see
/// [`embed_distance_geometry_v2_detail`]'s return type) — never silently dropped.
#[derive(Debug, Clone, Default)]
pub struct EmbedStats {
    /// Number of embedding attempts actually made (`1..=max_attempts`).
    pub attempts_used: usize,
    /// Per-cause failure counts across every failed attempt. Only populated when
    /// `EmbedParameters.track_failures` is true (see that field's docs); a genuine
    /// count map, never just a boolean "some attempts failed."
    pub failure_counts: HashMap<EmbedFailureCause, usize>,
    /// How many Gram-matrix eigenvalues (out of all `n`, not just the 3 used for
    /// coordinates) came out below `-1e-6` on the attempt that produced the returned
    /// result (or the last attempt made, if every attempt failed). A distance matrix
    /// built from approximate (sampled-within-bounds, not exactly Euclidean) distances
    /// routinely has some; this is reported for visibility rather than silently
    /// truncated away.
    pub negative_eigenvalues_beyond_embedding_dim: usize,
    /// Magnitude of the most negative such eigenvalue (0.0 if none).
    pub max_negative_eigenvalue_magnitude: f64,
    /// Whether the last attempt's smoothing invariant check passed (see
    /// [`EmbedFailureCause::BoundsSmoothingFailed`]).
    pub last_smoothing_invariants_ok: bool,
    /// Whether the attempt that produced the returned result used the
    /// `use_random_coords` initial-placement path instead of bounds-sampled MDS.
    pub used_random_coords: bool,
    /// How many caller-supplied [`DistanceBoundAdjustment`]s were actually written
    /// into the bound matrix before smoothing, on the attempt that produced the
    /// returned result. `0` for every call through the public
    /// [`embed_distance_geometry_v2`]/[`embed_distance_geometry_v2_detail`] API
    /// (which always passes an empty adjustment slice) — only Wave 2/3
    /// Coordinator integration (`pipeline_v2.rs`'s
    /// `embed_distance_geometry_v2_with_adjustments`) ever sets this above 0.
    pub adjustments_applied: usize,
    /// Whether `enforce_chirality`'s repair-before-refine step actually ran on the
    /// attempt that produced the returned result (i.e. `enforce_chirality` was set
    /// AND the molecule had declared stereo). `false` whenever `enforce_chirality`
    /// is `false` or the molecule declares no stereo, regardless of `enforce_chirality`.
    pub stereo_repair_attempted: bool,
}

fn record_failure(stats: &mut EmbedStats, params: &EmbedParameters, cause: EmbedFailureCause) {
    if params.track_failures {
        *stats.failure_counts.entry(cause).or_insert(0) += 1;
    }
}

/// Embed `mol` into 3D coordinates via stochastic distance geometry: bound-matrix
/// construction → triangle-inequality smoothing → stochastic metrization → classical
/// MDS (Gram matrix + Jacobi eigendecomposition) → bounds-driven refinement.
///
/// No force-field minimization runs here (see module docs) — this is the raw embedder
/// output the 3D Breakthrough Program's Wave 1 acceptance gate measures directly.
///
/// Discards the accumulated [`EmbedStats`]; use [`embed_distance_geometry_v2_detail`]
/// to see attempt counts, failure breakdowns, and eigenvalue diagnostics.
pub fn embed_distance_geometry_v2(
    mol: &Molecule,
    params: &EmbedParameters,
) -> Result<Coords3D, EmbedFailureCause> {
    match embed_distance_geometry_v2_detail(mol, params) {
        Ok((coords, _stats)) => Ok(coords),
        Err((cause, _stats)) => Err(cause),
    }
}

/// Same as [`embed_distance_geometry_v2`] but always returns [`EmbedStats`] alongside
/// the result — on success **and** on exhaustion of `max_attempts` — so
/// `track_failures` data and eigenvalue diagnostics are never lost on the failure path.
///
/// A thin wrapper over [`embed_distance_geometry_v2_with_adjustments`] with an empty
/// adjustment slice — this is what guarantees (structurally, not by convention) that
/// the public raw API's output with zero adjustments is byte-identical to before Wave
/// 2/3 Coordinator integration added that internal hook: there is exactly one control
/// flow, not two kept in sync by hand.
pub fn embed_distance_geometry_v2_detail(
    mol: &Molecule,
    params: &EmbedParameters,
) -> Result<(Coords3D, EmbedStats), (EmbedFailureCause, EmbedStats)> {
    match embed_distance_geometry_v2_with_adjustments(mol, params, &[]) {
        Ok(v) => Ok(v),
        Err((EmbedWithAdjustmentsFailure::Embed(cause), stats)) => Err((cause, stats)),
        Err((EmbedWithAdjustmentsFailure::InvalidAdjustment, _stats)) => {
            unreachable!(
                "no adjustments were passed by this wrapper, so InvalidAdjustment cannot occur"
            )
        }
    }
}

/// One proposed override of a single atom pair's `[lower, upper]` distance bound,
/// applied to the bound matrix immediately after [`build_bound_matrix`] and before
/// [`smooth_bounds`] runs (never the reverse — see
/// [`embed_distance_geometry_v2_with_adjustments`]'s doc). `pub(crate)`: this is a
/// Wave 2/3 Coordinator integration seam (`pipeline_v2.rs`, carrying Agent E's
/// `etkdg_knowledge::PairBoundAdjustment` — e.g. macrocycle 1-4 relaxation — into
/// Agent C's embedder), not part of this module's own public API.
#[derive(Debug, Clone, Copy)]
pub(crate) struct DistanceBoundAdjustment {
    pub atom1: AtomIdx,
    pub atom2: AtomIdx,
    pub lower: f64,
    pub upper: f64,
}

/// Failure mode for [`embed_distance_geometry_v2_with_adjustments`]: either a
/// malformed adjustment (checked once, up front, independent of any embedding
/// attempt/seed) or an ordinary embedding failure (same causes as the public API).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EmbedWithAdjustmentsFailure {
    /// An adjustment named an out-of-range or identical atom pair, or had a
    /// non-finite (`NaN`/`Inf`) bound, or `lower > upper`. Caught before any bound
    /// matrix is even touched, so a bad adjustment can never masquerade as
    /// [`EmbedFailureCause::BoundsConstructionFailed`] (that variant means the
    /// pre-existing bounds machinery itself produced an inconsistent bond bound,
    /// a structurally different problem from a caller-supplied override being bad).
    InvalidAdjustment,
    /// Same causes [`embed_distance_geometry_v2_detail`] can return.
    Embed(EmbedFailureCause),
}

/// Same as [`embed_distance_geometry_v2_detail`], but lets a caller (Wave 2/3
/// Coordinator's `pipeline_v2.rs`) inject pre-smoothing bound overrides — e.g. Agent
/// E's `macrocycle_14_bound_adjustments()` output — into the bound matrix this
/// module builds internally, without duplicating any of `dg_fft`'s bounds/smoothing/
/// Gram/eigendecomposition machinery.
///
/// Every adjustment is validated **once, up front, independent of any embedding
/// attempt** (index range, `atom1 != atom2`, both bounds finite, `lower <= upper`) —
/// a bad adjustment fails closed with [`EmbedWithAdjustmentsFailure::InvalidAdjustment`]
/// before any bound matrix is built, never partially applied.
///
/// Adjustments are written into the bound matrix immediately after
/// [`build_bound_matrix`] runs and *before* [`smooth_bounds`] runs, every attempt (the
/// matrix is rebuilt fresh per attempt, so the override must be re-applied each time
/// too). Critically, the smoothing-invariant check
/// ([`smoothing_preserves_invariants`]) compares smoothed bounds against the
/// **post-adjustment** matrix as its baseline, not the pre-adjustment one: the
/// invariant is a claim about what *smoothing* is allowed to do (only ever tighten),
/// and a caller-requested adjustment (e.g. a macrocycle 1-4 band widening a naive
/// single-trans-configuration pin) is a deliberate, disclosed relaxation that happens
/// *before* smoothing sees the matrix at all -- comparing against the pre-adjustment
/// matrix would make every such widening a spurious `BoundsSmoothingFailed`.
///
/// With `adjustments = &[]`, this is byte-identical to
/// [`embed_distance_geometry_v2_detail`]'s pre-existing behavior (in fact, that
/// function now delegates here) — the empty-adjustment case never touches the bound
/// matrix at all, so there is no drift to keep in sync.
pub(crate) fn embed_distance_geometry_v2_with_adjustments(
    mol: &Molecule,
    params: &EmbedParameters,
    adjustments: &[DistanceBoundAdjustment],
) -> Result<(Coords3D, EmbedStats), (EmbedWithAdjustmentsFailure, EmbedStats)> {
    let mut stats = EmbedStats::default();

    let n = mol.atom_count();
    for adj in adjustments {
        let i = adj.atom1.0 as usize;
        let j = adj.atom2.0 as usize;
        let well_formed = i < n
            && j < n
            && i != j
            && adj.lower.is_finite()
            && adj.upper.is_finite()
            && adj.lower >= 0.0
            && adj.lower <= adj.upper;
        if !well_formed {
            return Err((EmbedWithAdjustmentsFailure::InvalidAdjustment, stats));
        }
    }

    if mol_has_wildcard_atom(mol) {
        return Err((
            EmbedWithAdjustmentsFailure::Embed(EmbedFailureCause::InvalidTopology),
            stats,
        ));
    }
    // `enforce_chirality`'s real check-repair-or-retry logic runs per attempt
    // inside `try_embed_once` (see the module doc's "`enforce_chirality`" section)
    // -- there is no upfront fail-closed refusal anymore, since this module now
    // genuinely tries before giving up.

    if n > DG_MAX_ATOMS {
        return Err((
            EmbedWithAdjustmentsFailure::Embed(EmbedFailureCause::AtomLimitExceeded),
            stats,
        ));
    }
    if n == 0 {
        return Ok((Coords3D::new_zeroed(0), stats));
    }
    if n == 1 {
        let mut coords = Coords3D::new_zeroed(1);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        stats.attempts_used = 1;
        return Ok((coords, stats));
    }

    let start = Instant::now();
    let max_attempts = params.max_attempts.max(1);
    let mut last_cause = EmbedFailureCause::EigenEmbeddingFailed;

    for attempt in 0..max_attempts {
        stats.attempts_used = attempt + 1;

        if let Some(budget_ms) = params.timeout_ms
            && start.elapsed().as_millis() as u64 > budget_ms
        {
            last_cause = EmbedFailureCause::Timeout;
            record_failure(&mut stats, params, last_cause);
            break;
        }

        let attempt_seed = derive_attempt_seed(params.random_seed, attempt);
        match try_embed_once(mol, params, attempt_seed, &mut stats, adjustments) {
            Ok(coords) => return Ok((coords, stats)),
            Err(cause) => {
                last_cause = cause;
                record_failure(&mut stats, params, cause);
            }
        }
    }

    Err((EmbedWithAdjustmentsFailure::Embed(last_cause), stats))
}

/// How well does a **final, already-embedded** geometry satisfy the distance bounds
/// it was built from? Diagnostic, not part of the embedding call itself or the
/// acceptance gate: recomputes the same bound matrix `embed_distance_geometry_v2`
/// used internally (bounds construction is deterministic given `mol`, independent of
/// the random seed) and checks **every** atom pair -- not just bonded pairs -- against
/// `[lower, upper]`.
///
/// `refine_coords`'s SHAKE-like projection is a heuristic bounds-satisfaction pass,
/// not a guaranteed-convergent solver, so nonzero residual violation is expected;
/// this makes it visible and measurable rather than assumed away.
pub fn bounds_conformance(mol: &Molecule, coords: &Coords3D) -> BoundsConformance {
    let (lower0, upper0) = build_bound_matrix(mol);
    let mut lower = lower0;
    let mut upper = upper0;
    smooth_bounds(&mut lower, &mut upper);

    let n = mol.atom_count();
    let mut n_pairs = 0usize;
    let mut n_violations = 0usize;
    let mut max_rel_violation = 0.0_f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = coords
                .get(AtomIdx(i as u32))
                .distance(&coords.get(AtomIdx(j as u32)));
            let lo = lower[i][j];
            let hi = upper[i][j];
            n_pairs += 1;
            let rel = if d < lo - 1e-6 {
                Some((lo - d) / lo.max(1e-9))
            } else if hi.is_finite() && d > hi + 1e-6 {
                Some((d - hi) / hi.max(1e-9))
            } else {
                None
            };
            if let Some(rel) = rel {
                n_violations += 1;
                if rel > max_rel_violation {
                    max_rel_violation = rel;
                }
            }
        }
    }
    BoundsConformance {
        n_pairs,
        n_violations,
        max_rel_violation,
    }
}

/// Result of [`bounds_conformance`]: how many of `mol`'s atom pairs (all pairs, not
/// just bonded ones) land outside their own smoothed `[lower, upper]` bounds in a
/// given final geometry, and the worst relative violation among them.
#[derive(Debug, Clone, Copy, Default)]
pub struct BoundsConformance {
    pub n_pairs: usize,
    pub n_violations: usize,
    pub max_rel_violation: f64,
}

// ---------------------------------------------------------------------------
// Internals
// ---------------------------------------------------------------------------

fn mol_has_wildcard_atom(mol: &Molecule) -> bool {
    (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).wildcard)
}

fn mol_has_declared_stereo(mol: &Molecule) -> bool {
    if (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None) {
        return true;
    }
    mol.bonds()
        .any(|(_, bond)| matches!(bond.order, BondOrder::Up | BondOrder::Down))
}

/// Derive a per-attempt seed from the caller's base seed. Deterministic: the same
/// `(base, attempt)` pair always produces the same derived seed, so a fixed
/// `random_seed` reproduces the exact same sequence of attempts (and thus the exact
/// same final result) every time.
fn derive_attempt_seed(base: u64, attempt: usize) -> u64 {
    const GOLDEN: u64 = 0x9E37_79B9_7F4A_7C15;
    const OFFSET: u64 = 0xD1B5_4A32_D192_ED03;
    base ^ (attempt as u64).wrapping_mul(GOLDEN).wrapping_add(OFFSET)
}

fn try_embed_once(
    mol: &Molecule,
    params: &EmbedParameters,
    seed: u64,
    stats: &mut EmbedStats,
    adjustments: &[DistanceBoundAdjustment],
) -> Result<Coords3D, EmbedFailureCause> {
    let n = mol.atom_count();

    let (mut lower0, mut upper0) = build_bound_matrix(mol);

    // Caller-supplied overrides (e.g. Agent E's macrocycle 1-4 relaxation), already
    // validated as well-formed by `embed_distance_geometry_v2_with_adjustments`
    // before any attempt started -- write them into THIS attempt's freshly-built
    // matrix now, before the bonded-pair sanity check and smoothing below, so both
    // see the adjusted values as the baseline. Empty for every public-API call
    // ([`embed_distance_geometry_v2_detail`] always passes `&[]`), so this loop is a
    // no-op there.
    for adj in adjustments {
        let i = adj.atom1.0 as usize;
        let j = adj.atom2.0 as usize;
        lower0[i][j] = adj.lower;
        lower0[j][i] = adj.lower;
        upper0[i][j] = adj.upper;
        upper0[j][i] = adj.upper;
    }
    if !adjustments.is_empty() {
        stats.adjustments_applied = adjustments.len();
    }

    // A genuinely bonded pair with lower > upper before smoothing is a bounds-
    // construction bug, not a smoothing artifact -- catch it before smoothing muddies
    // the picture. Runs against the (possibly adjusted) `lower0`/`upper0` above, but
    // in practice adjustments only ever target genuine 1-4 (non-bonded) pairs (see
    // `bounds14.rs`'s own genuine-1-4 guard), so this never fires because of an
    // adjustment.
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        if lower0[i][j] > upper0[i][j] + INVARIANT_EPS {
            return Err(EmbedFailureCause::BoundsConstructionFailed);
        }
    }

    let mut lower = lower0.clone();
    let mut upper = upper0.clone();
    smooth_bounds(&mut lower, &mut upper);

    // Coordinator-requested invariant: smoothing must only TIGHTEN bounds (upper
    // never increases, lower never decreases relative to pre-smoothing) and must
    // never invert lower > upper.
    let invariants_ok = smoothing_preserves_invariants(&lower0, &upper0, &lower, &upper);
    stats.last_smoothing_invariants_ok = invariants_ok;
    if !invariants_ok {
        return Err(EmbedFailureCause::BoundsSmoothingFailed);
    }

    let max_finite_upper: f64 = upper
        .iter()
        .flat_map(|row| row.iter())
        .filter(|v| v.is_finite())
        .cloned()
        .fold(0.0_f64, f64::max);
    let fallback = (max_finite_upper * 4.0).max(10.0);

    let mut coords = if params.use_random_coords {
        stats.used_random_coords = true;
        random_initial_coords(n, fallback, seed)
    } else {
        stats.used_random_coords = false;
        let dist_matrix = sample_distance_matrix(n, &lower, &upper, fallback, seed);
        mds_embed(&dist_matrix, stats)?
    };

    center_coordinates(&mut coords);

    // Repair BEFORE refinement, not after: `repair_stereo`'s rigid-subtree
    // translation preserves bond lengths *within* the moved substituent exactly,
    // but can and does push its distances to the *rest* of the molecule outside
    // their smoothed bounds (measured empirically while designing this: e.g.
    // L-alanine's bounds_conformance violations went 1->4 when repair was applied
    // AFTER refine_coords). Applying it here lets the bounds-driven relaxation
    // below absorb that distortion the same way it absorbs everything else -- see
    // the module doc's "`enforce_chirality`" section.
    if params.enforce_chirality && mol_has_declared_stereo(mol) {
        stats.stereo_repair_attempted = true;
        coords = match repair_stereo(mol, &coords) {
            Ok(outcome) => outcome.coords,
            Err(failure) => failure.partial_coords,
        };
    }

    refine_coords(&mut coords, &lower, &upper, REFINE_ITERS);

    validate_final_coords(mol, &coords)?;

    // Never trust that refinement preserved the repair above -- it's a bounds-
    // driven relaxation with no chirality awareness of its own, so re-verify the
    // *actual returned* geometry, not the pre-refinement one.
    if params.enforce_chirality
        && mol_has_declared_stereo(mol)
        && !verify_stereo(mol, &coords).is_fully_satisfied()
    {
        return Err(EmbedFailureCause::StereoConstraintFailed);
    }

    Ok(coords)
}

/// Sample each pairwise target distance uniformly at random from within its smoothed
/// `[lower, upper]` interval, using a seeded stream (deterministic given `seed`).
///
/// ponytail: independent per-pair draws (not sequential/partial re-metrization with
/// re-smoothing after each draw) can, on paper, violate the triangle inequality for
/// some triples. This module treats the sampled matrix as an MDS *initial guess
/// only* — `refine_coords`'s bounds-driven projection (called by every caller of this
/// function) is what actually enforces the geometry, matching `dg_fft.rs`'s own
/// existing doc comment that classical MDS is "an initial guess only." Upgrade path:
/// sequential metrization with re-smoothing per draw, if a future corpus needs
/// tighter raw (pre-refinement) geometry than this gives.
fn sample_distance_matrix(
    n: usize,
    lower: &[Vec<f64>],
    upper: &[Vec<f64>],
    fallback: f64,
    seed: u64,
) -> Vec<Vec<f64>> {
    let mut prng = Prng::from_seed(seed);
    let mut dist = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let hi_bound = if upper[i][j].is_finite() {
                upper[i][j]
            } else {
                fallback
            };
            let lo_bound = lower[i][j].min(hi_bound);
            let hi_bound = hi_bound.max(lo_bound);
            let d = if hi_bound > lo_bound {
                lo_bound + (hi_bound - lo_bound) * prng.f64()
            } else {
                lo_bound
            };
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }
    dist
}

fn random_initial_coords(n: usize, radius_scale: f64, seed: u64) -> Coords3D {
    // Distinct sub-stream from the metrization sampler so `use_random_coords = true`
    // doesn't just replay the same bytes a bounds-sampling attempt would have used.
    let mut prng = Prng::from_seed(seed ^ 0xA5A5_A5A5_A5A5_A5A5);
    let box_size = radius_scale.max(3.0) * (n as f64).cbrt();
    let mut coords = Coords3D::new_zeroed(n);
    for i in 0..n {
        let x = (prng.f64() - 0.5) * 2.0 * box_size;
        let y = (prng.f64() - 0.5) * 2.0 * box_size;
        let z = (prng.f64() - 0.5) * 2.0 * box_size;
        coords.set(AtomIdx(i as u32), Point3::new(x, y, z));
    }
    coords
}

/// Classical MDS: distance matrix → Gram matrix → Jacobi eigendecomposition → top-3
/// positive eigenvalue/eigenvector coordinates. Reports (never silently truncates)
/// negative eigenvalues found beyond the 3 used for coordinates.
fn mds_embed(
    dist_matrix: &[Vec<f64>],
    stats: &mut EmbedStats,
) -> Result<Coords3D, EmbedFailureCause> {
    let n = dist_matrix.len();
    let gram = distance_to_gram_matrix(dist_matrix);
    let (eigenvalues, eigenvectors) = jacobi_eigendecompose(&gram);

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        eigenvalues[b]
            .partial_cmp(&eigenvalues[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let use_indices: Vec<usize> = order
        .iter()
        .copied()
        .filter(|&i| eigenvalues[i] > 1e-10)
        .take(3)
        .collect();

    if use_indices.is_empty() {
        return Err(EmbedFailureCause::EigenEmbeddingFailed);
    }

    // Negative-eigenvalue visibility: every eigenvalue below -EPS, not just ones
    // beyond the top 3 -- a sampled (non-exactly-Euclidean) distance matrix routinely
    // has some, and this must be visible in EmbedStats rather than silently dropped.
    let mut neg_count = 0usize;
    let mut max_neg_mag = 0.0_f64;
    for &ev in &eigenvalues {
        if ev < -NEGATIVE_EIGENVALUE_EPS {
            neg_count += 1;
            if -ev > max_neg_mag {
                max_neg_mag = -ev;
            }
        }
    }
    stats.negative_eigenvalues_beyond_embedding_dim = neg_count;
    stats.max_negative_eigenvalue_magnitude = max_neg_mag;

    let mut coords = Coords3D::new_zeroed(n);
    for i in 0..n {
        let coord = |dim: usize| -> f64 {
            use_indices
                .get(dim)
                .map_or(0.0, |&idx| eigenvalues[idx].sqrt() * eigenvectors[i][idx])
        };
        coords.set(AtomIdx(i as u32), Point3::new(coord(0), coord(1), coord(2)));
    }
    Ok(coords)
}

/// Verify the Coordinator-requested smoothing invariants: bounds only ever tighten
/// (never loosen), and smoothing never inverts `lower <= upper` into `lower > upper`.
fn smoothing_preserves_invariants(
    lower0: &[Vec<f64>],
    upper0: &[Vec<f64>],
    lower: &[Vec<f64>],
    upper: &[Vec<f64>],
) -> bool {
    let n = lower0.len();
    for i in 0..n {
        for j in 0..n {
            if upper[i][j] > upper0[i][j] + INVARIANT_EPS {
                return false; // smoothing must never LOOSEN an upper bound
            }
            if lower[i][j] < lower0[i][j] - INVARIANT_EPS {
                return false; // smoothing must never LOOSEN a lower bound
            }
            if lower[i][j] > upper[i][j] + INVARIANT_EPS {
                return false; // smoothing must never invert lower > upper
            }
        }
    }
    true
}

/// Reject non-finite coordinates and degenerate/coincident bonded-atom placements.
/// Never returns `Ok` for an all-zero or degenerate geometry -- there is no
/// zero-coordinates-on-failure path in this module.
fn validate_final_coords(mol: &Molecule, coords: &Coords3D) -> Result<(), EmbedFailureCause> {
    for i in 0..coords.atom_count() {
        let p = coords.get(AtomIdx(i as u32));
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            return Err(EmbedFailureCause::NonEuclideanDistanceMatrix);
        }
    }
    for (_, bond) in mol.bonds() {
        let d = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
        if d < MIN_ACCEPTABLE_BOND_DIST {
            return Err(EmbedFailureCause::EigenEmbeddingFailed);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn worst_bond(mol: &Molecule, coords: &Coords3D) -> f64 {
        mol.bonds()
            .map(|(_, bond)| coords.get(bond.atom1).distance(&coords.get(bond.atom2)))
            .fold(0.0_f64, f64::max)
    }

    /// KNOWN LIMITATION (documented, not fixed by this module -- see PR body):
    /// every 3-membered ring fails closed with `BoundsConstructionFailed`. Root
    /// cause is in `dg_fft::build_bound_matrix` (not this module's code): its
    /// angle-constraint loop treats every pair of a center atom's neighbors as a
    /// 1-3 (through-center) relationship and unconditionally tightens their bound
    /// using the *generic* ideal angle (~109.5°/120°, `ideal_bond_angle` has no
    /// notion of ring strain) -- but in a 3-membered ring, that "1-3" pair is
    /// *also* a direct 1-2 bonded pair (the ring closes one bond away), so the
    /// angle constraint's generic-angle-derived bound overwrites the correct,
    /// much shorter bond-length bound with a value inconsistent with it
    /// (concretely, for cyclopropane's ring-closing C-C pair: bond constraint
    /// gives upper ≈ 1.59 Å, the angle constraint then tightens lower to ≈ 2.41 Å
    /// using the generic ~109.5° angle instead of the real ~60° ring angle,
    /// producing `lower > upper` for the same pair). This module's own
    /// pre-smoothing sanity check (`try_embed_once`) catches exactly this and
    /// fails closed with a typed error -- correct behavior, just a real, disclosed
    /// gap for a common drug-discovery motif (cyclopropane/epoxide/aziridine
    /// rings), not silently mishandled.
    #[test]
    fn three_membered_rings_fail_closed_not_silently() {
        for smiles in [
            "C1CC1",         // cyclopropane
            "C1CO1",         // epoxide
            "C1CN1",         // aziridine
            "C1CS1",         // thiirane
            "C1CC1c1ccccc1", // cyclopropylbenzene
            "C1CC1C(=O)O",   // cyclopropanecarboxylic acid
        ] {
            let mol = parse(smiles).unwrap();
            let params = EmbedParameters::default();
            let err = embed_distance_geometry_v2(&mol, &params).unwrap_err();
            assert_eq!(
                err,
                EmbedFailureCause::BoundsConstructionFailed,
                "{smiles}: expected a typed BoundsConstructionFailed, not a silent success or a different error"
            );
        }
        // Controls: 4- and 5-membered rings are unaffected (this is specific to
        // 3-membered rings, not "any small ring").
        for smiles in ["C1CCC1", "C1CCCC1"] {
            let mol = parse(smiles).unwrap();
            let params = EmbedParameters::default();
            assert!(
                embed_distance_geometry_v2(&mol, &params).is_ok(),
                "{smiles}: 4/5-membered rings should still embed fine"
            );
        }
    }

    #[test]
    fn cyclopropane_exact_bound_contradiction_verified() {
        // Verifies the exact numbers cited in the doc comment above (and in the PR
        // body's Known Limitations) rather than trusting them from memory.
        let mol = parse("C1CC1").unwrap();
        let (lower, upper) = crate::dg_fft::build_bound_matrix(&mol);
        let mut found = false;
        for (_, bond) in mol.bonds() {
            let i = bond.atom1.0 as usize;
            let j = bond.atom2.0 as usize;
            if lower[i][j] > upper[i][j] {
                println!(
                    "cyclopropane ring-closing pair ({i},{j}): lower={:.3} upper={:.3}",
                    lower[i][j], upper[i][j]
                );
                assert!((lower[i][j] - 2.414).abs() < 0.01, "lower={}", lower[i][j]);
                assert!((upper[i][j] - 1.590).abs() < 0.01, "upper={}", upper[i][j]);
                found = true;
            }
        }
        assert!(
            found,
            "expected at least one bonded pair with lower > upper in cyclopropane's bound matrix"
        );
    }

    #[test]
    fn ethane_bond_length_reasonable() {
        let mol = parse("CC").unwrap();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).expect("ethane should embed");
        let d = coords.get(AtomIdx(0)).distance(&coords.get(AtomIdx(1)));
        assert!((d - 1.54).abs() < 0.2, "C-C distance {d}");
    }

    #[test]
    fn benzene_ring_closes_and_bonds_reasonable() {
        let mol = parse("c1ccccc1").unwrap();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).expect("benzene should embed");
        let worst = worst_bond(&mol, &coords);
        assert!(worst < 2.0, "benzene worst bond {worst}");
    }

    #[test]
    fn same_seed_reproducible() {
        let mol = parse("CCCCCCCCCC").unwrap(); // decane
        let params = EmbedParameters {
            random_seed: 12345,
            ..EmbedParameters::default()
        };
        let c1 = embed_distance_geometry_v2(&mol, &params).unwrap();
        let c2 = embed_distance_geometry_v2(&mol, &params).unwrap();
        for i in 0..mol.atom_count() {
            let p1 = c1.get(AtomIdx(i as u32));
            let p2 = c2.get(AtomIdx(i as u32));
            assert_eq!(
                p1, p2,
                "same seed must reproduce identical coords at atom {i}"
            );
        }
    }

    #[test]
    fn different_seed_gives_different_output() {
        let mol = parse("CCCCCCCCCC").unwrap(); // decane
        // Seeds 0/1 specifically -- under the old buggy `Prng::from_seed`
        // (`Self(seed | 1)`), 0 and 1 both mapped to the identical internal
        // state, so this exact pair is what would have caught that
        // regression. Seeds 1/2 (the original choice here) aren't an
        // aliased pair even under the old bug, so they couldn't have
        // detected it -- confirmed by independent verification during
        // review.
        let params_a = EmbedParameters {
            random_seed: 0,
            ..EmbedParameters::default()
        };
        let params_b = EmbedParameters {
            random_seed: 1,
            ..EmbedParameters::default()
        };
        let c1 = embed_distance_geometry_v2(&mol, &params_a).unwrap();
        let c2 = embed_distance_geometry_v2(&mol, &params_b).unwrap();
        let mut any_diff = false;
        for i in 0..mol.atom_count() {
            let p1 = c1.get(AtomIdx(i as u32));
            let p2 = c2.get(AtomIdx(i as u32));
            if p1 != p2 {
                any_diff = true;
            }
        }
        assert!(
            any_diff,
            "different seeds should not reproduce identical coords"
        );
    }

    #[test]
    fn track_failures_defaults_to_no_breakdown() {
        let mol = parse("CC").unwrap();
        let params = EmbedParameters::default(); // track_failures: false
        let (_, stats) = embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        assert!(stats.failure_counts.is_empty());
        assert_eq!(stats.attempts_used, 1);
    }

    #[test]
    fn track_failures_records_real_counts_on_forced_failure() {
        // O=C1CC[C@H]2CCC[C@H]12: both declared stereocenters are ring-fused with no
        // acyclic bridge substituent, so `repair_stereo` structurally cannot fix
        // either one -- empirically confirmed to fail all of seeds 0..8 with
        // max_attempts=1 (verified while writing this test, not assumed). Exercises
        // the genuine per-attempt failure path: every attempt is actually made and
        // recorded, unlike the old immediate-refusal behavior this replaced.
        let mol = parse("O=C1CC[C@H]2CCC[C@H]12").unwrap();
        let params = EmbedParameters {
            random_seed: 0,
            max_attempts: 3,
            enforce_chirality: true,
            track_failures: true,
            ..EmbedParameters::default()
        };
        let err = embed_distance_geometry_v2_detail(&mol, &params).unwrap_err();
        assert_eq!(err.0, EmbedFailureCause::StereoConstraintFailed);
        assert_eq!(
            err.1.attempts_used, 3,
            "every attempt must be genuinely made, not skipped"
        );
        assert_eq!(
            err.1
                .failure_counts
                .get(&EmbedFailureCause::StereoConstraintFailed),
            Some(&3),
            "each of the 3 attempts must be individually recorded"
        );
        assert!(err.1.stereo_repair_attempted);
    }

    #[test]
    fn enforce_chirality_repairs_and_succeeds_for_a_bridge_eligible_stereocenter() {
        // 2-butanol has one declared stereocenter with a plain acyclic substituent
        // (the ethyl group), so `repair_stereo` can always fix a wrong raw embedding
        // -- but `refine_coords`'s bounds-driven relaxation that runs after the
        // repair has no chirality awareness of its own and can occasionally flip a
        // correctly-repaired center back to `Violated` on a given attempt
        // (empirically confirmed for this exact molecule at random_seed=2,
        // max_attempts=1, while designing this: repair succeeds pre-refine, then
        // POST-REFINE is Violated again). That is exactly what the per-attempt
        // retry loop (`max_attempts`, distinct derived seed each try) exists to
        // absorb -- see the module doc's "`enforce_chirality`" section -- so this
        // asserts the retry loop delivers a satisfying geometry, not that any
        // single attempt does.
        let mol = parse("C[C@H](O)CC").unwrap();
        for seed in 0..20u64 {
            let params = EmbedParameters {
                random_seed: seed,
                max_attempts: 8,
                enforce_chirality: true,
                ..EmbedParameters::default()
            };
            let (coords, _stats) =
                embed_distance_geometry_v2_detail(&mol, &params).unwrap_or_else(|e| {
                    panic!("seed {seed} must succeed within 8 attempts, got {e:?}")
                });
            assert!(
                verify_stereo(&mol, &coords).is_fully_satisfied(),
                "seed {seed}: returned geometry must satisfy declared stereo"
            );
        }
    }

    #[test]
    fn enforce_chirality_false_is_unaffected_by_declared_stereo() {
        // Same molecule/seeds as the ring-fused failing case above, but with
        // enforce_chirality left at its default (false) -- must always succeed
        // (no stereo checking at all), proving the new logic is fully gated behind
        // the flag and doesn't leak into the default path.
        let mol = parse("O=C1CC[C@H]2CCC[C@H]12").unwrap();
        for seed in 0..8u64 {
            let params = EmbedParameters {
                random_seed: seed,
                max_attempts: 1,
                enforce_chirality: false,
                ..EmbedParameters::default()
            };
            let (_, stats) = embed_distance_geometry_v2_detail(&mol, &params)
                .unwrap_or_else(|e| panic!("seed {seed} must succeed, got {e:?}"));
            assert!(!stats.stereo_repair_attempted);
        }
    }

    #[test]
    fn enforce_chirality_noop_without_declared_stereo() {
        let mol = parse("CCCC").unwrap(); // butane, no stereo at all
        let params = EmbedParameters {
            enforce_chirality: true,
            ..EmbedParameters::default()
        };
        let (_, stats) = embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        assert!(!stats.stereo_repair_attempted);
    }

    #[test]
    fn atom_limit_exceeded_is_typed_not_silent() {
        // Build a long acyclic chain past DG_MAX_ATOMS to hit the typed limit path
        // without needing to actually construct 500+ Bonds by hand -- reuse the real
        // limit constant so this test tracks it if it ever changes.
        let smiles: String = "C".repeat(DG_MAX_ATOMS + 1);
        let mol = parse(&smiles).unwrap();
        let params = EmbedParameters::default();
        let err = embed_distance_geometry_v2(&mol, &params).unwrap_err();
        assert_eq!(err, EmbedFailureCause::AtomLimitExceeded);
    }

    #[test]
    fn use_random_coords_path_still_produces_valid_geometry() {
        let mol = parse("c1ccccc1").unwrap();
        let params = EmbedParameters {
            use_random_coords: true,
            ..EmbedParameters::default()
        };
        let (coords, stats) = embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        assert!(stats.used_random_coords);
        for i in 0..coords.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        }
    }

    #[test]
    fn wildcard_atom_is_invalid_topology() {
        let mol = parse("[*]C").unwrap();
        let params = EmbedParameters::default();
        let err = embed_distance_geometry_v2(&mol, &params).unwrap_err();
        assert_eq!(err, EmbedFailureCause::InvalidTopology);
    }

    #[test]
    fn empty_molecule_is_trivially_ok() {
        let mol = chematic_core::MoleculeBuilder::new().build();
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).unwrap();
        assert_eq!(coords.atom_count(), 0);
    }

    #[test]
    fn smoothing_invariants_hold_on_a_real_molecule() {
        // Direct check of the invariant helper against dg_fft's own bound matrix +
        // smoothing on a fused-ring molecule (naphthalene), matching the Coordinator's
        // requested "before vs after" check.
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let (lower0, upper0) = crate::dg_fft::build_bound_matrix(&mol);
        let mut lower = lower0.clone();
        let mut upper = upper0.clone();
        crate::dg_fft::smooth_bounds(&mut lower, &mut upper);
        assert!(smoothing_preserves_invariants(
            &lower0, &upper0, &lower, &upper
        ));
    }

    #[test]
    fn negative_eigenvalues_are_reported_not_silently_dropped() {
        // Not every molecule produces negative eigenvalues, but the field must always
        // be *present* and consistent (0 negative eigenvalues => magnitude 0.0).
        let mol = parse("c1ccc2ccccc2c1").unwrap(); // naphthalene: fused ring, good stress case
        let params = EmbedParameters::default();
        let (_, stats) = embed_distance_geometry_v2_detail(&mol, &params).unwrap();
        if stats.negative_eigenvalues_beyond_embedding_dim == 0 {
            assert_eq!(stats.max_negative_eigenvalue_magnitude, 0.0);
        } else {
            assert!(stats.max_negative_eigenvalue_magnitude > 0.0);
        }
    }

    #[test]
    fn bounds_conformance_reports_nonzero_residual_on_a_real_molecule() {
        // refine_coords is a heuristic, not a guaranteed-convergent solver -- confirm
        // bounds_conformance actually measures real (nonzero) residual violation on a
        // stress case, not just always reporting a trivial 0/0.
        let mol = parse("CN(C)CCOC(c1ccccc1)c1ccccc1").unwrap(); // diphenhydramine
        let params = EmbedParameters::default();
        let coords = embed_distance_geometry_v2(&mol, &params).unwrap();
        let bc = bounds_conformance(&mol, &coords);
        assert!(bc.n_pairs > 0);
        // Not asserting a specific violation count (that's this PR's honestly-reported
        // measured number, not a hardcoded expectation) -- just that the function
        // computes something real and internally consistent.
        assert!(bc.n_violations <= bc.n_pairs);
        if bc.n_violations == 0 {
            assert_eq!(bc.max_rel_violation, 0.0);
        } else {
            assert!(bc.max_rel_violation > 0.0);
        }
    }

    /// Stochastic gate (4-layer acceptance gate, layer 4): "uniqueness of derived
    /// per-attempt seed streams" -- if `max_attempts` allows multiple internal
    /// retries, each attempt must draw from a genuinely different seed, not
    /// accidentally repeat one (which would make a retry after a bad stochastic
    /// draw just redraw the same bad geometry). Checked directly against
    /// `derive_attempt_seed` rather than at the full-embedding level, since this is
    /// an internal-derivation invariant, not an observable-geometry one.
    #[test]
    fn derive_attempt_seed_is_distinct_across_attempts() {
        use std::collections::HashSet;
        for base in [0u64, 1, 42, 0xC0FF_EE42_D157_6E02, u64::MAX] {
            let seeds: HashSet<u64> = (0..32)
                .map(|attempt| derive_attempt_seed(base, attempt))
                .collect();
            assert_eq!(
                seeds.len(),
                32,
                "base seed {base:#x}: expected 32 distinct per-attempt seeds, got {}",
                seeds.len()
            );
        }
    }
}
