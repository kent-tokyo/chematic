//! Real (stochastic) distance-geometry conformer embedding — 3D Breakthrough Program,
//! Wave 1, Agent C. See `docs/rfcs/3d_breakthrough_master_plan.md` §1b/§4 and
//! `docs/rfcs/etkdg_3d_gap_rfc.md` for the diagnosis this module answers.
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
//! # Not wired into the default conformer path (stale note corrected 2026-08-11)
//!
//! `etkdg.rs`'s `generate_coords_etkdg` and `Mol.conformer_ensemble()`
//! (`generate_conformer_ensemble`/`generate_conformer_ensemble_with_config`) still
//! use the older `dg.rs` deterministic-midpoint embedder, not this module, and this
//! remains true today. What changed since this note was first written: `pipeline_v2.rs`
//! (Wave 2 → Wave 3 Coordinator Integration 1) now calls into this module directly
//! ([`embed_distance_geometry_v2_with_adjustments`]), and `pipeline_v2::embed_pipeline_v2`
//! itself **is** reachable from outside Rust-internal code — `Mol.embed_pipeline_v2()`
//! (chematic-py) and `embed_pipeline_v2_json` (chematic-wasm) both call it, an opt-in
//! surface distinct from (and not yet routed through) the default conformer path above.
//! `chematic-mcp` does not expose it. See `pipeline_v2.rs`'s own module doc for its
//! stage order and the `enforce_chirality`/`StereoPolicy` interaction.
//!
//! # `enforce_chirality` (added 2026-08-11, issue #291/#293; E/Z bound fix
//! 2026-08-11, issue #285)
//!
//! **Important, not obvious from the name**: a pairwise distance matrix is
//! reflection-invariant (a molecule and its mirror image have identical pairwise
//! distances — see `docs/rfcs/etkdg_3d_gap_rfc.md`'s Phase 3 "Correction" section), so
//! nothing in `build_bound_matrix`/`smooth_bounds`/the MDS embedding step above can
//! ever encode which **tetrahedral** (`@`/`@@`) enantiomer to prefer. **This does
//! NOT extend to declared E/Z**: cis and trans are not mirror images of each other,
//! they are two different scalar 1-4 separations for the same substituent pair, so a
//! distance bound *can* rule one out — see `apply_declared_ez_bounds` below, which
//! does exactly that. For tetrahedral centers the limitation is real and this
//! function doesn't touch it: `enforce_chirality` instead falls back to a repair-
//! after-embed strategy for those. Per attempt, after the raw MDS/random placement
//! (which operates on real coordinates, unlike the bound matrix, and so *can* carry a
//! tetrahedral chirality sign): apply [`stereo_constraints::repair_stereo`] to fix
//! whatever it can (bridge-eligible substituents only — ring-fused stereocenters
//! cannot be fixed this way, see that module's docs) BEFORE `refine_coords` runs, so
//! the bounds-driven relaxation absorbs whatever bond-length distortion the repair's
//! rigid-subtree translation introduced. Then, after `refine_coords`, independently
//! re-verify the *actual returned* geometry (never trust that refinement preserved
//! the repair). If that re-verify still finds a violation, `try_embed_once` runs one
//! more `repair_stereo` pass directly on the post-refinement geometry (no further
//! `refine_coords` call) as a safety net — measured directly (issue #291) for
//! ring-fused, multi-stereocenter molecules like testosterone/cholesterol,
//! `refine_coords`'s chirality-blind bound correction routinely undoes the
//! pre-refinement repair, and this second pass recovers it without reopening the
//! bounds-distortion problem the pre-refinement ordering exists to avoid, since
//! `repair_stereo`'s rigid-subtree translation preserves bond lengths exactly. Only
//! if *that* also fails to verify does the attempt fail with
//! [`EmbedFailureCause::StereoConstraintFailed`] and the existing per-attempt retry
//! loop tries a new seed. Molecules with no declared stereo are unaffected
//! (`repair_stereo`/`verify_stereo` are no-ops on them); `enforce_chirality: false`
//! (the default) is byte-identical to before this change — opt-in only, see
//! `ROADMAP.md`'s v0.14.0 S-tier item 1 for why the wider default-path decision
//! (issue #291's `Ignore`-policy population) was deliberately deferred, not folded
//! in here.
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
//! ## The `but2ene_Z` gap (issue #285), root-caused and fixed
//!
//! An earlier measurement (2026-08-11) found `but2ene_Z` (`C/C=C\C`, a plain acyclic
//! Z-alkene, no ring) failing 4/5 base seeds — NOT a retry-odds problem: the *raw*
//! (pre-repair) embedding was `Violated` for all of 10 tested derived seeds,
//! deterministically, contrasted against a larger/flexible molecule with the same
//! declared-E/Z shape (`cinnamic_acid_E`) whose raw sign genuinely split 5/10. The
//! mechanism was isolated by stage-by-stage instrumentation (bounds → sampling →
//! Gram matrix → eigendecomposition → MDS reconstruction → `refine_coords` →
//! verifier) and an earlier candidate explanation -- that `jacobi_eigendecompose`'s
//! eigenvector-sign outcome is a continuous function of the input Gram matrix that a
//! tightly-bounded molecule's samples never cross -- was **empirically refuted**: the
//! pre-refine scalar-triple-product sign genuinely varies seed to seed for
//! `but2ene_Z`, same as for the flexible molecule.
//!
//! The actual mechanism: `apply_vdw_bounds`'s generic non-bonded lower bound (sum of
//! Van der Waals radii; two carbons: 3.40 Å) was being applied to the declared-Z
//! alkene's own 1-4 substituent pair, whose analytically-correct cis separation
//! (using the exact `ideal_bond_length`/`ideal_bond_angle` model
//! `build_bond_angle_bounds` already uses, extended one bond further — see
//! `declared_1_4_distance`) is ≈2.88 Å for `but2ene_Z` — *below* that VDW floor. The
//! smoothed bound `[3.400, 4.186]` this produced was identical for every tested
//! molecule regardless of declared E/Z (VDW's contribution is generic, not
//! stereo-aware), and structurally excluded the correct cis geometry from ever being
//! sampled or reconstructed; `refine_coords`'s bounds-driven relaxation then pulled
//! every seed's differently-signed starting point toward the same (wrong, ~3.4-3.7 Å)
//! basin, which is what produced the *appearance* of deterministic sign-fixation. The
//! analytic trans separation (≈3.93 Å) sat comfortably inside that same generic
//! bound, which is why `but2ene_E` was unaffected (10/10 raw).
//!
//! Fixed by [`apply_declared_ez_bounds`]: for each declared E/Z double bond, compute
//! the analytic same-side/opposite-side 1-4 distance for its
//! [`stereo_constraints::build_stereo_constraints`]-normalized substituent pair, and
//! intersect a ±0.1 Å window around it into the bound matrix *before* the VDW loop
//! runs (not after — see that function's own doc for why the ordering matters: VDW's
//! non-bonded assumption doesn't hold for a pair whose separation is fixed by nearer,
//! more specific declared stereochemistry, the same exemption bonded/1-3 pairs
//! already get). `enforce_chirality: false` never calls this function. Deterministic
//! reflection injection (post-embed) was considered and rejected: reflection is a
//! distance-preserving isometry, so it cannot fix a distance-*magnitude*
//! infeasibility — it only flips the dihedral-sign classification `verify_double_bond`
//! checks, which could make a still-wrong-distance geometry read as `Satisfied`. Gram-
//! matrix perturbation was also rejected: the true cis geometry is provably outside
//! the (pre-fix) feasible region, so no amount of noise samples it more often than
//! chance — stereo assignment shouldn't depend on getting lucky.
//!
//! **Measured (2026-08-11, `max_attempts: 8`, the 29 stereo-bearing molecules of
//! `scripts/etkdg_vs_rdkit_gap.py::CORPUS`, `random_seed` swept over 0..5 — see
//! `ez_bounds_29_corpus_regression_but2ene_z_fixed_nothing_else_broken` for the exact
//! before/after sets, independently re-measured on unmodified `main` for the "before"
//! side). This measures the embedder alone (no UFF minimization, unlike issue #291's
//! `embed_pipeline_v2`-level 18/29 figure — the two numbers are not directly
//! comparable):** 26/29 succeed across the 5 base seeds tested (up from 25/29 before
//! this fix) — `but2ene_Z` now passes 5/5 (was 1/5), and every molecule that fully
//! passed before this fix still fully passes after (zero newly-broken molecules,
//! verified directly, not inferred). The only remaining recurring (5/5) failures:
//! - `testosterone`, `cholesterol` — expected: ring-fused declared **tetrahedral**
//!   stereocenters, the documented `NoBridgeEligibleSubstituent` case above,
//!   unaffected by this fix (declared-scoped to E/Z only — see `apply_declared_ez_
//!   bounds`'s doc; a tetrahedral chiral-volume-penalty approach is separately-scoped
//!   future work, not this PR).
//!
//! `cinnamic_acid_E` continues to show the same pre-existing, unrelated flexible-
//! molecule retry variance (4/5, same both before and after this fix) — expected, see
//! `ez_bounds_cinnamic_acid_e_retry_loop_still_resolves_flexible_variance`.

#[allow(unused_imports)]
use crate::fmath::DetMath;
use std::collections::HashMap;

use crate::clock::Instant;

use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule};

use crate::coords::{Coords3D, Point3};
use crate::dg_fft::{
    DG_MAX_ATOMS, apply_vdw_bounds, build_bond_angle_bounds, build_bound_matrix,
    center_coordinates, distance_to_gram_matrix, ideal_bond_angle, ideal_bond_length,
    jacobi_eigendecompose, refine_coords, smooth_bounds,
};
use crate::prng::Prng;
use crate::stereo_constraints::{build_stereo_constraints, repair_stereo, verify_stereo};

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
    /// Has no effect unless `enforce_chirality` is also `true`. When both are
    /// set and the molecule has declared stereo, every attempt embeds a
    /// temporary copy of `mol` with every implicit hydrogen materialized as
    /// a real, explicit atom (`chematic_chem::add_hydrogens`), then truncates
    /// the returned coordinates back to the caller's original atom count —
    /// the external coordinate contract (one entry per atom of the molecule
    /// passed in) is unchanged either way.
    ///
    /// Exists because `repair_tetrahedral_center`'s substituent-reflection
    /// repair has no coordinate to move for an *implicit* H — a declared
    /// tetrahedral center whose only non-ring substituent is an implicit H
    /// (e.g. a ring-fusion carbon in a steroid) is unrepairable as-is
    /// (`RepairRejectionReason::NoBridgeEligibleSubstituent`). A real,
    /// explicit H atom is terminal and non-ring, so it becomes a valid
    /// reflection candidate with no new geometry mechanism needed.
    ///
    /// Roughly doubles atom count for an all-implicit-H molecule, so
    /// `refine_coords`'s O(n²)-per-iteration cost rises accordingly (up to
    /// ~4x) — opt-in, not automatic, for exactly this reason. Molecules with
    /// several *simultaneous* ring-fused declared stereocenters (e.g.
    /// cholesterol) may also need a substantially higher `max_attempts` to
    /// reliably converge than a single easier stereocenter would (measured:
    /// low single-digit attempts typically suffice for one such center, but
    /// cholesterol's three needed ~100 attempts for a reliable draw) — this
    /// field does not change `max_attempts`'s default; the caller absorbs
    /// that cost explicitly, same as `enforce_chirality`'s existing
    /// retry-on-violation behavior already works today. `false` (default)
    /// is byte-identical to this field's behavior before it existed.
    /// **Consumed.**
    pub materialize_implicit_h_for_chirality: bool,
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
            materialize_implicit_h_for_chirality: false,
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
    /// A structurally inconsistent bound was produced before smoothing even ran (e.g.
    /// `lower > upper`) — a bug in bounds construction, not a smoothing or sampling
    /// issue. Two distinct sources: an actually-bonded pair (a `dg_fft::
    /// build_bond_angle_bounds` defect), or -- `enforce_chirality` only -- a declared-
    /// E/Z 1-4 pair whose stereo-derived window has no overlap with its existing
    /// bond/angle bound (see `apply_declared_ez_bounds`'s doc; not expected for a
    /// chemically sane molecule, but reported as this variant, not silently ignored,
    /// if it happens).
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

    let original_n = mol.atom_count();

    // Opt-in: embed a temporary H-expanded copy so `repair_tetrahedral_center`
    // has a real, movable substituent for ring-fused declared stereocenters
    // whose only non-ring "substituent" is an implicit H (see
    // `EmbedParameters::materialize_implicit_h_for_chirality`'s own doc
    // comment). `add_hydrogens` keeps every heavy atom at its original index
    // (new H atoms are strictly appended), so every downstream index
    // (`adjustments`, bound-matrix sizing) stays valid unchanged against the
    // larger molecule -- only the final coordinates need truncating back
    // down before returning to the caller.
    let expanded_mol;
    let mol: &Molecule = if params.materialize_implicit_h_for_chirality
        && params.enforce_chirality
        && mol_has_declared_stereo(mol)
    {
        expanded_mol = chematic_chem::add_hydrogens(mol);
        &expanded_mol
    } else {
        mol
    };

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
            Ok(coords) => {
                let coords = if n != original_n {
                    truncate_coords(&coords, original_n)
                } else {
                    coords
                };
                return Ok((coords, stats));
            }
            Err(cause) => {
                last_cause = cause;
                record_failure(&mut stats, params, cause);
            }
        }
    }

    Err((EmbedWithAdjustmentsFailure::Embed(last_cause), stats))
}

/// Drop every atom from index `n` onward -- used to map a temporarily
/// H-expanded embed (see `materialize_implicit_h_for_chirality`) back onto
/// the caller's original atom count. Heavy atoms keep their original index
/// under `add_hydrogens`, so a plain prefix copy is exact, not an
/// approximation. `pub(crate)`: also reused by `pipeline_v2.rs`'s own
/// H-expanded-geometry path (issue #291), which needs to truncate at its own
/// stage 11/12 boundary rather than immediately after embed.
pub(crate) fn truncate_coords(coords: &Coords3D, n: usize) -> Coords3D {
    let mut out = Coords3D::new_zeroed(n);
    for i in 0..n {
        out.set(AtomIdx(i as u32), coords.get(AtomIdx(i as u32)));
    }
    out
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

/// `pub(crate)`: also reused by `pipeline_v2.rs`'s own H-expanded-geometry
/// path (issue #291) to gate expansion the same way this module already
/// gates `materialize_implicit_h_for_chirality`'s own effect.
pub(crate) fn mol_has_declared_stereo(mol: &Molecule) -> bool {
    if (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None) {
        return true;
    }
    mol.bonds()
        .any(|(_, bond)| matches!(bond.order, BondOrder::Up | BondOrder::Down))
}

/// Tolerance (Å) around the declared-E/Z analytic 1-4 distance, matching
/// `dg_fft::build_bond_angle_bounds`'s own angle-derived (1-3) bound tolerance --
/// this bound is one bond further out on the same generic bond-length/bond-angle
/// model, not a different precision claim.
const EZ_BOUND_TOLERANCE: f64 = 0.1;

/// Analytic 1-4 distance between `sub1` (on `end1`'s side) and `sub2` (on `end2`'s
/// side) of the planar `sub1-end1=end2-sub2` fragment, given whether they're declared
/// same-side or opposite-side. Uses the exact same `ideal_bond_length`/
/// `ideal_bond_angle` model `dg_fft::build_bond_angle_bounds` already uses for its
/// 1-2/1-3 bounds -- this is that same law-of-cosines-chain construction extended one
/// bond further, not a new geometric model.
///
/// Places `end1` at the origin and `end2` along +x (so the `end1`-`end2` bond has
/// direction 0 rad); `sub1` sits at angle `ideal_bond_angle(mol, end1)`
/// counterclockwise from the `end1`->`end2` ray (that angle's definition, measured at
/// `end1`); `sub2` sits at `ideal_bond_angle(mol, end2)` from the `end2`->`end1` ray,
/// on the same side as `sub1` (+y) if `same_side`, the opposite side (-y) otherwise.
fn declared_1_4_distance(
    mol: &Molecule,
    end1: AtomIdx,
    end2: AtomIdx,
    sub1: AtomIdx,
    sub2: AtomIdx,
    same_side: bool,
) -> f64 {
    let d_sub1 = ideal_bond_length(mol, sub1, end1);
    let d_ends = ideal_bond_length(mol, end1, end2);
    let d_sub2 = ideal_bond_length(mol, end2, sub2);
    let angle1 = ideal_bond_angle(mol, end1);
    let angle2 = ideal_bond_angle(mol, end2);

    let sub1_pos = (d_sub1 * angle1.dcos(), d_sub1 * angle1.dsin());
    let sign = if same_side { 1.0 } else { -1.0 };
    let sub2_pos = (
        d_ends - d_sub2 * angle2.dcos(),
        sign * d_sub2 * angle2.dsin(),
    );

    ((sub1_pos.0 - sub2_pos.0).powi(2) + (sub1_pos.1 - sub2_pos.1).powi(2)).sqrt()
}

/// Add declared-E/Z-derived 1-4 distance bounds to the (bond/angle-only, pre-VDW)
/// bound matrix, one pair per declared E/Z double bond, intersected with (never
/// overwriting) whatever bound already exists for that pair. Only called from
/// `try_embed_once` when `enforce_chirality` is set -- `enforce_chirality: false`
/// never calls this, keeping that path byte-identical to before this function
/// existed.
///
/// # Why this is sound for E/Z but not tetrahedral chirality
///
/// A pairwise distance matrix is reflection-invariant -- a molecule and its mirror
/// image have identical pairwise distances -- so nothing here can ever encode which
/// tetrahedral (`@`/`@@`) enantiomer to prefer (see the module doc's opening
/// paragraph; that limitation is real and unaffected by this function). Declared E/Z
/// is different in kind, not just degree: cis and trans are not mirror images of each
/// other, they are two genuinely different scalar 1-4 separations for the same atom
/// pair (empirically confirmed for `C/C=C\C`: analytic cis 1-4 ≈ 2.88 Å vs. analytic
/// trans 1-4 ≈ 3.93 Å, using this exact function -- see the PR body). A distance bound
/// *can* rule one of them out, so this function does that directly at the bound-
/// construction stage, rather than repairing/retrying after the fact.
///
/// # Why this must run before the VDW loop, not after
///
/// `apply_vdw_bounds` unconditionally raises a non-bonded pair's lower bound to the
/// sum of Van der Waals radii (e.g. two carbons: 3.40 Å) -- this is exactly what makes
/// a declared-cis small alkene's raw embedding structurally unable to reach the
/// correct ~2.88 Å geometry today (see the PR body's root-cause diagnosis). Calling
/// this function on the bond/angle-only bounds (`build_bond_angle_bounds`'s output,
/// before `apply_vdw_bounds` runs) lets `apply_vdw_bounds`'s own existing guard
/// (`vdw_sum <= upper[i][j]`) correctly skip this pair once its declared-E/Z upper
/// bound is already tighter than the generic VDW floor -- the same exemption pattern
/// bonded and 1-3 pairs already get from that loop, just reached one pair further out.
/// Intersecting *after* `apply_vdw_bounds` instead would make the VDW floor win for
/// every declared-cis case (their true separation is routinely below it), producing an
/// empty intersection and `BoundsConstructionFailed` on every one -- technically
/// fail-closed, but never actually delivering the corrected geometry this exists for.
///
/// # Normalization
///
/// Reuses `build_stereo_constraints`'s already-normalized `DoubleBondConstraint`
/// (`end1`/`end2`/`sub1`/`sub2`/`same_side`) rather than reading `/`/`\` markers
/// directly -- the same central-double-bond + stereo-defining-substituent-pair +
/// declared-relation shape `verify_double_bond`/`repair_double_bond` already use, so
/// this can never disagree with what `enforce_chirality`'s own verify/repair step
/// considers declared.
fn apply_declared_ez_bounds(
    mol: &Molecule,
    lower: &mut [Vec<f64>],
    upper: &mut [Vec<f64>],
) -> Result<(), EmbedFailureCause> {
    let constraints = build_stereo_constraints(mol);
    for c in &constraints.double_bond {
        let (i, j) = (c.sub1.0 as usize, c.sub2.0 as usize);
        if i == j {
            continue; // degenerate (shared substituent); not expected, nothing to bound
        }
        let dist = declared_1_4_distance(mol, c.end1, c.end2, c.sub1, c.sub2, c.same_side);
        let stereo_lo = (dist - EZ_BOUND_TOLERANCE).max(0.0);
        let stereo_hi = dist + EZ_BOUND_TOLERANCE;

        let new_lo = lower[i][j].max(stereo_lo);
        let new_hi = upper[i][j].min(stereo_hi);
        if new_lo > new_hi + INVARIANT_EPS {
            return Err(EmbedFailureCause::BoundsConstructionFailed);
        }
        lower[i][j] = new_lo;
        lower[j][i] = new_lo;
        upper[i][j] = new_hi;
        upper[j][i] = new_hi;
    }
    Ok(())
}

/// Derive a per-attempt seed from the caller's base seed. Deterministic: the same
/// `(base, attempt)` pair always produces the same derived seed, so a fixed
/// `random_seed` reproduces the exact same sequence of attempts (and thus the exact
/// same final result) every time.
///
/// `pub(crate)`, not private: reused as-is by `ensemble_v2`'s outer
/// multi-conformer loop (A2) rather than inventing a second derivation scheme
/// for the same "derive N distinct, deterministic seeds from one base seed"
/// need. Still not part of this crate's public API.
pub(crate) fn derive_attempt_seed(base: u64, attempt: usize) -> u64 {
    const GOLDEN: u64 = 0x9E37_79B9_7F4A_7C15;
    const OFFSET: u64 = 0xD1B5_4A32_D192_ED03;
    base ^ (attempt as u64).wrapping_mul(GOLDEN).wrapping_add(OFFSET)
}

type BoundMatrix = Vec<Vec<f64>>;

fn build_attempt_bounds(
    mol: &Molecule,
    enforce_chirality: bool,
) -> Result<(BoundMatrix, BoundMatrix), EmbedFailureCause> {
    if !enforce_chirality {
        return Ok(build_bound_matrix(mol));
    }
    let (mut lower, mut upper) = build_bond_angle_bounds(mol);
    apply_declared_ez_bounds(mol, &mut lower, &mut upper)?;
    apply_vdw_bounds(mol, &mut lower, &mut upper);
    Ok((lower, upper))
}

fn apply_bound_adjustments(
    lower: &mut [Vec<f64>],
    upper: &mut [Vec<f64>],
    adjustments: &[DistanceBoundAdjustment],
    stats: &mut EmbedStats,
) {
    for adjustment in adjustments {
        let i = adjustment.atom1.0 as usize;
        let j = adjustment.atom2.0 as usize;
        lower[i][j] = adjustment.lower;
        lower[j][i] = adjustment.lower;
        upper[i][j] = adjustment.upper;
        upper[j][i] = adjustment.upper;
    }
    if !adjustments.is_empty() {
        stats.adjustments_applied = adjustments.len();
    }
}

fn validate_bond_bounds(
    mol: &Molecule,
    lower: &[Vec<f64>],
    upper: &[Vec<f64>],
) -> Result<(), EmbedFailureCause> {
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        if lower[i][j] > upper[i][j] + INVARIANT_EPS {
            return Err(EmbedFailureCause::BoundsConstructionFailed);
        }
    }
    Ok(())
}

fn smooth_attempt_bounds(
    lower0: &[Vec<f64>],
    upper0: &[Vec<f64>],
    stats: &mut EmbedStats,
) -> Result<(BoundMatrix, BoundMatrix), EmbedFailureCause> {
    let mut lower = lower0.to_vec();
    let mut upper = upper0.to_vec();
    smooth_bounds(&mut lower, &mut upper);
    let invariants_ok = smoothing_preserves_invariants(lower0, upper0, &lower, &upper);
    stats.last_smoothing_invariants_ok = invariants_ok;
    if !invariants_ok {
        return Err(EmbedFailureCause::BoundsSmoothingFailed);
    }
    Ok((lower, upper))
}

fn fallback_distance(upper: &[Vec<f64>]) -> f64 {
    let max_finite_upper = upper
        .iter()
        .flat_map(|row| row.iter())
        .filter(|value| value.is_finite())
        .cloned()
        .fold(0.0_f64, f64::max);
    (max_finite_upper * 4.0).max(10.0)
}

fn initial_attempt_coords(
    atom_count: usize,
    lower: &[Vec<f64>],
    upper: &[Vec<f64>],
    params: &EmbedParameters,
    seed: u64,
    stats: &mut EmbedStats,
) -> Result<Coords3D, EmbedFailureCause> {
    let fallback = fallback_distance(upper);
    let mut coords = if params.use_random_coords {
        stats.used_random_coords = true;
        random_initial_coords(atom_count, fallback, seed)
    } else {
        stats.used_random_coords = false;
        let distances = sample_distance_matrix(atom_count, lower, upper, fallback, seed);
        mds_embed(&distances, stats)?
    };
    center_coordinates(&mut coords);
    Ok(coords)
}

fn repair_declared_stereo(mol: &Molecule, coords: Coords3D, stats: &mut EmbedStats) -> Coords3D {
    stats.stereo_repair_attempted = true;
    match repair_stereo(mol, &coords) {
        Ok(outcome) => outcome.coords,
        Err(failure) => failure.partial_coords,
    }
}

fn finalize_attempt_stereo(
    mol: &Molecule,
    mut coords: Coords3D,
    enforce_chirality: bool,
) -> Result<Coords3D, EmbedFailureCause> {
    if !enforce_chirality
        || !mol_has_declared_stereo(mol)
        || verify_stereo(mol, &coords).is_fully_satisfied()
    {
        return Ok(coords);
    }
    coords = match repair_stereo(mol, &coords) {
        Ok(outcome) => outcome.coords,
        Err(failure) => failure.partial_coords,
    };
    validate_final_coords(mol, &coords)?;
    if verify_stereo(mol, &coords).is_fully_satisfied() {
        Ok(coords)
    } else {
        Err(EmbedFailureCause::StereoConstraintFailed)
    }
}

fn try_embed_once(
    mol: &Molecule,
    params: &EmbedParameters,
    seed: u64,
    stats: &mut EmbedStats,
    adjustments: &[DistanceBoundAdjustment],
) -> Result<Coords3D, EmbedFailureCause> {
    let (mut lower0, mut upper0) = build_attempt_bounds(mol, params.enforce_chirality)?;
    apply_bound_adjustments(&mut lower0, &mut upper0, adjustments, stats);
    validate_bond_bounds(mol, &lower0, &upper0)?;
    let (lower, upper) = smooth_attempt_bounds(&lower0, &upper0, stats)?;
    let mut coords = initial_attempt_coords(mol.atom_count(), &lower, &upper, params, seed, stats)?;

    // Repair BEFORE refinement, not after: `repair_stereo`'s rigid-subtree
    // translation preserves bond lengths *within* the moved substituent exactly,
    // but can and does push its distances to the *rest* of the molecule outside
    // their smoothed bounds (measured empirically while designing this: e.g.
    // L-alanine's bounds_conformance violations went 1->4 when repair was applied
    // AFTER refine_coords). Applying it here lets the bounds-driven relaxation
    // below absorb that distortion the same way it absorbs everything else -- see
    // the module doc's "`enforce_chirality`" section.
    if params.enforce_chirality && mol_has_declared_stereo(mol) {
        coords = repair_declared_stereo(mol, coords, stats);
    }

    refine_coords(&mut coords, &lower, &upper, REFINE_ITERS);

    validate_final_coords(mol, &coords)?;

    finalize_attempt_stereo(mol, coords, params.enforce_chirality)
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
    let box_size = radius_scale.max(3.0) * (n as f64).dcbrt();
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
#[path = "distance_geometry_v2_tests.rs"]
mod embedded_tests;
