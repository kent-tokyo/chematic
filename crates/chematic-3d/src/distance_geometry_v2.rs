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

    let sub1_pos = (d_sub1 * angle1.cos(), d_sub1 * angle1.sin());
    let sign = if same_side { 1.0 } else { -1.0 };
    let sub2_pos = (d_ends - d_sub2 * angle2.cos(), sign * d_sub2 * angle2.sin());

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

fn try_embed_once(
    mol: &Molecule,
    params: &EmbedParameters,
    seed: u64,
    stats: &mut EmbedStats,
    adjustments: &[DistanceBoundAdjustment],
) -> Result<Coords3D, EmbedFailureCause> {
    let n = mol.atom_count();

    // `enforce_chirality`-only: insert declared-E/Z 1-4 bounds between the bond/angle
    // bounds and the VDW floor (see `apply_declared_ez_bounds`'s doc for why that
    // ordering, not before-or-after as one call, matters). `enforce_chirality: false`
    // takes the untouched `build_bound_matrix` path, byte-identical to before this
    // function existed.
    let (mut lower0, mut upper0) = if params.enforce_chirality {
        let (mut lower0, mut upper0) = build_bond_angle_bounds(mol);
        apply_declared_ez_bounds(mol, &mut lower0, &mut upper0)?;
        apply_vdw_bounds(mol, &mut lower0, &mut upper0);
        (lower0, upper0)
    } else {
        build_bound_matrix(mol)
    };

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
        // Safety net for centers `refine_coords` un-repairs: for ring-fused
        // multi-stereocenter molecules (e.g. testosterone, cholesterol) the
        // pre-refinement repair above is routinely undone by refine_coords's
        // chirality-blind bound correction -- measured directly, every seed
        // tried failed here without this net. Repairing again on the
        // POST-refine geometry, with no further refine_coords call, resolves
        // it: `repair_stereo`'s rigid-subtree translation preserves bond
        // lengths exactly, so re-running `validate_final_coords` below (not
        // full `refine_coords`) is enough to catch any resulting distortion.
        // Only reached when the pre-refinement repair didn't already survive --
        // an attempt that was already `is_fully_satisfied()` here never enters this
        // branch, so its own returned coordinates are byte-identical to before this
        // change. That does NOT mean overall behavior is unchanged: an attempt that
        // used to fail this check now often succeeds instead, so which seed/attempt
        // the outer retry loop settles on, how many attempts it consumes, and
        // whether a molecule succeeds at all can all change -- by design, that's
        // the fix (see `cinnamic_acid_E`/`chembl_tier_b_0168`'s tests below, both
        // measurably 4/5 -> 5/5).
        coords = match repair_stereo(mol, &coords) {
            Ok(outcome) => outcome.coords,
            Err(failure) => failure.partial_coords,
        };
        validate_final_coords(mol, &coords)?;
        if !verify_stereo(mol, &coords).is_fully_satisfied() {
            return Err(EmbedFailureCause::StereoConstraintFailed);
        }
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

    /// FORMERLY a known limitation (`dg_fft::build_bound_matrix`'s angle-constraint
    /// loop treated every pair of a center atom's neighbors as a 1-3
    /// (through-center) relationship and unconditionally tightened their bound with
    /// the *generic* ideal angle (~109.5°/120°) -- but in a 3-membered ring, that
    /// "1-3" pair is *also* a direct 1-2 bonded pair, so the generic-angle bound
    /// overwrote the correct, much shorter bond-length bound with a contradictory
    /// one). Fixed by skipping the angle-derived bound for any neighbor pair that
    /// is itself directly bonded (`dg_fft::build_bond_angle_bounds`): a 3-membered
    /// ring's three bond-length constraints already fully determine its shape, so
    /// nothing is lost. See `cyclopropane_ring_closing_bond_uses_bond_length_bound_not_angle`
    /// for the exact numbers this replaced.
    #[test]
    fn three_membered_rings_embed_successfully() {
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
            let coords = embed_distance_geometry_v2(&mol, &params)
                .unwrap_or_else(|e| panic!("{smiles}: expected a successful embed, got {e:?}"));
            let worst = worst_bond(&mol, &coords);
            assert!(worst < 2.0, "{smiles}: worst bond length {worst}");
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
    fn cyclopropane_ring_closing_bond_uses_bond_length_bound_not_angle() {
        // The ring-closing C-C pair must keep the tight bond-length bound (upper
        // ≈ 1.59 Å) rather than being overwritten by the generic ~109.5°-angle-
        // derived bound (which previously forced lower ≈ 2.414 Å > upper, the
        // exact contradiction `three_membered_rings_embed_successfully` used to
        // fail closed on).
        let mol = parse("C1CC1").unwrap();
        let (lower, upper) = crate::dg_fft::build_bound_matrix(&mol);
        let mut found = false;
        for (_, bond) in mol.bonds() {
            let i = bond.atom1.0 as usize;
            let j = bond.atom2.0 as usize;
            assert!(
                lower[i][j] <= upper[i][j],
                "bonded pair ({i},{j}): lower={} > upper={}",
                lower[i][j],
                upper[i][j]
            );
            assert!((upper[i][j] - 1.590).abs() < 0.01, "upper={}", upper[i][j]);
            found = true;
        }
        assert!(found, "cyclopropane has no bonds?");
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

    // -------------------------------------------------------------------------
    // `apply_declared_ez_bounds` (issue #285 but2ene_Z root-cause fix)
    // -------------------------------------------------------------------------

    /// Primary gate: the RAW (`max_attempts: 1`, no retry) embedding itself must land
    /// in the declared-E/Z-compatible region, for every SMILES notation of the same
    /// declared stereochemistry, across many seeds. Before this fix, but2ene_Z's raw
    /// embedding was deterministically wrong-signed for all of 10 tested seeds (see
    /// the PR body's root-cause diagnosis) -- `max_attempts` retry could not help,
    /// because the bound matrix itself excluded the correct geometry. A test that only
    /// checked "does `embed_distance_geometry_v2_detail` eventually succeed" (its
    /// default `max_attempts: 8`) would not catch a regression back to that state, so
    /// this asserts every individual raw attempt is already correct.
    #[test]
    fn ez_bounds_but2ene_all_notation_variants_raw_embed_all_seeds_satisfied() {
        let cases: &[(&str, &str)] = &[
            ("Z, C/C=C\\C", r"C/C=C\C"),
            ("Z, C\\C=C/C", r"C\C=C/C"),
            ("E, C/C=C/C", "C/C=C/C"),
            ("E, C\\C=C\\C", r"C\C=C\C"),
        ];
        for (label, smiles) in cases.iter().copied() {
            let mol = parse(smiles).unwrap();
            for seed in 0..10u64 {
                let params = EmbedParameters {
                    random_seed: seed,
                    max_attempts: 1,
                    enforce_chirality: true,
                    ..EmbedParameters::default()
                };
                let (coords, _stats) = embed_distance_geometry_v2_detail(&mol, &params)
                    .unwrap_or_else(|e| {
                        panic!("but2ene ({label}) seed {seed}: raw embed must succeed, got {e:?}")
                    });
                assert!(
                    verify_stereo(&mol, &coords).is_fully_satisfied(),
                    "but2ene ({label}) seed {seed}: raw embed must already satisfy declared E/Z"
                );
            }
        }
    }

    /// Same raw-embed-all-seeds gate, generalized across the corpus's other small/
    /// rigid declared-E/Z alkenes (`scripts/etkdg_vs_rdkit_gap.py::CORPUS`'s
    /// `alkene_ez` group minus `cinnamic_acid_E`, see the next test) -- confirms the
    /// fix isn't a but2ene-specific coincidence.
    #[test]
    fn ez_bounds_rigid_alkene_corpus_raw_embed_all_seeds_satisfied() {
        let cases: &[(&str, &str)] = &[
            ("but2ene_E", "C/C=C/C"),
            ("but2ene_Z", r"C/C=C\C"),
            ("chloropropene_E", "C(/C=C/C)Cl"),
            ("chloropropene_Z", r"C(/C=C\C)Cl"),
            ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1"),
            ("pent2ene_E", "CC/C=C/C"),
            ("pent2ene_Z", r"CC/C=C\C"),
        ];
        for (name, smiles) in cases.iter().copied() {
            let mol = parse(smiles).unwrap();
            for seed in 0..10u64 {
                let params = EmbedParameters {
                    random_seed: seed,
                    max_attempts: 1,
                    enforce_chirality: true,
                    ..EmbedParameters::default()
                };
                let (coords, _stats) = embed_distance_geometry_v2_detail(&mol, &params)
                    .unwrap_or_else(|e| {
                        panic!("{name} seed {seed}: raw embed must succeed, got {e:?}")
                    });
                assert!(
                    verify_stereo(&mol, &coords).is_fully_satisfied(),
                    "{name} seed {seed}: raw embed must already satisfy declared E/Z"
                );
            }
        }
    }

    /// `cinnamic_acid_E` is the corpus's one alkene with genuine, pre-existing
    /// per-seed conformational variability (measured before this fix, see the module
    /// doc's original but2ene_Z diagnosis and the PR body: 3/10 raw sign-flip rate,
    /// larger/flexible molecule, NOT the small-rigid-bound-infeasibility mechanism
    /// this fix addresses) -- raw embedding is not expected to be all-seed-green for
    /// it, and this fix does not change that (its declared configuration is E, whose
    /// analytic 1-4 distance already sat inside the pre-existing VDW-derived bound).
    /// What must still hold: the ordinary retry loop (`max_attempts: 8`, the default)
    /// resolves it at the SAME rate as before this fix (independently measured on
    /// unmodified `main`, seeds 0..5: `[true, true, true, false, true]`, 4/5) -- a
    /// control against this PR accidentally making the flexible-molecule case worse
    /// while fixing the rigid one. Asserting the exact pre-fix count, not just
    /// "mostly succeeds", so a regression down to e.g. 3/5 or 2/5 is still caught.
    #[test]
    fn ez_bounds_cinnamic_acid_e_retry_loop_still_resolves_flexible_variance() {
        let mol = parse("OC(=O)/C=C/c1ccccc1").unwrap();
        let mut passes = 0usize;
        for seed in 0..5u64 {
            let params = EmbedParameters {
                random_seed: seed,
                max_attempts: 8,
                enforce_chirality: true,
                ..EmbedParameters::default()
            };
            if embed_distance_geometry_v2_detail(&mol, &params).is_ok() {
                passes += 1;
            }
        }
        // Was 4/5 before issue #291's post-refinement repair safety net
        // (`try_embed_once`'s second `repair_stereo` pass when the pre-refinement
        // repair doesn't survive `refine_coords`) -- that net recovers the one
        // seed that used to fail here too, independent of ring-fused centers.
        assert_eq!(
            passes, 5,
            "cinnamic_acid_E: expected 5/5 seeds to resolve within 8 attempts \
             (post-refinement repair safety net, issue #291), got {passes}/5"
        );
    }

    /// Issue #285's two named fixtures (`chembl_tier_b_0126`/`chembl_tier_b_0168`,
    /// each with one declared E double bond alongside two declared tetrahedral
    /// centers, one ring-fused): confirms this fix generalizes to a realistic
    /// drug-like molecule, not just an isolated 4-atom toy alkene, and doesn't break
    /// the standard retry loop's ability to resolve the (unrelated, unaffected by
    /// this PR) tetrahedral-center variance those molecules still have. Exact counts
    /// independently measured on unmodified `main`, seeds 0..5, `max_attempts: 8`:
    /// `chembl_tier_b_0126` 5/5, `chembl_tier_b_0168` 4/5 (`[true, true, true, true,
    /// false]`) -- asserted exactly, not just "mostly succeeds", so a regression is
    /// still caught. Re-measured after issue #291's post-refinement repair safety
    /// net (`try_embed_once`'s second `repair_stereo` pass): `chembl_tier_b_0168`'s
    /// one previously-failing seed now recovers too, so both are 5/5.
    #[test]
    fn ez_bounds_chembl_tier_b_0126_0168_retry_loop_succeeds() {
        let cases: &[(&str, &str, usize)] = &[
            (
                "chembl_tier_b_0126",
                "CC(=O)/C=C/CC1C(=O)N2[C@@H](C(=O)O)C(C)(C)S(=O)(=O)[C@@H]12",
                5,
            ),
            (
                "chembl_tier_b_0168",
                "CC(=O)/C=C/CC1C(=O)N2[C@@H](C(=O)O)C(C)(C)S(=O)(=O)[C@H]12",
                5,
            ),
        ];
        for (name, smiles, expected_passes) in cases.iter().copied() {
            let mol = parse(smiles).unwrap();
            let mut passes = 0usize;
            for seed in 0..5u64 {
                let params = EmbedParameters {
                    random_seed: seed,
                    max_attempts: 8,
                    enforce_chirality: true,
                    ..EmbedParameters::default()
                };
                if embed_distance_geometry_v2_detail(&mol, &params).is_ok() {
                    passes += 1;
                }
            }
            assert_eq!(
                passes, expected_passes,
                "{name}: expected {expected_passes}/5 seeds to resolve within 8 attempts \
                 (unchanged pre-existing tetrahedral-center variance), got {passes}/5"
            );
        }
    }

    /// Full re-measure of the 29 stereo-bearing molecules from
    /// `scripts/etkdg_vs_rdkit_gap.py::CORPUS` (same protocol as the module doc's
    /// original `enforce_chirality` measurement: 5 base seeds, `max_attempts: 8`),
    /// checked against the exact before-this-fix pass/fail set (independently
    /// re-measured on unmodified `main` while designing this fix, see the PR body) --
    /// not a re-derived guess. Two things must both hold:
    /// - `but2ene_Z` moves from failing (1/5 base seeds before) to fully passing
    ///   (5/5) -- the fix actually fires on the corpus, not just the isolated test
    ///   fixtures above.
    /// - Every molecule that fully passed (5/5) *before* this fix still fully passes
    ///   *after* -- zero newly-broken molecules. `testosterone`/`cholesterol` (ring-
    ///   fused tetrahedral centers) and `cinnamic_acid_E` (flexible-molecule variance)
    ///   are excluded from that "before" set on purpose: this PR is declared-scoped to
    ///   E/Z bounds only (see the module doc), so their pre-existing, unrelated
    ///   failure/partial-pass modes are expected to be unaffected, not fixed here.
    #[test]
    fn ez_bounds_29_corpus_regression_but2ene_z_fixed_nothing_else_broken() {
        const CORPUS_29: &[(&str, &str)] = &[
            ("l_alanine", "N[C@@H](C)C(=O)O"),
            ("d_alanine", "N[C@H](C)C(=O)O"),
            ("l_serine", "N[C@@H](CO)C(=O)O"),
            ("l_threonine", "C[C@H](O)[C@@H](N)C(=O)O"),
            ("2_butanol_R", "C[C@H](O)CC"),
            ("2_butanol_S", "C[C@@H](O)CC"),
            ("2_chlorobutane_R", "C[C@H](Cl)CC"),
            ("ibuprofen_S", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"),
            ("naproxen_S", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"),
            ("menthol", "C[C@@H]1CC[C@@H](C(C)C)C[C@H]1O"),
            ("chfclbr_R", "[C@H](F)(Cl)Br"),
            ("chfclbr_S", "[C@@H](F)(Cl)Br"),
            ("quaternary_1_R", "[C@](F)(Cl)(Br)I"),
            ("quaternary_1_S", "[C@@](F)(Cl)(Br)I"),
            ("quaternary_2_R", "[C@](C)(N)(O)F"),
            ("quaternary_2_S", "[C@@](C)(N)(O)F"),
            ("but2ene_E", "C/C=C/C"),
            ("but2ene_Z", r"C/C=C\C"),
            ("chloropropene_E", "C(/C=C/C)Cl"),
            ("chloropropene_Z", r"C(/C=C\C)Cl"),
            ("cinnamic_acid_E", "OC(=O)/C=C/c1ccccc1"),
            ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1"),
            ("pent2ene_E", "CC/C=C/C"),
            ("pent2ene_Z", r"CC/C=C\C"),
            (
                "penicillin_core",
                "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O",
            ),
            (
                "testosterone",
                "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
            ),
            (
                "cholesterol",
                "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
            ),
            (
                "atorvastatin_fragment",
                "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
            ),
            ("gly_ala_gly", "NCC(=O)N[C@@H](C)C(=O)NCC(=O)O"),
        ];

        // Independently re-measured on unmodified `main` (before this fix) with the
        // identical protocol below -- not a re-derived guess. Matches the module
        // doc's documented 25-27/29 range (this run: 25/29).
        const FULLY_PASSED_BEFORE_FIX: &[&str] = &[
            "2_butanol_R",
            "2_butanol_S",
            "2_chlorobutane_R",
            "atorvastatin_fragment",
            "but2ene_E",
            "chfclbr_R",
            "chfclbr_S",
            "chloropropene_E",
            "chloropropene_Z",
            "cinnamic_acid_Z",
            "d_alanine",
            "gly_ala_gly",
            "ibuprofen_S",
            "l_alanine",
            "l_serine",
            "l_threonine",
            "menthol",
            "naproxen_S",
            "penicillin_core",
            "pent2ene_E",
            "pent2ene_Z",
            "quaternary_1_R",
            "quaternary_1_S",
            "quaternary_2_R",
            "quaternary_2_S",
        ];

        let mut fully_passed_after: Vec<&str> = Vec::new();
        let mut but2ene_z_seeds_passed = 0usize;
        for (name, smiles) in CORPUS_29.iter().copied() {
            let mol = parse(smiles).unwrap();
            let mut passes = 0usize;
            for base_seed in 0..5u64 {
                let params = EmbedParameters {
                    random_seed: base_seed,
                    max_attempts: 8,
                    enforce_chirality: true,
                    ..EmbedParameters::default()
                };
                if embed_distance_geometry_v2_detail(&mol, &params).is_ok() {
                    passes += 1;
                }
            }
            if name == "but2ene_Z" {
                but2ene_z_seeds_passed = passes;
            }
            if passes == 5 {
                fully_passed_after.push(name);
            }
        }

        assert_eq!(
            but2ene_z_seeds_passed, 5,
            "but2ene_Z must now pass all 5 base seeds (was 1/5 before this fix)"
        );
        for &name in FULLY_PASSED_BEFORE_FIX {
            assert!(
                fully_passed_after.contains(&name),
                "{name} fully passed before this fix and must still fully pass after -- \
                 regression introduced by the declared-E/Z bound change"
            );
        }
    }

    /// Issue #291 residual (testosterone, one of the two ring-fused molecules
    /// `repair_tetrahedral_center`'s bridge-eligibility check alone cannot fix, since
    /// every real neighbor of the affected centers is a ring atom): confirms
    /// `materialize_implicit_h_for_chirality` combined with the post-refinement
    /// repair safety net in `try_embed_once` above actually delivers *correct*
    /// geometry, not just an internally-reported "satisfied". Checked against
    /// `chematic_chem::assign_cip` (declared) vs `stereo3d::assign_stereo_from_3d`
    /// (perceived from the *returned* coordinates) -- a genuinely different code
    /// path than this module's own `verify_stereo`, so it can't be fooled the same
    /// way a bug in `verify_stereo`'s own input (e.g. the now-fixed `add_hydrogens`
    /// stereo-order loss) could fool `verify_stereo` alone.
    ///
    /// Without the post-refinement safety net, every one of these 5 seeds failed
    /// with `StereoConstraintFailed` (measured directly while diagnosing this) even
    /// though `repair_stereo` alone, called directly on a raw embed, reached 6/6 on
    /// every seed -- `refine_coords`'s chirality-blind bound correction was undoing
    /// the pre-refinement repair every time.
    #[test]
    fn materialize_implicit_h_for_chirality_fixes_testosterone_with_correct_geometry() {
        let mol = parse("C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O").unwrap();
        let declared = chematic_chem::assign_cip(&mol);
        assert_eq!(
            declared.assignments.len(),
            6,
            "sanity: 6 declared stereocenters"
        );

        let mut passes = 0usize;
        for seed in 0..5u64 {
            let params = EmbedParameters {
                random_seed: seed,
                max_attempts: 8,
                enforce_chirality: true,
                materialize_implicit_h_for_chirality: true,
                ..EmbedParameters::default()
            };
            let coords = match embed_distance_geometry_v2_detail(&mol, &params) {
                Ok((coords, _)) => coords,
                Err(_) => continue,
            };
            // `assign_stereo_from_3d` only assigns centers with exactly 4 heavy-atom
            // neighbors (no implicit H) -- skip declared centers it can't independently
            // check (e.g. any with an implicit H), don't treat that as a mismatch.
            let perceived = crate::stereo3d::assign_stereo_from_3d(&mol, &coords);
            for &(idx, code) in &declared.assignments {
                if let Some(perceived_code) = perceived.get(idx) {
                    assert_eq!(
                        perceived_code, code,
                        "testosterone seed={seed}: atom {idx:?} declared {code:?} but \
                         3D-perceived {perceived_code:?} -- embed reported success with \
                         wrong chirality"
                    );
                }
            }
            passes += 1;
        }
        assert!(
            passes >= 4,
            "testosterone: expected at least 4/5 seeds to embed successfully with \
             materialize_implicit_h_for_chirality (post-refinement safety net), got {passes}/5"
        );
    }
}
