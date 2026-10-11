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

#[allow(unused_imports)]
use crate::fmath::DetMath;
use chematic_core::{AtomIdx, Molecule};

use crate::clock::Instant;
use crate::coords::Coords3D;
use crate::dg_fft::ideal_bond_length;
use crate::distance_geometry_v2::{
    self, BoundsConformance, DistanceBoundAdjustment, EmbedFailureCause, EmbedParameters,
    EmbedStats, EmbedWithAdjustmentsFailure, bounds_conformance, truncate_coords,
};
use crate::etkdg_knowledge::{
    PairBoundAdjustment, TorsionKnowledgeConfig, TorsionKnowledgeError, TorsionKnowledgeReport,
    TorsionKnowledgeSource, TorsionOptimizationConfig, TorsionOptimizationReport, TorsionPotential,
    build_torsion_knowledge, evaluate_torsion_energy, is_bridge_bond,
    macrocycle_14_bound_adjustments, optimize_torsions,
};
use crate::minimize::{
    ForceFieldBridgeError, ForceFieldPolicy, MAX_SANE_BOND_LENGTH, MinimizeConfig, Mmff94Minimizer,
    PolicyMinimizeResult, minimize_with_policy_gated_using,
    minimize_with_policy_gated_with_constraint,
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

include!("pipeline_v2/types.rs");
include!("pipeline_v2/execute.rs");
include!("pipeline_v2/validation.rs");

#[cfg(test)]
#[path = "pipeline_v2_tests.rs"]
mod embedded_tests;
