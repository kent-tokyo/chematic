# pipeline v2 vs RDKit ETKDGv3 — Wave 1 independent 3D benchmark

Measurement-only. No pipeline v2 or force-field algorithm code was changed to produce these numbers. Historical numbers are NOT reused -- everything below was regenerated fresh against this repo's current `main` in this session. All tables below are auto-generated from `validation/results/pipeline_v2_vs_rdkit_aggregate.json` by this script; the aggregate JSON is the source of truth if anything here looks stale.

## Corpus

- Tier A (curated stress): 65 molecules, sha256 `8d4f3f3a70b8ae00...`
- Tier B (fixed drug-like, ChEMBL-derived): 200 molecules, sha256 `0059fcbcb862663b...`
- Total: 265 molecules

## Atom mapping

- Checked: 265, verified matching: 265, unavailable/mismatched: 0

## Environment record (reproducibility)

- `benchmark_session`: release_grade_remeasure_v0_15_0_494d634
- `benchmark_commit`: 494d63437dfa6a5b840d669b7715367e0b7eb986
- `benchmark_date`: 2026-08-15
- `benchmark_branch`: docs/pipeline-v2-vs-rdkit-remeasure-494d634
- `common_scorer_blob_sha`: e7d33e9802cec01cc68814e977d410611e45b5b2
- `tier_a_manifest_sha256`: 6a478ea0f5d4ef067a4d1739e77a7209e8f76ecaa837e3f487c723dd6f465d6b
- `tier_b_manifest_sha256`: b3cde3fedcc68391ba3d3cbae228acd3057cadea6e6fc17499592b23bdc7550a
- `rdkit_version`: 2026.03.4
- `rust_version`: rustc 1.97.0 (2d8144b78 2026-07-07)
- `python_version`: 3.13.6
- `os_arch`: aarch64-apple-darwin

## Coverage and usable geometry (explicit denominators)

`usable_coverage` = independently-sound successes / total inputs for that arm -- the fraction of the *whole corpus* that arm turns into a geometry this benchmark's own independent scorer (not the pipeline's internal judgment) certifies sound. `sound_given_success` = independently-sound / successes only (the old, incomplete framing -- kept for context, never presented alone).

| Engine | Arm | total | success | indep. sound | sound_given_success | usable_coverage | typed_failure | unsupported | timeout | internal_error |
|---|---|---|---|---|---|---|---|---|---|---|
| chematic | chematic_pipeline_v2_no_ff | 265 | 254 | 254 | 100.0% | 95.8% | 11 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_dreiding | 265 | 254 | 254 | 100.0% | 95.8% | 11 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_uff_only | 265 | 250 | 250 | 100.0% | 94.3% | 15 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_strict | 265 | 241 | 241 | 100.0% | 90.9% | 12 | 12 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 254 | 254 | 100.0% | 95.8% | 11 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_strict_repair | 265 | 231 | 231 | 100.0% | 87.2% | 22 | 12 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 265 | 244 | 244 | 100.0% | 92.1% | 21 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 265 | 241 | 241 | 100.0% | 90.9% | 12 | 12 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 265 | 253 | 253 | 100.0% | 95.5% | 11 | 0 | 1 | 0 |
| chematic | chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 265 | 241 | 241 | 100.0% | 90.9% | 12 | 12 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 265 | 254 | 254 | 100.0% | 95.8% | 11 | 0 | 0 | 0 |
| chematic | chematic_legacy_etkdg | 265 | 265 | 248 | 93.6% | 93.6% | 0 | 0 | 0 | 0 |
| rdkit | rdkit_etkdgv3_raw | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |
| rdkit | rdkit_etkdgv3_uff | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |
| rdkit | rdkit_etkdgv3_mmff94 | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |
| rdkit | rdkit_etkdgv3_best_of_n | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |

**mmff94_strict, spelled out per the fix request:** 241/241 successful outputs are independently sound, but only 241/265 of the *total corpus* ends up as a usable geometry under this arm -- the rest is the 24-molecule MMFF94 coverage gap (issue #227), governed by this arm's own bond+angle coverage gate (`mmff94_strict` never gated stretch-bend, even before Priority 2B -- the Dfsb fallback's periodic-row default resolves every stretch-bend term in production now, 0 final-unresolved measured this run, see the Bonded-term coverage gate section below; NOT the same as this arm's *output* being unaffected -- Dfsb changes energy/gradient unconditionally for every MMFF94 arm, which can shift convergence/success outcomes even where gate eligibility doesn't move, see that section's own note on the one verified status change this run), not a geometry-quality problem.

## Common heavy-atom geometry quality (same independent scorer, both engines)

Applied identically to chematic's and RDKit's already-saved heavy-atom coordinates (`crates/chematic-3d/examples/pipeline_v2_vs_rdkit_common_scorer.rs`) -- ideal bond length from `Element::covalent_radius()`, never chematic-3d's own `pub(crate)` thresholds. RDKit's coordinates are heavy-atom-only by construction (the oracle script never exports its `AddHs`-added hydrogens).

| Engine | Arm | n scored | all finite | mean bond>15% | mean bond>50% | molecules w/ clash | molecules w/ coincident atoms | independently sound |
|---|---|---|---|---|---|---|---|---|
| chematic | chematic_pipeline_v2_no_ff | 254 | 100.0% | 2.8% | 0.0% | 3 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_dreiding | 254 | 100.0% | 3.0% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_uff_only | 250 | 100.0% | 0.5% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict | 241 | 100.0% | 0.9% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback | 254 | 100.0% | 0.8% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict_repair | 231 | 100.0% | 0.7% | 0.0% | 3 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 244 | 100.0% | 0.7% | 0.0% | 3 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 241 | 100.0% | 0.9% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 253 | 100.0% | 0.8% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 241 | 100.0% | 0.9% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 254 | 100.0% | 0.8% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_legacy_etkdg | 265 | 100.0% | 54.8% | 23.7% | 225 | 0 | 93.6% |
| rdkit | rdkit_etkdgv3_raw | 264 | 100.0% | 0.6% | 0.0% | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_uff | 264 | 100.0% | 0.0% | 0.0% | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_mmff94 | 264 | 100.0% | 0.9% | 0.0% | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_best_of_n | 264 | 100.0% | 0.0% | 0.0% | 0 | 0 | 100.0% |

This common scorer checks for exactly-coincident atom pairs (distance < 1e-3 Å), which the original ad-hoc legacy scorer did not -- 0/265 legacy outputs have ≥1 coincident atom pair and are NOT independently sound under this stricter, shared check. 11/11 pipeline_v2 arms are 100% independently sound this run (matching their own internal `final_validation.sound`); see the table above for any arm below 100%.

## Stereo preservation (same judge -- chematic's own `verify_stereo` -- applied to both engines)

**Methodology, read before the numbers**: the 9 `Ignore`-policy arms below (including the 4 bonded-term-gate arms added in Priority 2, which only change the coverage gate's scope, not stereo policy) reflect raw distance-geometry-embedding output -- `Ignore` never repairs a violated stereocenter, so those rows are NOT chematic's best achievable stereo correctness. Starting Priority 1 (v0.11.0 re-benchmark), 2 `StereoPolicy::RepairAndVerify` arms (`chematic_pipeline_v2_mmff94_strict_repair` / `..._with_uff_fallback_repair`) ARE exercised and shown below -- read those rows, not the Ignore rows, for chematic's best achievable stereo number under MMFF94. Their lower `declared`/`molecules w/ declared stereo` counts vs. the matching Ignore arm reflect fewer molecules reaching success at all under RepairAndVerify (see the RepairAndVerify effectiveness section below for the paired-arm accounting), not a smaller stereo-bearing subset by construction. RDKit's numbers use `enforceChirality=True` for real -- verified here with the identical judge, not assumed.

| Engine | Arm | molecules w/ declared stereo | declared | satisfied | violated | unevaluable | satisfaction rate |
|---|---|---|---|---|---|---|---|
| chematic | chematic_pipeline_v2_no_ff | 83 | 146 | 82 | 64 | 0 | 56.2% |
| chematic | chematic_pipeline_v2_dreiding | 83 | 146 | 89 | 57 | 0 | 61.0% |
| chematic | chematic_pipeline_v2_uff_only | 80 | 140 | 88 | 52 | 0 | 62.9% |
| chematic | chematic_pipeline_v2_mmff94_strict | 83 | 146 | 82 | 64 | 0 | 56.2% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback | 83 | 146 | 82 | 64 | 0 | 56.2% |
| chematic | chematic_pipeline_v2_mmff94_strict_repair | 73 | 112 | 112 | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 73 | 112 | 112 | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 83 | 146 | 82 | 64 | 0 | 56.2% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 83 | 146 | 82 | 64 | 0 | 56.2% |
| chematic | chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 83 | 146 | 82 | 64 | 0 | 56.2% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 83 | 146 | 82 | 64 | 0 | 56.2% |
| chematic | chematic_legacy_etkdg | 90 | 170 | 107 | 63 | 0 | 62.9% |
| rdkit | rdkit_etkdgv3_raw | 90 | 170 | 170 | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_uff | 90 | 170 | 170 | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_mmff94 | 90 | 170 | 170 | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_best_of_n | 90 | 170 | 170 | 0 | 0 | 100.0% |

`violated` encompasses both tetrahedral inversion and E/Z flipping (both fail the declared-direction check `verify_stereo` performs) -- the shared judge does not currently distinguish these as separate sub-categories, so this report doesn't either, rather than fabricate a split it can't measure.

## Workflow comparison vs. common heavy-atom output comparison

**Workflow comparison** (each library's own recommended, practical usage: RDKit with `AddHs`, chematic's implicit-H pipeline as-is): this is what the Performance section's wall-clock numbers below measure. Not an algorithm-only, hydrogen-representation-controlled comparison.

**Common heavy-atom output comparison**: the geometry-quality and stereo tables above restrict to heavy atoms only on both sides, via the identical scorer/judge, so differing internal hydrogen treatment cannot bias the output-quality numbers.

An RDKit `AddHs=false` auxiliary arm was NOT added this round (would meaningfully grow the arm matrix) -- performance numbers below should be read as workflow-level, not algorithm-only apples-to-apples.

## Performance

### Process-level performance: NOT RUN this round -- the chematic arm matrix has grown from the `1bc1b63`-era 6 to 12 (2 RepairAndVerify arms added in Priority 1, 4 bonded-term-gate arms added in Priority 2), so the stored `1bc1b63`-era process-level file would no longer be measuring the same binary and was deliberately excluded rather than presented as if comparable. In-process per-(molecule, arm) timing below is the primary comparable metric this round. Re-run `scripts/pipeline_v2_vs_rdkit_process_level_perf.sh` in a follow-up if the whole-corpus process-level figure is needed against the new arm matrix.

### In-process per-(molecule, arm) timing (secondary)

_In-process wall-clock per (molecule, arm) call within a single long-running process -- NOT process-isolated. Secondary metric; see performance_process_level for the primary comparison._

#### chematic

| Arm | n | p50 (ms) | p95 (ms) | p99 (ms) | max (ms) |
|---|---|---|---|---|---|
| chematic_pipeline_v2_no_ff | 265 | 4.0 | 25.8 | 37.1 | 48 |
| chematic_pipeline_v2_dreiding | 265 | 134.0 | 1036.8 | 1288.8 | 1692 |
| chematic_pipeline_v2_uff_only | 265 | 209.0 | 1462.0 | 2280.4 | 3508 |
| chematic_pipeline_v2_mmff94_strict | 265 | 1047.0 | 6452.4 | 9022.7 | 13032 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 1068.0 | 6383.4 | 9475.7 | 11623 |
| chematic_pipeline_v2_mmff94_strict_repair | 265 | 942.0 | 6340.6 | 10127.2 | 17213 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 265 | 958.0 | 6480.8 | 11366.8 | 17476 |
| chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 265 | 1047.0 | 6672.8 | 10735.3 | 17244 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 265 | 1065.0 | 6643.2 | 10758.6 | 21946 |
| chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 265 | 1062.0 | 6646.8 | 11812.2 | 19154 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 265 | 1067.0 | 6357.8 | 12855.3 | 19345 |
| chematic_legacy_etkdg | 265 | 1.0 | 4.0 | 6.4 | 11 |

#### RDKit

| Arm | n | p50 (ms) | p95 (ms) | p99 (ms) | max (ms) |
|---|---|---|---|---|---|
| rdkit_etkdgv3_raw | 265 | 9.0 | 78.0 | 144.2 | 346 |
| rdkit_etkdgv3_uff | 265 | 19.0 | 113.0 | 189.4 | 413 |
| rdkit_etkdgv3_mmff94 | 265 | 22.0 | 120.0 | 197.0 | 440 |
| rdkit_etkdgv3_best_of_n | 265 | 190.0 | 1089.2 | 2062.0 | 3102 |

## Cyclopentane RDKit crash — scoped ablation

**Classification: `nondefault_small_ring_torsion_only`**

12/60 trials crashed. Crashing configs (`useSmallRingTorsions`, `enforceChirality`): ['(True, False)', '(True, True)']. Crashing seeds: [4, 20260801]. Crashes under RDKit's own default config (`useSmallRingTorsions=False`): 0.

In plain terms: this crash requires the non-default `useSmallRingTorsions=True`, occurs during `EmbedMolecule` itself (before any force-field stage runs), and only reproduces for a subset of tested seeds -- **not** a general "RDKit crashes on cyclopentane" finding, and not reproducible under RDKit's own ETKDGv3 defaults in this ablation. Minimal repro: `scripts/pipeline_v2_vs_rdkit_cyclopentane_crash_ablation.py`.

## Force-field coverage (chematic MMFF94 arms)

- chematic_pipeline_v2_mmff94_with_uff_fallback: n=254, fallback_rate=5.1%, converged_rate=21.7%
- chematic_pipeline_v2_mmff94_strict: n=241, fallback_rate=0.0%, converged_rate=22.4%
- chematic_pipeline_v2_mmff94_strict_repair: n=231, fallback_rate=0.0%, converged_rate=22.5%
- chematic_pipeline_v2_mmff94_with_uff_fallback_repair: n=244, fallback_rate=5.3%, converged_rate=21.7%
- chematic_pipeline_v2_mmff94_strict_stretch_bend_gated: n=241, fallback_rate=0.0%, converged_rate=22.4%
- chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated: n=253, fallback_rate=5.1%, converged_rate=21.7%
- chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated: n=241, fallback_rate=0.0%, converged_rate=22.4%
- chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated: n=254, fallback_rate=5.1%, converged_rate=21.7%

## Stage funnel (per-arm denominator hierarchy)

Real `pipeline_v2` execution order (`crates/chematic-3d/src/pipeline_v2.rs` `PipelineStage` enum + its actual call sequence): embed (`DistanceGeometry`) -> torsion optimization -> **stereo verify/repair** -> force-field minimization -> final stereo verify -> final geometry validation. Stereo repair happens *before* force-field minimization, not after -- the columns below follow that real order, not an assumed embed-then-FF-then-stereo sequence. A row is counted under an `_attempted`/`_reached` column if its `failure_stage` is that stage or later (or it succeeded outright); under a `_succeeded`/`_verified` column only if `failure_stage` is strictly later than that stage (or it succeeded outright) -- a row that failed AT a stage reached it but did not pass it, so `ff_attempted` and `ff_succeeded` are genuinely different counts, not the same check twice. Never collapsed into a single success rate -- see `feedback_fallback_pooling_measurement_error`: `mmff94_strict` and `mmff94_with_uff_fallback` are reported as fully separate rows, never blended.

| Arm | attempted | embed_succeeded | stereo_repair_reached | ff_attempted | ff_succeeded | final_stereo_verified | final_validation_passed |
|---|---|---|---|---|---|---|---|
| chematic_pipeline_v2_no_ff | 265 | 254 | 254 | 254 | 254 | 254 | 254 |
| chematic_pipeline_v2_dreiding | 265 | 254 | 254 | 254 | 254 | 254 | 254 |
| chematic_pipeline_v2_uff_only | 265 | 254 | 254 | 254 | 250 | 250 | 250 |
| chematic_pipeline_v2_mmff94_strict | 265 | 254 | 254 | 254 | 241 | 241 | 241 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 254 | 254 | 254 | 254 | 254 | 254 |
| chematic_pipeline_v2_mmff94_strict_repair | 265 | 254 | 254 | 244 | 231 | 231 | 231 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 265 | 254 | 254 | 244 | 244 | 244 | 244 |
| chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 265 | 254 | 254 | 254 | 241 | 241 | 241 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 265 | 254 | 254 | 254 | 253 | 253 | 253 |
| chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 265 | 254 | 254 | 254 | 241 | 241 | 241 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 265 | 254 | 254 | 254 | 254 | 254 | 254 |
| chematic_legacy_etkdg | 265 | n/a | n/a | n/a | n/a | n/a | 265 |

Note: `chematic_legacy_etkdg` does not run through `pipeline_v2` at all (separate `generate_coords_etkdg` entry point, no `PipelineStage` tracking) -- its row reports `attempted`/`final_validation_passed` only; the intermediate columns are `n/a` rather than a fabricated 0 or a misleading 265 (a naive reuse of the success-implies-passed-every-stage rule above would have printed 265 for every column here, which would misrepresent a code path that never runs those stages at all).

## RepairAndVerify effectiveness (paired-arm comparison, Priority 1 new arms)

Each `StereoPolicy::RepairAndVerify` arm is a genuinely independent arm (not a config edit to the pre-existing `Ignore`-policy arm of the same `ForceFieldPolicy`), paired here per-molecule against its Ignore counterpart. `repair_time_delta` is `(repair-arm elapsed_ms - Ignore-arm elapsed_ms)` per molecule -- a paired-arm difference, not a directly instrumented repair-stage timer (pipeline_v2 does not currently expose one).

**Why `before-mismatch`/`repair attempted`/`repair succeeded` are identical between the two repair arms below**: verified directly (not assumed) -- `stereo_before_violations` matches per-molecule, 1:1, across `chematic_pipeline_v2_mmff94_strict_repair` and `chematic_pipeline_v2_mmff94_with_uff_fallback_repair` for all 265 molecules. This is structural, not a bug: stereo verify/repair runs BEFORE force-field minimization in `pipeline_v2`'s real execution order (see the stage-funnel note above), so both arms see the identical pre-FF geometry and make identical repair decisions -- the two `ForceFieldPolicy` values can only diverge afterward. `after-mismatch` DOES differ (11 vs. 13) because `final_stereo` is measured after FF minimization, where the two force fields' behavior can differ.

| Ignore arm | Repair arm | n compared | excluded (incomparable) | before-mismatch | repair attempted | repair succeeded | outcome unavailable | after-mismatch | geometry pairs | geometry degraded | time delta median (ms) | time delta p95 (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chematic_pipeline_v2_mmff94_strict | chematic_pipeline_v2_mmff94_strict_repair | 254 | 11 | 46 | 64 | 52 | 0 | 0 | 231 | 0 | 0 | 212 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 254 | 11 | 46 | 64 | 52 | 0 | 0 | 244 | 0 | 0 | 122 |

_Note on reading `after-mismatch` next to the stereo-preservation table's "100% satisfaction" figure below: they are measured over different populations. The 100% satisfaction rate is computed only over the arm's *successful* rows (131/229); `after-mismatch` here is computed over all 254 comparable rows, including ones that failed after repair (e.g. in force-field minimization) -- so a non-zero after-mismatch count and a 100%-among-successes satisfaction rate are not in conflict, they answer different denominators. See the Stage funnel table above for the exact per-arm counts at each stage._

## Bonded-term coverage gate (Priority 2 / 2B / Stage 1B, issue #227)

**Priority 2** added `gate_mmff94_stretch_bend` (`PipelineV2Config`/`minimize_with_policy_gated`) — an opt-in gate refusing `Mmff94BondAngleStrict`/`Mmff94WithUffFallback` on a missing stretch-bend term — plus 4 benchmark arms exercising a 3-stage comparison (legacy -> stretch-bend-gated -> complete-bonded-term-gated, for both `mmff94_strict` and `mmff94_with_uff_fallback`; "complete-bonded-term", not "complete MMFF94" -- vdW/partial-charge coverage are never gated by any arm here). That round's own diagnostic audit found the single largest missing-term bucket (StretchBend, 2,107 instances total: 1,680 genuine `table_gap` + 427 `routing_bug_candidate`) was **100% coverable** by porting a small, pinned-RDKit-commit-verified 29-row periodic-table-row fallback table (`MMFFDfsbCollection`'s real RDKit equivalent) into production.

**Priority 2B (this round) ships that port.** `chematic_ff::mmff94_stbn` now tries the existing specific/generic MMFF-type table first (unchanged, always wins if it has a row), and on failure falls back to the ported Dfsb table — **unconditionally, not behind any opt-in flag** (this is a production accuracy fix, not a diagnostic feature; it applies to every MMFF94 policy's energy/gradient calculation, and to the coverage gate the same way). The `gate_mmff94_stretch_bend`/`gate_mmff94_torsion_oop` *strict-refusal* gates from Priority 2 are unaffected by this and remain independent opt-ins, still `false` by default — Priority 2B only changes what counts as "covered" underneath those gates, not whether the gates themselves are on.

**Coverage parity achieved (0/2,107 final-unresolved), but this is NOT the same as parameter-selection parity.** Two structurally different outcomes hide behind the same "resolved" status:

- **1,680 instances** were genuine table gaps (absent at *every* classification code chematic-ff's tables define) -- Dfsb resolving these matches RDKit's own real behavior exactly. This IS the case Dfsb was built to close.
- **427 instances** were routing-bug candidates (a real, correctly-typed parameter already exists at a *different* classification code than the one this molecule's context computed) that Dfsb *also* happens to resolve. Coverage is achieved, but chematic is now using RDKit's generic periodic-row default instead of the specific parameter a correctly-routed classification would have used -- **masked, not fixed**. Before this fix, these instances were reported as "missing"; after, they silently look identical to genuinely-resolved instances unless this breakdown is consulted.

This distinction is preserved in `mmff94_term_coverage_audit.rs`'s own output specifically so it doesn't disappear: the audit emits a row whenever the TYPE-ONLY lookup misses, regardless of whether Dfsb then rescues it, with a `dfsb_resolved` field and the original `present_at_different_classification` discriminator both preserved on the same row. **Not fixed in this PR** (would require investigating `angle_type_for`'s classification logic, a different root cause than the Dfsb port -- tracked as follow-up work, not silently dropped).

### Missing-term sub-classification (fresh re-run, `mmff94_term_coverage_audit.rs`)

Per-term-instance classification across the 265-molecule corpus, using the TYPE-ONLY lookup (`mmff94_stbn_type_only` for StretchBend) -- independent of whether production `mmff94_stbn`'s Dfsb fallback ultimately resolves a given instance (see the coverage-vs-parameter-selection-parity note above for StretchBend specifically). `routing_bug_candidate` = this exact atom-type tuple has a table row at a *different* classification code than the one this molecule's context computed -- a candidate for an `angle_type_for`/`torsion_type_for`/`bond_type_for` classification bug, not necessarily a genuine table gap. `table_gap` = absent at *every* classification code chematic-ff's tables define. `Oop` is listed explicitly even at 0 -- omitting a measured-zero term kind would be indistinguishable from "not measured", which it is not.

| Term kind | total missing instances | routing_bug_candidate | table_gap |
|---|---|---|---|
| Bond | 84 | 79 (94.0%) | 5 (6.0%) |
| Angle | 374 | 277 (74.1%) | 97 (25.9%) |
| Torsion | 1121 | 1107 (98.8%) | 14 (1.2%) |
| Oop | 0 | 0 (n/a) | 0 (n/a) |
| StretchBend | 2107 | 427 (20.3%) | 1680 (79.7%) |

For Bond/Angle/Torsion/Oop, `table_gap` is not further sub-classified -- chematic-ff implements neither MMFF94 equivalence-class substitution nor empirical-rule (e.g. Badger's-rule bond) estimation at all for these term kinds (`Mmff94NumericTypeInfo.equivalence_levels` carries real MMFF94 equivalence data but has zero readers anywhere in the codebase, verified, not assumed) -- deferred, not fabricated. Per this round's explicit scope decision, this PR does not touch those routing-bug candidates either (nor StretchBend's own 427 masked routing candidates above), to keep a single root cause (Dfsb port only).

### Legacy -> stretch-bend -> complete-bonded-term (3-stage paired comparison)

Each stage's arm is a genuinely independent arm (never a config edit to a previous stage's arm), compared per-molecule against the immediately preceding stage. For `mmff94_strict` (pure gate, no fallback), widening the gate can only ever turn a prior success into a failure, never the reverse -- verified as a hard invariant at generation time (both stage transitions), not just a display column. This is NOT a hard invariant for `mmff94_with_uff_fallback`: see the note below the table.

| Earlier stage | Later stage | n compared | earlier success | later success | newly failing |
|---|---|---|---|---|---|
| chematic_pipeline_v2_mmff94_strict | chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 265 | 241 | 241 | 0 |
| chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 265 | 241 | 241 | 0 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 265 | 254 | 253 | 1 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 265 | 253 | 254 | 0 |

`chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated` newly-failing molecules (1): chembl_tier_b_0030

`chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated` **also has 1 molecule(s) that flip the other way** (earlier stage fails, later stage succeeds) -- independently verified against one of two recognized explanation categories (`GATE_WIDENING_EXPLANATIONS` in `scripts/gen_pipeline_v2_vs_rdkit_report.py`), never asserted away. Any case matching neither category fails report generation instead of being silently accepted.

- **identical_coverage_timing_variance** (1): chembl_tier_b_0030. The later row succeeded via the SAME (non-fallback) force field, with the coverage dimension this stage newly gates on measured as fully covered (zero missing) -- the gate flag is a mechanical no-op for these molecules, so both stages run identical minimization on identical geometry; the timeout/success flip is consistent with wall-clock scheduling variance around the shared 20000ms boundary, corroborated by the later row's own `elapsed_ms` being a substantial fraction of that budget (not an arbitrary fast, unrelated success). First observed on `chembl_tier_b_0030` during A1's 265-molecule re-audit (see `docs/rfcs/a1_conformer_benchmark_failure_ledger.md`, Finding 3).

**Stretch-bend-gated -> complete-bonded-term is 0 newly-failing for every policy this round** -- every molecule that already survives the stretch-bend gate also has complete torsion+OOP coverage in this specific 265-molecule corpus. This is empirical, not structural (torsion has 1,121 missing instances measured above; they evidently concentrate on molecules that already fail the stretch-bend gate, in this corpus) -- a different, larger, or differently-composed corpus could show a non-zero delta at this stage. Practical effect for this run: the stretch-bend-gated and complete-bonded-term-gated success counts are numerically identical here, so the corrected, narrower name (`..._stretch_bend_gated`, not "true complete-term") only matters for what the number *means*, not for its value on this particular corpus.

**On the legacy `mmff94_strict`/`mmff94_with_uff_fallback` arms' success counts changing at all (148->149, and a similar 1-2 molecule shift seen in earlier rounds): this is NOT structurally guaranteed to be zero, and is not claimed to be.** `mmff94_strict` never gated stretch-bend, before or after Priority 2B -- bond+angle *gate eligibility* is unchanged -- but Priority 2B changes stretch-bend's contribution to every MMFF94 policy's energy AND finite-difference gradient *unconditionally*, for every molecule that reaches minimization under any policy. That can change minimizer convergence, iteration count, final residual force, and therefore final soundness/success -- in principle for better or worse, not just gate-count-preserving by construction. What this round actually measured, with a real per-molecule diff against a baseline saved *before* re-running (not reconstructed after the fact): 0 soundness regressions among molecules sound in both runs, and exactly 1 status change on `mmff94_strict` -- `chembl_tier_b_0166` (elapsed_ms 20530 -> 16221, status timeout -> success). `embed_seed` governs geometry/RNG determinism but not real-time scheduling, so a molecule sitting near the `total_timeout_ms=20000` boundary is a plausible site for this kind of flip regardless of cause -- but a ~4.3s drop is a substantial, consistent-direction change, not obviously pure machine-load noise, and is reported here as verified-but-not-fully-explained rather than asserted to be "known jitter" without checking. The same molecule ID was *also* the timeout-boundary case in Priority 2's own `mmff94_with_uff_fallback` measurement -- a recurring boundary case across multiple rounds, consistent with a genuinely ~20s-class molecule under this policy family, sensitive to any change in computation.

**Scope of "adopted" this round**: the Dfsb periodic-row fallback itself (Priority 2B) IS now unconditional production behavior for every MMFF94 policy's energy/gradient calculation and coverage measurement -- not gated, not opt-in. What remains opt-in and `false` by default is the *strict-refusal* gate on top of that coverage (`gate_mmff94_stretch_bend`/`gate_mmff94_torsion_oop`). This changes what energy/gradient every MMFF94 arm computes, not just what a gate refuses -- see the paragraph above for why that is a real, if empirically small this round, source of output change even for arms whose gate eligibility never touches stretch-bend.

## Ring-torsion FailClosed probe

1 row(s) -- demonstrates `RingTorsionApplicationPolicy::FailClosed`'s documented behavior. Not folded into any of the 12 main arms' coverage numbers (those use `DiagnosticOnly`).

## Reference geometry subset

Status: **insufficient_evidence**. No experimentally-determined reference conformers were available for this benchmark round. RMSD-vs-reference, best-of-N RMSD, torsion fingerprint deviation, and duplicate-conformer-rate metrics are NOT computed here -- reported as insufficient evidence, not fabricated.

## Known issues filed from this benchmark

- MMFF94 coverage gap (24/265 not successful under mmff94_strict, PR #236/#238/#239/#241 fixes already reflected in this run): https://github.com/kent-tokyo/chematic/issues/227

## Data integrity

- Unclassified rows: 0 (hard-gated at 0 by the report generator)
- chematic rows sha256: `afd098159d6e9754...`
- RDKit rows sha256: `4c4e02540ede2d59...`
- All integrity gates (row-count, unclassified, atom-mapping, missing/mismatched coords, non-finite coords, common-scorer coverage, denominator self-consistency) passed at generation time -- see `run_integrity_gates` in this script.

## Conclusions

Classified per class/metric — no single overall win/loss score.

| Metric | Classification | Basis |
|---|---|---|
| Coverage — no_ff/dreiding/uff_only/mmff94_with_uff_fallback vs. RDKit | Roughly comparable | chematic 94.3%-95.8% success vs. RDKit 99.6% |
| Coverage — mmff94_strict | RDKit-favor (chematic gap, issue #227 filed) | 90.9% success, 12/265 unsupported |
| Common heavy-atom geometry — pipeline_v2 force-field arms | Chematic strength on soundness | 100% independently-sound across dreiding/uff_only/mmff94 arms, matching pipeline-internal judgment |
| Common heavy-atom geometry — legacy etkdg | Known gap, refined this round | 14/265 legacy outputs have coincident atoms under the stricter common scorer (not caught by the original Wave 1 ad-hoc check); the already-documented clash-rate gap stands |
| Stereo preservation (same judge, `Ignore`) | RDKit-favor | RDKit 100.0% satisfaction vs. chematic 62.9% under `StereoPolicy::Ignore` -- not chematic's best achievable number, see next row |
| Stereo preservation (same judge, `RepairAndVerify`, new this round) | Parity with RDKit among successes, coverage gap remains the real cost | mmff94_strict_repair 100.0%, mmff94_with_uff_fallback_repair 100.0% satisfaction among molecules that reached success under RepairAndVerify (both match RDKit's 100% on that subset) -- but RepairAndVerify also reduces the success *count* vs. the matching Ignore arm (fewer molecules reach final success at all when repair is required to pass); see the RepairAndVerify effectiveness section for the exact paired accounting |
| Bonded-term coverage gate, mmff94_strict -> mmff94_strict_stretch_bend_gated (new this round) | Real coverage gap surfaced, widening the gate is a real cost | 241 earlier-stage successes -> 241 under the later stage's gate (0 newly fail, 0.0% of earlier-stage successes) -- see the Bonded-term coverage gate section for the term-kind sub-classification and full molecule list |
| Bonded-term coverage gate, mmff94_strict_stretch_bend_gated -> mmff94_strict_complete_bonded_term_gated (new this round) | Real coverage gap surfaced, widening the gate is a real cost | 241 earlier-stage successes -> 241 under the later stage's gate (0 newly fail, 0.0% of earlier-stage successes) -- see the Bonded-term coverage gate section for the term-kind sub-classification and full molecule list |
| Bonded-term coverage gate, mmff94_with_uff_fallback -> mmff94_with_uff_fallback_stretch_bend_gated (new this round) | Real coverage gap surfaced, widening the gate is a real cost | 254 earlier-stage successes -> 253 under the later stage's gate (1 newly fail, 0.4% of earlier-stage successes) -- see the Bonded-term coverage gate section for the term-kind sub-classification and full molecule list |
| Bonded-term coverage gate, mmff94_with_uff_fallback_stretch_bend_gated -> mmff94_with_uff_fallback_complete_bonded_term_gated (new this round) | Real coverage gap surfaced, widening the gate is a real cost | 253 earlier-stage successes -> 254 under the later stage's gate (0 newly fail, 0.0% of earlier-stage successes) -- see the Bonded-term coverage gate section for the term-kind sub-classification and full molecule list |
| Force-field convergence rate | RDKit-favor, and an input to Priority 3 (Stage 1C) | chematic mmff94_with_uff_fallback 21.7% converged within 200 iterations, yet 254/265 of that arm's runs pass final validation regardless -- i.e. most successful outputs did NOT converge within 200 iterations and still passed geometry validation. Either `force_field_converged` is narrower than "produced a usable geometry" (an iteration-budget artifact, not necessarily a quality problem), or this is a real gap worth diagnosing -- Priority 3's MinimizationFailed root-causing (CatastrophicBondBlowup vs. ExcessiveResidualForce) is the next stage that should resolve which; corroborates open issues #185/#188 |
| Known crashes | RDKit has a narrowly-scoped one; chematic none found this round | cyclopentane crash classified `nondefault_small_ring_torsion_only` -- non-default config, seed-dependent, not RDKit's own default behavior |
| Unsupported chemistry | RDKit-favor | chematic mmff94_strict 24/265 unsupported (issue #227); RDKit's 4 arms show 0 unsupported_chemistry rows |
| Reference-geometry accuracy / torsion fingerprint / conformer diversity | Insufficient evidence | not measured this round, not fabricated |
| Overall "does chematic beat RDKit" | Not claimed | per this program's explicit rule -- findings are class/metric-specific |

