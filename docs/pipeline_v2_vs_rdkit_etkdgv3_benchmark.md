# pipeline v2 vs RDKit ETKDGv3 — Wave 1 independent 3D benchmark

Measurement-only. No pipeline v2 or force-field algorithm code was changed to produce these numbers. Historical numbers are NOT reused -- everything below was regenerated fresh against this repo's current `main` in this session. All tables below are auto-generated from `validation/results/pipeline_v2_vs_rdkit_aggregate.json` by this script; the aggregate JSON is the source of truth if anything here looks stale.

## Corpus

- Tier A (curated stress): 65 molecules, sha256 `8d4f3f3a70b8ae00...`
- Tier B (fixed drug-like, ChEMBL-derived): 200 molecules, sha256 `1a1698f5444c1b6e...`
- Total: 265 molecules

## Atom mapping

- Checked: 265, verified matching: 265, unavailable/mismatched: 0

## Environment record (reproducibility)

- `benchmark_session`: priority2_stretch_bend_gate_227_rerun_review_fixes
- `benchmark_commit`: 1a1366b (+ this commit's review-driven blocker fixes: real complete-bonded-term arms, Dfsb sub-classification, strengthened timeout-rescue assertion, WASM/Python API compat)
- `benchmark_date`: 2026-08-05
- `benchmark_branch`: fix/mmff94-stretch-bend-gate-227-priority2
- `common_scorer_blob_sha`: 49a8854c0ff493e22b109f69d893d08662e45326
- `tier_a_manifest_sha256`: 8d4f3f3a70b8ae00c56cdbea81e398093ef25d896f4f252e46448632632166b8
- `tier_b_manifest_sha256`: 1a1698f5444c1b6e32da4eaa896cc347c7f051fdf71044bf53b6a644f0b9af77
- `rdkit_version`: 2026.03.3
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
| chematic | chematic_pipeline_v2_mmff94_strict | 265 | 148 | 148 | 100.0% | 55.8% | 12 | 104 | 1 | 0 |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 252 | 252 | 100.0% | 95.1% | 13 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_strict_repair | 265 | 131 | 131 | 100.0% | 49.4% | 33 | 101 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 265 | 229 | 229 | 100.0% | 86.4% | 36 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 265 | 37 | 37 | 100.0% | 14.0% | 12 | 216 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 265 | 250 | 250 | 100.0% | 94.3% | 15 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 265 | 37 | 37 | 100.0% | 14.0% | 12 | 216 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 265 | 250 | 250 | 100.0% | 94.3% | 15 | 0 | 0 | 0 |
| chematic | chematic_legacy_etkdg | 265 | 265 | 251 | 94.7% | 94.7% | 0 | 0 | 0 | 0 |
| rdkit | rdkit_etkdgv3_raw | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |
| rdkit | rdkit_etkdgv3_uff | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |
| rdkit | rdkit_etkdgv3_mmff94 | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |
| rdkit | rdkit_etkdgv3_best_of_n | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |

**mmff94_strict, spelled out per the fix request:** 148/148 successful outputs are independently sound, but only 148/265 of the *total corpus* ends up as a usable geometry under this arm -- the rest is the 117-molecule MMFF94 coverage gap (issue #227, 2,107 stretch-bend terms still ungated by *this* arm (a new, independent `gate_mmff94_stretch_bend=true` opt-in exists as of Priority 2/Stage 1B -- see that section below -- but is not adopted as this arm's default)), not a geometry-quality problem.

## Common heavy-atom geometry quality (same independent scorer, both engines)

Applied identically to chematic's and RDKit's already-saved heavy-atom coordinates (`crates/chematic-3d/examples/pipeline_v2_vs_rdkit_common_scorer.rs`) -- ideal bond length from `Element::covalent_radius()`, never chematic-3d's own `pub(crate)` thresholds. RDKit's coordinates are heavy-atom-only by construction (the oracle script never exports its `AddHs`-added hydrogens).

| Engine | Arm | n scored | all finite | mean bond>15% | mean bond>50% | molecules w/ clash | molecules w/ coincident atoms | independently sound |
|---|---|---|---|---|---|---|---|---|
| chematic | chematic_pipeline_v2_no_ff | 254 | 100.0% | 2.8% | 0.0% | 3 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_dreiding | 254 | 100.0% | 3.0% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_uff_only | 250 | 100.0% | 0.5% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict | 148 | 100.0% | 0.7% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback | 252 | 100.0% | 0.7% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict_repair | 131 | 100.0% | 0.5% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 229 | 100.0% | 0.4% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 37 | 100.0% | 0.0% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 250 | 100.0% | 0.5% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 37 | 100.0% | 0.0% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 250 | 100.0% | 0.5% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_legacy_etkdg | 265 | 100.0% | 49.8% | 15.5% | 229 | 14 | 94.7% |
| rdkit | rdkit_etkdgv3_raw | 264 | 100.0% | 0.6% | 0.0% | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_uff | 264 | 100.0% | 0.0% | 0.0% | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_mmff94 | 264 | 100.0% | 0.9% | 0.0% | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_best_of_n | 264 | 100.0% | 0.0% | 0.0% | 0 | 0 | 100.0% |

This common scorer checks for exactly-coincident atom pairs (distance < 1e-3 Å), which the original ad-hoc legacy scorer did not -- 14/265 legacy outputs have ≥1 coincident atom pair and are NOT independently sound under this stricter, shared check. 11/11 pipeline_v2 arms are 100% independently sound this run (matching their own internal `final_validation.sound`); see the table above for any arm below 100%.

## Stereo preservation (same judge -- chematic's own `verify_stereo` -- applied to both engines)

**Methodology, read before the numbers**: the 9 `Ignore`-policy arms below (including the 4 bonded-term-gate arms added in Priority 2, which only change the coverage gate's scope, not stereo policy) reflect raw distance-geometry-embedding output -- `Ignore` never repairs a violated stereocenter, so those rows are NOT chematic's best achievable stereo correctness. Starting Priority 1 (v0.11.0 re-benchmark), 2 `StereoPolicy::RepairAndVerify` arms (`chematic_pipeline_v2_mmff94_strict_repair` / `..._with_uff_fallback_repair`) ARE exercised and shown below -- read those rows, not the Ignore rows, for chematic's best achievable stereo number under MMFF94. Their lower `declared`/`molecules w/ declared stereo` counts vs. the matching Ignore arm reflect fewer molecules reaching success at all under RepairAndVerify (see the RepairAndVerify effectiveness section below for the paired-arm accounting), not a smaller stereo-bearing subset by construction. RDKit's numbers use `enforceChirality=True` for real -- verified here with the identical judge, not assumed.

| Engine | Arm | molecules w/ declared stereo | declared | satisfied | violated | unevaluable | satisfaction rate |
|---|---|---|---|---|---|---|---|
| chematic | chematic_pipeline_v2_no_ff | 83 | 146 | 82 | 64 | 0 | 56.2% |
| chematic | chematic_pipeline_v2_dreiding | 83 | 146 | 89 | 57 | 0 | 61.0% |
| chematic | chematic_pipeline_v2_uff_only | 80 | 140 | 88 | 52 | 0 | 62.9% |
| chematic | chematic_pipeline_v2_mmff94_strict | 63 | 103 | 55 | 48 | 0 | 53.4% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback | 82 | 143 | 81 | 62 | 0 | 56.6% |
| chematic | chematic_pipeline_v2_mmff94_strict_repair | 46 | 64 | 64 | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 59 | 85 | 85 | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 9 | 18 | 10 | 8 | 0 | 55.6% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 80 | 140 | 88 | 52 | 0 | 62.9% |
| chematic | chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 9 | 18 | 10 | 8 | 0 | 55.6% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 80 | 140 | 88 | 52 | 0 | 62.9% |
| chematic | chematic_legacy_etkdg | 90 | 170 | 81 | 86 | 3 | 47.6% |
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
| chematic_pipeline_v2_no_ff | 265 | 6.0 | 45.4 | 81.3 | 132 |
| chematic_pipeline_v2_dreiding | 265 | 191.0 | 1712.0 | 2241.5 | 2680 |
| chematic_pipeline_v2_uff_only | 265 | 339.0 | 2522.4 | 3561.1 | 7455 |
| chematic_pipeline_v2_mmff94_strict | 265 | 46.0 | 9593.4 | 13553.7 | 20530 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 835.0 | 9988.4 | 13326.5 | 18710 |
| chematic_pipeline_v2_mmff94_strict_repair | 265 | 26.0 | 9605.4 | 13848.8 | 17601 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 265 | 610.0 | 9587.4 | 13213.8 | 18287 |
| chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 265 | 9.0 | 155.4 | 455.9 | 3134 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 265 | 349.0 | 2398.6 | 3643.8 | 6185 |
| chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 265 | 10.0 | 119.4 | 479.1 | 2721 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 265 | 336.0 | 2405.6 | 3463.5 | 7830 |
| chematic_legacy_etkdg | 265 | 1.0 | 6.8 | 9.7 | 31 |

#### RDKit

| Arm | n | p50 (ms) | p95 (ms) | p99 (ms) | max (ms) |
|---|---|---|---|---|---|
| rdkit_etkdgv3_raw | 265 | 11.0 | 78.8 | 140.7 | 332 |
| rdkit_etkdgv3_uff | 265 | 20.0 | 113.6 | 189.7 | 409 |
| rdkit_etkdgv3_mmff94 | 265 | 22.0 | 121.0 | 196.7 | 421 |
| rdkit_etkdgv3_best_of_n | 265 | 192.0 | 1080.4 | 2046.6 | 3047 |

## Cyclopentane RDKit crash — scoped ablation

**Classification: `nondefault_small_ring_torsion_only`**

12/60 trials crashed. Crashing configs (`useSmallRingTorsions`, `enforceChirality`): ['(True, False)', '(True, True)']. Crashing seeds: [4, 20260801]. Crashes under RDKit's own default config (`useSmallRingTorsions=False`): 0.

In plain terms: this crash requires the non-default `useSmallRingTorsions=True`, occurs during `EmbedMolecule` itself (before any force-field stage runs), and only reproduces for a subset of tested seeds -- **not** a general "RDKit crashes on cyclopentane" finding, and not reproducible under RDKit's own ETKDGv3 defaults in this ablation. Minimal repro: `scripts/pipeline_v2_vs_rdkit_cyclopentane_crash_ablation.py`.

## Force-field coverage (chematic MMFF94 arms)

- chematic_pipeline_v2_mmff94_with_uff_fallback: n=252, fallback_rate=40.9%, converged_rate=17.1%
- chematic_pipeline_v2_mmff94_strict: n=148, fallback_rate=0.0%, converged_rate=24.3%
- chematic_pipeline_v2_mmff94_strict_repair: n=131, fallback_rate=0.0%, converged_rate=26.0%
- chematic_pipeline_v2_mmff94_with_uff_fallback_repair: n=229, fallback_rate=42.8%, converged_rate=17.5%
- chematic_pipeline_v2_mmff94_strict_stretch_bend_gated: n=37, fallback_rate=0.0%, converged_rate=73.0%
- chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated: n=250, fallback_rate=85.2%, converged_rate=14.0%
- chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated: n=37, fallback_rate=0.0%, converged_rate=73.0%
- chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated: n=250, fallback_rate=85.2%, converged_rate=14.0%

## Stage funnel (per-arm denominator hierarchy)

Real `pipeline_v2` execution order (`crates/chematic-3d/src/pipeline_v2.rs` `PipelineStage` enum + its actual call sequence): embed (`DistanceGeometry`) -> torsion optimization -> **stereo verify/repair** -> force-field minimization -> final stereo verify -> final geometry validation. Stereo repair happens *before* force-field minimization, not after -- the columns below follow that real order, not an assumed embed-then-FF-then-stereo sequence. A row is counted under an `_attempted`/`_reached` column if its `failure_stage` is that stage or later (or it succeeded outright); under a `_succeeded`/`_verified` column only if `failure_stage` is strictly later than that stage (or it succeeded outright) -- a row that failed AT a stage reached it but did not pass it, so `ff_attempted` and `ff_succeeded` are genuinely different counts, not the same check twice. Never collapsed into a single success rate -- see `feedback_fallback_pooling_measurement_error`: `mmff94_strict` and `mmff94_with_uff_fallback` are reported as fully separate rows, never blended.

| Arm | attempted | embed_succeeded | stereo_repair_reached | ff_attempted | ff_succeeded | final_stereo_verified | final_validation_passed |
|---|---|---|---|---|---|---|---|
| chematic_pipeline_v2_no_ff | 265 | 254 | 254 | 254 | 254 | 254 | 254 |
| chematic_pipeline_v2_dreiding | 265 | 254 | 254 | 254 | 254 | 254 | 254 |
| chematic_pipeline_v2_uff_only | 265 | 254 | 254 | 254 | 250 | 250 | 250 |
| chematic_pipeline_v2_mmff94_strict | 265 | 254 | 254 | 254 | 148 | 148 | 148 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 254 | 254 | 254 | 252 | 252 | 252 |
| chematic_pipeline_v2_mmff94_strict_repair | 265 | 254 | 254 | 244 | 142 | 131 | 131 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 265 | 254 | 254 | 244 | 242 | 229 | 229 |
| chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 265 | 254 | 254 | 254 | 37 | 37 | 37 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 265 | 254 | 254 | 254 | 250 | 250 | 250 |
| chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 265 | 254 | 254 | 254 | 37 | 37 | 37 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 265 | 254 | 254 | 254 | 250 | 250 | 250 |
| chematic_legacy_etkdg | 265 | n/a | n/a | n/a | n/a | n/a | 265 |

Note: `chematic_legacy_etkdg` does not run through `pipeline_v2` at all (separate `generate_coords_etkdg` entry point, no `PipelineStage` tracking) -- its row reports `attempted`/`final_validation_passed` only; the intermediate columns are `n/a` rather than a fabricated 0 or a misleading 265 (a naive reuse of the success-implies-passed-every-stage rule above would have printed 265 for every column here, which would misrepresent a code path that never runs those stages at all).

## RepairAndVerify effectiveness (paired-arm comparison, Priority 1 new arms)

Each `StereoPolicy::RepairAndVerify` arm is a genuinely independent arm (not a config edit to the pre-existing `Ignore`-policy arm of the same `ForceFieldPolicy`), paired here per-molecule against its Ignore counterpart. `repair_time_delta` is `(repair-arm elapsed_ms - Ignore-arm elapsed_ms)` per molecule -- a paired-arm difference, not a directly instrumented repair-stage timer (pipeline_v2 does not currently expose one).

**Why `before-mismatch`/`repair attempted`/`repair succeeded` are identical between the two repair arms below**: verified directly (not assumed) -- `stereo_before_violations` matches per-molecule, 1:1, across `chematic_pipeline_v2_mmff94_strict_repair` and `chematic_pipeline_v2_mmff94_with_uff_fallback_repair` for all 265 molecules. This is structural, not a bug: stereo verify/repair runs BEFORE force-field minimization in `pipeline_v2`'s real execution order (see the stage-funnel note above), so both arms see the identical pre-FF geometry and make identical repair decisions -- the two `ForceFieldPolicy` values can only diverge afterward. `after-mismatch` DOES differ (11 vs. 13) because `final_stereo` is measured after FF minimization, where the two force fields' behavior can differ.

| Ignore arm | Repair arm | n compared | excluded (incomparable) | before-mismatch | repair attempted | repair succeeded | outcome unavailable | after-mismatch | geometry pairs | geometry degraded | time delta median (ms) | time delta p95 (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chematic_pipeline_v2_mmff94_strict | chematic_pipeline_v2_mmff94_strict_repair | 254 | 11 | 46 | 64 | 52 | 0 | 11 | 131 | 0 | -2 | 1147 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | chematic_pipeline_v2_mmff94_with_uff_fallback_repair | 254 | 11 | 46 | 64 | 52 | 0 | 13 | 229 | 0 | 0 | 646 |

_Note on reading `after-mismatch` next to the stereo-preservation table's "100% satisfaction" figure below: they are measured over different populations. The 100% satisfaction rate is computed only over the arm's *successful* rows (131/229); `after-mismatch` here is computed over all 254 comparable rows, including ones that failed after repair (e.g. in force-field minimization) -- so a non-zero after-mismatch count and a 100%-among-successes satisfaction rate are not in conflict, they answer different denominators. See the Stage funnel table above for the exact per-arm counts at each stage._

## Bonded-term coverage gate (Priority 2 / Stage 1B, issue #227)

Stretch-bend cross terms (Halgren MMFF.V eq. 4) were historically never gated by `Mmff94BondAngleStrict`/`Mmff94WithUffFallback` at all -- `Mmff94CoverageReport` had no field for them, and `stretch_bend_energy` (chematic-ff's `mmff94_minimizer`) silently contributes zero energy for an uncovered term instead of erroring. This PR adds `gate_mmff94_stretch_bend` (`PipelineV2Config`/`minimize_with_policy_gated`), an independent opt-in with the same shape as the pre-existing `gate_mmff94_torsion_oop` -- and 4 new benchmark arms exercising a real 3-stage comparison (legacy -> stretch-bend-gated -> complete-bonded-term-gated, for both `mmff94_strict` and `mmff94_with_uff_fallback`). **Review-driven correction**: an earlier version of this PR only added the stretch-bend-gated stage and mislabeled its result as "true complete-term coverage" -- it left `gate_mmff94_torsion_oop` at its default `false`, so torsion/OOP (1,121/0 missing instances respectively, measured below) were never actually gated. The complete-bonded-term stage below fixes this by gating stretch-bend AND torsion AND OOP together. Still not "complete MMFF94" -- vdW and partial-charge coverage are never gated by any arm in this benchmark.

### Missing-term sub-classification (fresh re-run, `mmff94_term_coverage_audit.rs`)

Per-term-instance classification across the 265-molecule corpus, superseding the stale pre-#236/#238/#239/#241 numbers in `docs/mmff94_coverage_gap_227_audit.md` (that document cited **~6,900** missing stretch-bend terms; the fresh count below is **substantially lower**, most likely because the same typing fixes that closed the bond/angle gap in Priority 1 also improved stretch-bend's shared angle-type classification path). `routing_bug_candidate` = this exact atom-type tuple has a table row at a *different* classification code than the one this molecule's context computed -- a candidate for an `angle_type_for`/`torsion_type_for`/`bond_type_for` classification bug, not necessarily a genuine table gap. `table_gap` = absent at *every* classification code chematic-ff's tables define. `Oop` is listed explicitly even at 0 -- omitting a measured-zero term kind would be indistinguishable from "not measured", which it is not.

| Term kind | total missing instances | routing_bug_candidate | table_gap |
|---|---|---|---|
| Bond | 84 | 79 (94.0%) | 5 (6.0%) |
| Angle | 374 | 277 (74.1%) | 97 (25.9%) |
| Torsion | 1121 | 1107 (98.8%) | 14 (1.2%) |
| Oop | 0 | 0 (n/a) | 0 (n/a) |
| StretchBend | 2107 | 427 (20.3%) | 1680 (79.7%) |

For Bond/Angle/Torsion/Oop, `table_gap` is not further sub-classified this round -- chematic-ff implements neither MMFF94 equivalence-class substitution nor empirical-rule (e.g. Badger's-rule bond) estimation at all for these term kinds (`Mmff94NumericTypeInfo.equivalence_levels` carries real MMFF94 equivalence data but has zero readers anywhere in the codebase, verified, not assumed) -- deferred, not fabricated.

**StretchBend's `table_gap` IS further sub-classified below** -- a diagnostic-only (never wired into production chematic-ff) port of RDKit's REAL stretch-bend resolution path, verified against the pinned RDKit commit (`scripts/mmff94_provenance/PROVENANCE.md`): `MMFFMolProperties::getMMFFStretchBendParams` tries the specific/generic MMFF-type lookup first (matches chematic's existing `mmff94_stbn` chain structurally), and on failure falls back to `MMFFDfsbCollection::getMMFFDfsbParams` -- a small (29-row), periodic-table-row-keyed default table (`defaultMMFFDfsb`), ported verbatim here. **Confirmed: no equivalence-class (`eqLevel`) step exists anywhere in RDKit's real stretch-bend path** -- `eqLevel` is used only by RDKit's angle/torsion/OOP fallback functions -- so this Dfsb port is RDKit's *complete* residual fallback story for stretch-bend, not a partial one.

| table_gap total | dfsb_default_resolvable | truly_unsupported (under RDKit's own algorithm) |
|---|---|---|
| 1680 | 1680 (100.0%) | 0 (0.0%) |

**100% of the 1680-instance StretchBend table_gap would be resolved by porting this 29-row Dfsb table into chematic-ff's production stretch-bend resolution** (verified with negative controls -- the port correctly returns `false` for out-of-table periodic-row combinations and for the table's one all-zero row, not a blanket `true`). This reframes the earlier "primarily a genuine table-completeness gap, next step unclear" conclusion: the largest missing-term bucket across all 5 term kinds is concretely, narrowly closable -- a small, well-defined, low-risk follow-up PR (port the table into production, not build new inference machinery), not resolved in *this* PR.

### Legacy -> stretch-bend -> complete-bonded-term (3-stage paired comparison)

Each stage's arm is a genuinely independent arm (never a config edit to a previous stage's arm), compared per-molecule against the immediately preceding stage. For `mmff94_strict` (pure gate, no fallback), widening the gate can only ever turn a prior success into a failure, never the reverse -- verified as a hard invariant at generation time (both stage transitions), not just a display column. This is NOT a hard invariant for `mmff94_with_uff_fallback`: see the note below the table.

| Earlier stage | Later stage | n compared | earlier success | later success | newly failing |
|---|---|---|---|---|---|
| chematic_pipeline_v2_mmff94_strict | chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | 265 | 148 | 37 | 111 |
| chematic_pipeline_v2_mmff94_strict_stretch_bend_gated | chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated | 265 | 37 | 37 | 0 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | 265 | 252 | 250 | 2 |
| chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated | chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated | 265 | 250 | 250 | 0 |

`chematic_pipeline_v2_mmff94_strict_stretch_bend_gated` newly-failing molecules (111): aspirin, atorvastatin_fragment, chloropropene_E, chloropropene_Z, cinnamic_acid_E, cinnamic_acid_Z, cubane, d_alanine, diphenhydramine, gly_ala_gly, ibuprofen, ibuprofen_S, l_alanine, l_serine, l_threonine, naproxen_S, paracetamol, penicillin_core, purine, testosterone, chembl_tier_b_0012, chembl_tier_b_0013, chembl_tier_b_0014, chembl_tier_b_0015, chembl_tier_b_0036, chembl_tier_b_0037, chembl_tier_b_0038, chembl_tier_b_0039, chembl_tier_b_0040, chembl_tier_b_0041, chembl_tier_b_0042, chembl_tier_b_0043, chembl_tier_b_0044, chembl_tier_b_0045, chembl_tier_b_0046, chembl_tier_b_0047, chembl_tier_b_0048, chembl_tier_b_0049, chembl_tier_b_0051, chembl_tier_b_0063, chembl_tier_b_0064, chembl_tier_b_0065, chembl_tier_b_0066, chembl_tier_b_0067, chembl_tier_b_0073, chembl_tier_b_0074, chembl_tier_b_0075, chembl_tier_b_0076, chembl_tier_b_0078, chembl_tier_b_0079, chembl_tier_b_0081, chembl_tier_b_0083, chembl_tier_b_0089, chembl_tier_b_0090, chembl_tier_b_0091, chembl_tier_b_0092, chembl_tier_b_0093, chembl_tier_b_0094, chembl_tier_b_0095, chembl_tier_b_0096, chembl_tier_b_0097, chembl_tier_b_0101, chembl_tier_b_0102, chembl_tier_b_0108, chembl_tier_b_0109, chembl_tier_b_0110, chembl_tier_b_0111, chembl_tier_b_0112, chembl_tier_b_0113, chembl_tier_b_0114, chembl_tier_b_0115, chembl_tier_b_0116, chembl_tier_b_0117, chembl_tier_b_0124, chembl_tier_b_0125, chembl_tier_b_0126, chembl_tier_b_0129, chembl_tier_b_0130, chembl_tier_b_0131, chembl_tier_b_0132, chembl_tier_b_0133, chembl_tier_b_0134, chembl_tier_b_0135, chembl_tier_b_0136, chembl_tier_b_0137, chembl_tier_b_0138, chembl_tier_b_0139, chembl_tier_b_0140, chembl_tier_b_0141, chembl_tier_b_0142, chembl_tier_b_0143, chembl_tier_b_0144, chembl_tier_b_0146, chembl_tier_b_0147, chembl_tier_b_0148, chembl_tier_b_0149, chembl_tier_b_0150, chembl_tier_b_0153, chembl_tier_b_0154, chembl_tier_b_0155, chembl_tier_b_0156, chembl_tier_b_0157, chembl_tier_b_0158, chembl_tier_b_0160, chembl_tier_b_0162, chembl_tier_b_0163, chembl_tier_b_0167, chembl_tier_b_0168, chembl_tier_b_0169, chembl_tier_b_0174, chembl_tier_b_0193

`chematic_pipeline_v2_mmff94_with_uff_fallback_stretch_bend_gated` newly-failing molecules (2): chembl_tier_b_0014, chembl_tier_b_0041

**Stretch-bend-gated -> complete-bonded-term is 0 newly-failing for both policies this round** (37->37, 250->250) -- every molecule that already survives the stretch-bend gate also has complete torsion+OOP coverage in this specific 265-molecule corpus. This is empirical, not structural (torsion has 1,121 missing instances measured above; they evidently concentrate on molecules that already fail the stretch-bend gate, in this corpus) -- a different, larger, or differently-composed corpus could show a non-zero delta at this stage. Practical effect for this run: the stretch-bend-gated and complete-bonded-term-gated success counts are numerically identical here, so the corrected, narrower name (`..._stretch_bend_gated`, not "true complete-term") only matters for what the number *means*, not for its value on this particular corpus.

**Run-to-run note on `mmff94_with_uff_fallback`'s legacy success count**: an earlier draft of this measurement (same code, same seed, different wall-clock conditions) found 251/265 legacy successes with 1 timeout-rescue exception; this run finds 252/265 with zero exceptions. Root-caused, not hand-waved: `chembl_tier_b_0166`'s legacy MMFF94 attempt sits right at the `total_timeout_ms=20000` boundary (18710ms in this run vs. 21378ms in the earlier one) -- `embed_seed` governs geometry/RNG determinism but NOT wall-clock scheduling, so a molecule this close to a real-time budget can legitimately flip between timeout and success across runs depending on machine load, independent of any code change. 252/265 also matches the legacy `mmff94_with_uff_fallback` count already on `main` (this PR's parent, #248) -- this run, not the earlier draft, is the representative one. Not a hidden source of non-determinism -- flagged here so a future re-run landing on 251 or 253 isn't mistaken for a regression.

Not adopted as the new default this round -- the legacy `gate_mmff94_torsion_oop=false, gate_mmff94_stretch_bend=false` arms remain the primary `mmff94_strict`/`mmff94_with_uff_fallback` numbers reported elsewhere in this document, per the plan's explicit staged-measurement instruction (measure each delta transparently before formalizing a new gate default, not silently change the existing success count).

## Ring-torsion FailClosed probe

1 row(s) -- demonstrates `RingTorsionApplicationPolicy::FailClosed`'s documented behavior. Not folded into any of the 12 main arms' coverage numbers (those use `DiagnosticOnly`).

## Reference geometry subset

Status: **insufficient_evidence**. No experimentally-determined reference conformers were available for this benchmark round. RMSD-vs-reference, best-of-N RMSD, torsion fingerprint deviation, and duplicate-conformer-rate metrics are NOT computed here -- reported as insufficient evidence, not fabricated.

## Known issues filed from this benchmark

- MMFF94 coverage gap (117/265 not successful under mmff94_strict, PR #236/#238/#239/#241 fixes already reflected in this run): https://github.com/kent-tokyo/chematic/issues/227

## Data integrity

- Unclassified rows: 0 (hard-gated at 0 by the report generator)
- chematic rows sha256: `59b609a4b51f2f5a...`
- RDKit rows sha256: `6325a07151a10dae...`
- All integrity gates (row-count, unclassified, atom-mapping, missing/mismatched coords, non-finite coords, common-scorer coverage, denominator self-consistency) passed at generation time -- see `run_integrity_gates` in this script.

## Conclusions

Classified per class/metric — no single overall win/loss score.

| Metric | Classification | Basis |
|---|---|---|
| Coverage — no_ff/dreiding/uff_only/mmff94_with_uff_fallback vs. RDKit | Roughly comparable | chematic 94.3%-95.8% success vs. RDKit 99.6% |
| Coverage — mmff94_strict | RDKit-favor (chematic gap, issue #227 filed) | 55.8% success, 104/265 unsupported |
| Common heavy-atom geometry — pipeline_v2 force-field arms | Chematic strength on soundness | 100% independently-sound across dreiding/uff_only/mmff94 arms, matching pipeline-internal judgment |
| Common heavy-atom geometry — legacy etkdg | Known gap, refined this round | 14/265 legacy outputs have coincident atoms under the stricter common scorer (not caught by the original Wave 1 ad-hoc check); the already-documented clash-rate gap stands |
| Stereo preservation (same judge, `Ignore`) | RDKit-favor | RDKit 100.0% satisfaction vs. chematic 62.9% under `StereoPolicy::Ignore` -- not chematic's best achievable number, see next row |
| Stereo preservation (same judge, `RepairAndVerify`, new this round) | Parity with RDKit among successes, coverage gap remains the real cost | mmff94_strict_repair 100.0%, mmff94_with_uff_fallback_repair 100.0% satisfaction among molecules that reached success under RepairAndVerify (both match RDKit's 100% on that subset) -- but RepairAndVerify also reduces the success *count* vs. the matching Ignore arm (fewer molecules reach final success at all when repair is required to pass); see the RepairAndVerify effectiveness section for the exact paired accounting |
| Bonded-term coverage gate, mmff94_strict -> mmff94_strict_stretch_bend_gated (new this round) | Real coverage gap surfaced, widening the gate is a real cost | 148 earlier-stage successes -> 37 under the later stage's gate (111 newly fail, 75.0% of earlier-stage successes) -- see the Bonded-term coverage gate section for the term-kind sub-classification and full molecule list |
| Bonded-term coverage gate, mmff94_strict_stretch_bend_gated -> mmff94_strict_complete_bonded_term_gated (new this round) | Real coverage gap surfaced, widening the gate is a real cost | 37 earlier-stage successes -> 37 under the later stage's gate (0 newly fail, 0.0% of earlier-stage successes) -- see the Bonded-term coverage gate section for the term-kind sub-classification and full molecule list |
| Bonded-term coverage gate, mmff94_with_uff_fallback -> mmff94_with_uff_fallback_stretch_bend_gated (new this round) | Real coverage gap surfaced, widening the gate is a real cost | 252 earlier-stage successes -> 250 under the later stage's gate (2 newly fail, 0.8% of earlier-stage successes) -- see the Bonded-term coverage gate section for the term-kind sub-classification and full molecule list |
| Bonded-term coverage gate, mmff94_with_uff_fallback_stretch_bend_gated -> mmff94_with_uff_fallback_complete_bonded_term_gated (new this round) | Real coverage gap surfaced, widening the gate is a real cost | 250 earlier-stage successes -> 250 under the later stage's gate (0 newly fail, 0.0% of earlier-stage successes) -- see the Bonded-term coverage gate section for the term-kind sub-classification and full molecule list |
| Force-field convergence rate | RDKit-favor, and an input to Priority 3 (Stage 1C) | chematic mmff94_with_uff_fallback 17.1% converged within 200 iterations, yet 252/265 of that arm's runs pass final validation regardless -- i.e. most successful outputs did NOT converge within 200 iterations and still passed geometry validation. Either `force_field_converged` is narrower than "produced a usable geometry" (an iteration-budget artifact, not necessarily a quality problem), or this is a real gap worth diagnosing -- Priority 3's MinimizationFailed root-causing (CatastrophicBondBlowup vs. ExcessiveResidualForce) is the next stage that should resolve which; corroborates open issues #185/#188 |
| Known crashes | RDKit has a narrowly-scoped one; chematic none found this round | cyclopentane crash classified `nondefault_small_ring_torsion_only` -- non-default config, seed-dependent, not RDKit's own default behavior |
| Unsupported chemistry | RDKit-favor | chematic mmff94_strict 117/265 unsupported (issue #227); RDKit's 4 arms show 0 unsupported_chemistry rows |
| Reference-geometry accuracy / torsion fingerprint / conformer diversity | Insufficient evidence | not measured this round, not fabricated |
| Overall "does chematic beat RDKit" | Not claimed | per this program's explicit rule -- findings are class/metric-specific |

