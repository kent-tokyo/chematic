# pipeline v2 vs RDKit ETKDGv3 — Wave 1 independent 3D benchmark

Measurement-only. No pipeline v2 or force-field algorithm code was changed to produce these numbers. Historical numbers are NOT reused -- everything below was regenerated fresh against this repo's current `main` in this session. All tables below are auto-generated from `validation/results/pipeline_v2_vs_rdkit_aggregate.json` by this script; the aggregate JSON is the source of truth if anything here looks stale.

## Corpus

- Tier A (curated stress): 65 molecules, sha256 `8d4f3f3a70b8ae00...`
- Tier B (fixed drug-like, ChEMBL-derived): 200 molecules, sha256 `1a1698f5444c1b6e...`
- Total: 265 molecules

## Atom mapping

- Checked: 265, verified matching: 265, unavailable/mismatched: 0

## Coverage and usable geometry (explicit denominators)

`usable_coverage` = independently-sound successes / total inputs for that arm -- the fraction of the *whole corpus* that arm turns into a geometry this benchmark's own independent scorer (not the pipeline's internal judgment) certifies sound. `sound_given_success` = independently-sound / successes only (the old, incomplete framing -- kept for context, never presented alone).

| Engine | Arm | total | success | indep. sound | sound_given_success | usable_coverage | typed_failure | unsupported | timeout | internal_error |
|---|---|---|---|---|---|---|---|---|---|---|
| chematic | chematic_pipeline_v2_no_ff | 265 | 254 | 254 | 100.0% | 95.8% | 11 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_dreiding | 265 | 254 | 254 | 100.0% | 95.8% | 11 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_uff_only | 265 | 250 | 250 | 100.0% | 94.3% | 15 | 0 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_strict | 265 | 38 | 38 | 100.0% | 14.3% | 11 | 216 | 0 | 0 |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 250 | 250 | 100.0% | 94.3% | 15 | 0 | 0 | 0 |
| chematic | chematic_legacy_etkdg | 265 | 265 | 251 | 94.7% | 94.7% | 0 | 0 | 0 | 0 |
| rdkit | rdkit_etkdgv3_raw | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |
| rdkit | rdkit_etkdgv3_uff | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |
| rdkit | rdkit_etkdgv3_mmff94 | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |
| rdkit | rdkit_etkdgv3_best_of_n | 265 | 264 | 264 | 100.0% | 99.6% | 0 | 0 | 0 | 1 |

**mmff94_strict, spelled out per the fix request:** 38/38 successful outputs are independently sound, but only 38/265 of the *total corpus* ends up as a usable geometry under this arm -- the rest is the 216-molecule MMFF94 parameter coverage gap (issue #227), not a geometry-quality problem.

## Common heavy-atom geometry quality (same independent scorer, both engines)

Applied identically to chematic's and RDKit's already-saved heavy-atom coordinates (`crates/chematic-3d/examples/pipeline_v2_vs_rdkit_common_scorer.rs`) -- ideal bond length from `Element::covalent_radius()`, never chematic-3d's own `pub(crate)` thresholds. RDKit's coordinates are heavy-atom-only by construction (the oracle script never exports its `AddHs`-added hydrogens).

| Engine | Arm | n scored | all finite | mean bond>15% | mean bond>50% | molecules w/ clash | molecules w/ coincident atoms | independently sound |
|---|---|---|---|---|---|---|---|---|
| chematic | chematic_pipeline_v2_no_ff | 254 | 100.0% | 2.8% | 0.0% | 3 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_dreiding | 254 | 100.0% | 3.0% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_uff_only | 250 | 100.0% | 0.5% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_strict | 38 | 100.0% | 0.7% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback | 250 | 100.0% | 0.6% | 0.0% | 0 | 0 | 100.0% |
| chematic | chematic_legacy_etkdg | 265 | 100.0% | 49.8% | 15.5% | 229 | 14 | 94.7% |
| rdkit | rdkit_etkdgv3_raw | 264 | 100.0% | 0.6% | 0.0% | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_uff | 264 | 100.0% | 0.0% | 0.0% | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_mmff94 | 264 | 100.0% | 0.9% | 0.0% | 0 | 0 | 100.0% |
| rdkit | rdkit_etkdgv3_best_of_n | 264 | 100.0% | 0.0% | 0.0% | 0 | 0 | 100.0% |

Note (correction vs. the original Wave 1 report): the legacy `etkdg` arm was previously reported as 100% sound. This common scorer additionally checks for exactly-coincident atom pairs (distance < 1e-3 Å), which the original ad-hoc legacy scorer did not -- 14/265 legacy outputs have ≥1 coincident atom pair and are NOT independently sound under this stricter, shared check. All 5 pipeline_v2 arms remain 100% independently sound (matching their own internal `final_validation.sound`).

## Stereo preservation (same judge -- chematic's own `verify_stereo` -- applied to both engines)

**Methodology, read before the numbers**: chematic's arms below were benchmarked with `StereoPolicy::Ignore` (deliberate Wave 1 choice, to keep coverage/geometry metrics free of stereo-driven failures). `Ignore` never repairs a violated stereocenter -- so these numbers reflect raw distance-geometry-embedding output, NOT chematic's best achievable stereo correctness (`StereoPolicy::RepairAndVerify`, not exercised this round). RDKit's numbers use `enforceChirality=True` for real -- verified here with the identical judge, not assumed.

| Engine | Arm | molecules w/ declared stereo | declared | satisfied | violated | unevaluable | satisfaction rate |
|---|---|---|---|---|---|---|---|
| chematic | chematic_pipeline_v2_no_ff | 83 | 146 | 82 | 64 | 0 | 56.2% |
| chematic | chematic_pipeline_v2_dreiding | 83 | 146 | 89 | 57 | 0 | 61.0% |
| chematic | chematic_pipeline_v2_uff_only | 80 | 140 | 88 | 52 | 0 | 62.9% |
| chematic | chematic_pipeline_v2_mmff94_strict | 19 | 36 | 17 | 19 | 0 | 47.2% |
| chematic | chematic_pipeline_v2_mmff94_with_uff_fallback | 80 | 140 | 88 | 52 | 0 | 62.9% |
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

### Process-level (primary comparison — separate OS process per run, sequential, 5 runs each)

_Whole-corpus (265 molecules), separate-process wall-clock per run, sequential (never concurrent with another run or a build). Includes process startup (Rust binary startup / Python+RDKit import) -- not steady-state-only._

| Engine | runs | median total (s) | min (s) | max (s) | stdev (s) | coeff. of variation |
|---|---|---|---|---|---|---|
| chematic | 5 | 312.8 | 303.7 | 615.3 | 136.15 | 0.366 |
| rdkit | 5 | 142.9 | 124.5 | 188.9 | 29.39 | 0.192 |

Whole-corpus median: chematic is ~2.2x slower than RDKit -- **substantially smaller** than the ~11x seen on the force-field-heavy arms alone (see in-process table below). This whole-corpus figure blends all 6 chematic arms (including the very fast `no_ff`/`legacy` arms) with all 4 RDKit arms; it is not in conflict with the per-arm figure, it answers a different question ("run the whole benchmark once" vs. "run this one force-field arm"). chematic's first run (615.3s) is a likely system-contention outlier relative to the other 4 (~304-320s, tight cluster) -- reported as-measured, not excluded, but flagged rather than silently averaged in as if typical; machine load average was already elevated (~6 on a 10-core machine) before this measurement began, from other concurrent activity on the same machine.

### In-process per-(molecule, arm) timing (secondary)

_In-process wall-clock per (molecule, arm) call within a single long-running process -- NOT process-isolated. Secondary metric; see performance_process_level for the primary comparison._

#### chematic

| Arm | n | p50 (ms) | p95 (ms) | p99 (ms) | max (ms) |
|---|---|---|---|---|---|
| chematic_pipeline_v2_no_ff | 265 | 9.0 | 55.8 | 93.3 | 119 |
| chematic_pipeline_v2_dreiding | 265 | 306.0 | 2063.2 | 2701.2 | 3040 |
| chematic_pipeline_v2_uff_only | 265 | 523.0 | 3436.8 | 4754.0 | 7937 |
| chematic_pipeline_v2_mmff94_strict | 265 | 18.0 | 157.4 | 1196.7 | 4140 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 591.0 | 3283.4 | 4548.7 | 7146 |
| chematic_legacy_etkdg | 265 | 2.0 | 9.0 | 14.0 | 17 |

#### RDKit

| Arm | n | p50 (ms) | p95 (ms) | p99 (ms) | max (ms) |
|---|---|---|---|---|---|
| rdkit_etkdgv3_raw | 265 | 25.0 | 210.6 | 371.6 | 803 |
| rdkit_etkdgv3_uff | 265 | 48.0 | 292.0 | 495.4 | 977 |
| rdkit_etkdgv3_mmff94 | 265 | 54.0 | 323.4 | 559.1 | 1024 |
| rdkit_etkdgv3_best_of_n | 265 | 489.0 | 2730.0 | 5059.8 | 7327 |

## Cyclopentane RDKit crash — scoped ablation

**Classification: `nondefault_small_ring_torsion_only`**

12/60 trials crashed. Crashing configs (`useSmallRingTorsions`, `enforceChirality`): ['(True, False)', '(True, True)']. Crashing seeds: [4, 20260801]. Crashes under RDKit's own default config (`useSmallRingTorsions=False`): 0.

In plain terms: this crash requires the non-default `useSmallRingTorsions=True`, occurs during `EmbedMolecule` itself (before any force-field stage runs), and only reproduces for a subset of tested seeds -- **not** a general "RDKit crashes on cyclopentane" finding, and not reproducible under RDKit's own ETKDGv3 defaults in this ablation. Minimal repro: `scripts/pipeline_v2_vs_rdkit_cyclopentane_crash_ablation.py`.

## Force-field coverage (chematic MMFF94 arms)

- chematic_pipeline_v2_mmff94_with_uff_fallback: n=250, fallback_rate=84.8%, converged_rate=14.0%
- chematic_pipeline_v2_mmff94_strict: n=38, fallback_rate=0.0%, converged_rate=65.8%

## Ring-torsion FailClosed probe

1 row(s) -- demonstrates `RingTorsionApplicationPolicy::FailClosed`'s documented behavior. Not folded into the 6 main arms' coverage numbers (those use `DiagnosticOnly`).

## Reference geometry subset

Status: **insufficient_evidence**. No experimentally-determined reference conformers were available for this benchmark round. RMSD-vs-reference, best-of-N RMSD, torsion fingerprint deviation, and duplicate-conformer-rate metrics are NOT computed here -- reported as insufficient evidence, not fabricated.

## Known issues filed from this benchmark

- MMFF94 coverage gap (216/265 unsupported, incl. plain benzene): https://github.com/kent-tokyo/chematic/issues/227

## Data integrity

- Unclassified rows: 0 (hard-gated at 0 by the report generator)
- chematic rows sha256: `585889615aae9002...`
- RDKit rows sha256: `ba24a7c64df1c350...`
- All integrity gates (row-count, unclassified, atom-mapping, missing/mismatched coords, non-finite coords, common-scorer coverage, denominator self-consistency) passed at generation time -- see `run_integrity_gates` in this script.

## Conclusions

Classified per class/metric — no single overall win/loss score.

| Metric | Classification | Basis |
|---|---|---|
| Coverage — no_ff/dreiding/uff_only/mmff94_with_uff_fallback vs. RDKit | Roughly comparable | chematic 94.3%-95.8% success vs. RDKit 99.6% |
| Coverage — mmff94_strict | RDKit-favor (chematic gap, issue #227 filed) | 14.3% success, 216/265 unsupported |
| Common heavy-atom geometry — pipeline_v2 force-field arms | Chematic strength on soundness | 100% independently-sound across dreiding/uff_only/mmff94 arms, matching pipeline-internal judgment |
| Common heavy-atom geometry — legacy etkdg | Known gap, refined this round | 14/265 legacy outputs have coincident atoms under the stricter common scorer (not caught by the original Wave 1 ad-hoc check); the already-documented clash-rate gap stands |
| Stereo preservation (same judge) | RDKit-favor, methodology caveat applies | RDKit 100.0% satisfaction vs. chematic 62.9% under `StereoPolicy::Ignore` (no repair attempted this round -- not chematic's best achievable number) |
| Force-field convergence rate | RDKit-favor | chematic mmff94_with_uff_fallback 14.0% converged within 200 iterations; corroborates open issues #185/#188 |
| Performance (process-level, whole corpus) | RDKit-favor | median 312.8s (chematic) vs. 142.9s (RDKit) for the full 265-molecule x arms run |
| Known crashes | RDKit has a narrowly-scoped one; chematic none found this round | cyclopentane crash classified `nondefault_small_ring_torsion_only` -- non-default config, seed-dependent, not RDKit's own default behavior |
| Unsupported chemistry | RDKit-favor | chematic mmff94_strict 216/265 unsupported (issue #227); RDKit's 4 arms show 0 unsupported_chemistry rows |
| Reference-geometry accuracy / torsion fingerprint / conformer diversity | Insufficient evidence | not measured this round, not fabricated |
| Overall "does chematic beat RDKit" | Not claimed | per this program's explicit rule -- findings are class/metric-specific |

