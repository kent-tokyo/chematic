# pipeline v2 vs RDKit ETKDGv3 — Wave 1 independent 3D benchmark

Measurement-only. No pipeline v2 or force-field algorithm code was changed to produce these numbers. Historical numbers (e.g. a previously-cited 25.9% geometrically-valid rate) are NOT reused here -- everything below was regenerated fresh against this repo's current `main` in this session.

## Corpus

- Tier A (curated stress): 65 molecules, sha256 `8d4f3f3a70b8ae00...`
- Tier B (fixed drug-like, ChEMBL-derived): 200 molecules, sha256 `1a1698f5444c1b6e...`
- Total: 265 molecules

## Atom mapping

- Checked (both engines produced a heavy-atom element sequence): 265
- Verified matching (per-index element-symbol equality): 265
- Unavailable/mismatched: 0

## Coverage — chematic (6 arms)

| Arm | n | success | typed_failure | unsupported_chemistry | timeout | internal_error |
|---|---|---|---|---|---|---|
| chematic_pipeline_v2_no_ff | 265 | 254 | 11 | 0 | 0 | 0 |
| chematic_pipeline_v2_dreiding | 265 | 254 | 11 | 0 | 0 | 0 |
| chematic_pipeline_v2_uff_only | 265 | 250 | 15 | 0 | 0 | 0 |
| chematic_pipeline_v2_mmff94_strict | 265 | 38 | 11 | 216 | 0 | 0 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 250 | 15 | 0 | 0 | 0 |
| chematic_legacy_etkdg | 265 | 265 | 0 | 0 | 0 | 0 |

## Coverage — RDKit ETKDGv3 (4 arms)

| Arm | n | success | oracle_failure | unsupported_chemistry | internal_error |
|---|---|---|---|---|---|
| rdkit_etkdgv3_raw | 265 | 264 | 0 | 0 | 1 |
| rdkit_etkdgv3_uff | 265 | 264 | 0 | 0 | 1 |
| rdkit_etkdgv3_mmff94 | 265 | 264 | 0 | 0 | 1 |
| rdkit_etkdgv3_best_of_n | 265 | 264 | 0 | 0 | 1 |

## Geometry validity (chematic)

Independent per-engine measurement -- not a cross-engine RMSD (no reference geometry available, see below).

| Arm | n success | all finite | sound | mean bond-viol >15% | mean bond-viol >50% | mean clashes | molecules w/ clash |
|---|---|---|---|---|---|---|---|
| chematic_pipeline_v2_no_ff | 254 | 100.0% | 100.0% | 12.6% | 0.0% | 0.03 | 3 |
| chematic_pipeline_v2_dreiding | 254 | 100.0% | 100.0% | 0.9% | 0.0% | 0.00 | 0 |
| chematic_pipeline_v2_uff_only | 250 | 100.0% | 100.0% | 0.0% | 0.0% | 0.00 | 0 |
| chematic_pipeline_v2_mmff94_strict | 38 | 100.0% | 100.0% | 0.0% | 0.0% | 0.00 | 0 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | 250 | 100.0% | 100.0% | 0.0% | 0.0% | 0.00 | 0 |
| chematic_legacy_etkdg | 265 | 100.0% | 100.0% | 49.8% | 15.5% | 26.93 | 236 |

## Stereo (chematic, molecules with declared stereo only)

| Arm | n w/ declared stereo | mean declared | mean satisfied | mean violations | mean unevaluable | total repaired | total repair failed |
|---|---|---|---|---|---|---|---|
| chematic_pipeline_v2_no_ff | 83 | 1.76 | 0.99 | 0.77 | 0.00 | 0 | 0 |
| chematic_pipeline_v2_dreiding | 83 | 1.76 | 1.07 | 0.69 | 0.00 | 0 | 0 |
| chematic_pipeline_v2_uff_only | 80 | 1.75 | 1.10 | 0.65 | 0.00 | 0 | 0 |
| chematic_pipeline_v2_mmff94_strict | 19 | 1.89 | 0.89 | 1.00 | 0.00 | 0 | 0 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | 80 | 1.75 | 1.10 | 0.65 | 0.00 | 0 | 0 |
| chematic_legacy_etkdg | 0 | n/a | n/a | n/a | n/a | n/a | n/a |

## Force-field coverage (chematic MMFF94 arms)

- chematic_pipeline_v2_mmff94_with_uff_fallback: n=250, fallback_rate=84.8%, converged_rate=14.0%
- chematic_pipeline_v2_mmff94_strict: n=38, fallback_rate=0.0%, converged_rate=65.8%

## Performance

_In-process wall-clock timing per (molecule, arm) call within a single long-running process -- NOT separate-process-isolated. p50/p95/p99/max reported per arm; process-level variance (repeated whole-process runs) was NOT measured this round._

### chematic

| Arm | n | p50 (ms) | p95 (ms) | p99 (ms) | max (ms) |
|---|---|---|---|---|---|
| chematic_pipeline_v2_no_ff | 265 | 9.0 | 55.8 | 93.3 | 119 |
| chematic_pipeline_v2_dreiding | 265 | 306.0 | 2063.2 | 2701.2 | 3040 |
| chematic_pipeline_v2_uff_only | 265 | 523.0 | 3436.8 | 4754.0 | 7937 |
| chematic_pipeline_v2_mmff94_strict | 265 | 18.0 | 157.4 | 1196.7 | 4140 |
| chematic_pipeline_v2_mmff94_with_uff_fallback | 265 | 591.0 | 3283.4 | 4548.7 | 7146 |
| chematic_legacy_etkdg | 265 | 2.0 | 9.0 | 14.0 | 17 |

### RDKit

| Arm | n | p50 (ms) | p95 (ms) | p99 (ms) | max (ms) |
|---|---|---|---|---|---|
| rdkit_etkdgv3_raw | 265 | 25.0 | 210.6 | 371.6 | 803 |
| rdkit_etkdgv3_uff | 265 | 48.0 | 292.0 | 495.4 | 977 |
| rdkit_etkdgv3_mmff94 | 265 | 54.0 | 323.4 | 559.1 | 1024 |
| rdkit_etkdgv3_best_of_n | 265 | 489.0 | 2730.0 | 5059.8 | 7327 |

## Ring-torsion FailClosed probe

1 row(s) -- demonstrates `RingTorsionApplicationPolicy::FailClosed`'s documented behavior on the dedicated `known_fail_closed_case` fixture. Not folded into the 6 main arms' coverage numbers (those use `DiagnosticOnly`, see the dump executable's own comment for why).

## Reference geometry subset

Status: **insufficient_evidence**. No experimentally-determined reference conformers were available for this benchmark round. RMSD-vs-reference, best-of-N RMSD, torsion fingerprint deviation, and duplicate-conformer-rate metrics are NOT computed here -- reported as insufficient evidence, not fabricated against a synthetic or absent reference.

## Data integrity

- Unclassified rows: 0 (must be 0; see aggregate JSON if not)
- chematic rows sha256: `585889615aae9002...`
- RDKit rows sha256: `ba24a7c64df1c350...`

## Notable findings (not folded into the tables above)

- **RDKit itself crashes on plain cyclopentane under this benchmark's config.** All
  4 RDKit arms show exactly 1 `internal_error` — the same molecule
  (`cyclopentane`, Tier A) on every arm: `RuntimeError: Invariant Violation —
  bad direction in linearSearch` (`Code/Numerics/Optimizer/BFGSOpt.h:224`),
  RDKit 2026.03.3. Independently reproduced directly (not just observed in
  the dump) with `useSmallRingTorsions=True` + the seed this benchmark uses.
  Not a chematic finding to claim credit for — reported because a failure is
  a failure regardless of which engine it belongs to, and because it
  explains 100% of RDKit's `internal_error` bucket across all 4 arms.
- **`Mmff94BondAngleStrict` has severe real-world parameter coverage gaps.**
  216/265 (81.5%) of the corpus lands in `unsupported_chemistry` under that
  arm — including plain benzene (missing aromatic-ring angle/torsion MMFF94
  parameters). `Mmff94WithUffFallback` recovers nearly all of that
  (250/265 success) — but 84.8% of its successes are actually silent-to-the-
  caller UFF fallbacks (`force_field_fallback: true`), not real MMFF94 runs.
- **Force-field convergence (200-iteration cap) is a real bottleneck.**
  `Mmff94WithUffFallback`'s reported `force_field_converged` rate is only
  14.0% across its 250 successes — most of that arm's runs are the UFF
  fallback path, and most of those don't converge within 200 iterations
  (soundness/finiteness is still fine — `sound: true` in the geometry table
  — but the minimizer isn't reaching its own convergence criterion).
  Consistent with this repo's existing open issues #185/#188 (UFF minimizer
  blow-up/non-monotonic behavior) — not new evidence against those issues,
  but independent corroboration from a different measurement path.
- **The legacy `etkdg` entry point "succeeds" on 100% of the corpus but
  produces geometrically poor structures most of the time.** 236/265 (89.1%)
  of its outputs have at least one non-bonded clash (mean 26.93 clashes per
  molecule), and 15.5% of its bonds deviate from ideal length by >50% on
  average — the legacy path is infallible (always returns *a* geometry) but
  that is not the same as returning a *good* one. Contrast: every pipeline
  v2 arm shows 0% mean bond-violation >50% and near-zero clash rates.
- **Chematic's force-field arms are meaningfully slower than RDKit's
  equivalents on wall-clock time**, in this in-process, non-isolated
  measurement: `uff_only` p50 523ms vs. RDKit `etkdgv3_uff` p50 48ms (~11x);
  `mmff94_with_uff_fallback` p50 591ms vs. RDKit `etkdgv3_mmff94` p50 54ms
  (~11x). Notably, RDKit's `best_of_n` (10 full embed+UFF-optimize cycles)
  has a p50 of 489ms — comparable to chematic doing a *single* uff_only
  embed. This is a real, substantial performance gap, not a rounding
  difference.
- **`no_ff`/`dreiding` arms (95.8% success) and force-field arms
  (94.3%/94.3%) show non-trivial, non-force-field-related typed-failure
  rates too** (11-15 molecules) — not investigated further in this
  measurement-only round; see the aggregate JSON's `coverage_by_class` for
  the exact failing molecules/categories.
- **Atom mapping held for 100% of the corpus** (265/265) — every molecule's
  heavy-atom element sequence matched exactly between chematic's parse and
  RDKit's post-`AddHs` heavy atoms, confirming the "same SMILES string, same
  atom order in both engines" assumption this benchmark's cross-referencing
  depends on, rather than assuming it.

## Conclusions

Classified per class/metric — no single overall win/loss score.

| Metric | Classification | Basis |
|---|---|---|
| Coverage — `no_ff`/`dreiding`/`uff_only`/`mmff94_with_uff_fallback` vs. RDKit's 4 arms | **Roughly comparable** | chematic 94.3-95.8% success vs. RDKit 99.6% (264/265) — RDKit slightly ahead, within a plausible margin given differing failure taxonomies (chematic's failures are typed pipeline stages; RDKit's are embed/FF exceptions) |
| Coverage — `mmff94_strict` | **RDKit-favor** (chematic gap) | 14.3% success — real, measured MMFF94 bond/angle parameter coverage gap, not a benchmark artifact (confirmed via `MissingParameters` cause on plain benzene) |
| Geometry validity — force-field-refined arms (`dreiding`/`uff_only`/`mmff94_with_uff_fallback`) | **Roughly comparable / chematic strength on soundness** | 100% all-finite, 100% sound, ≤0.9% mean bond-violation>15%, 0% mean bond-violation>50%, near-zero clashes — no reference geometry exists to compare absolute accuracy against RDKit's own output, but chematic's internal soundness gate is consistently met |
| Geometry validity — legacy `etkdg` | **Chematic gap (known, not new)** | 89.1% of outputs have ≥1 clash; this is the *documented, already-diagnosed* legacy path (`docs/etkdg_3d_gap_rfc.md`), not pipeline v2 — the gap this whole 3D Breakthrough Program exists to close |
| Stereo preservation (declared-stereo subset) | **Roughly comparable, evidence-limited** | ~55-63% mean satisfaction rate across arms with `StereoPolicy::Ignore` (non-gating); no RDKit-side stereo-preservation metric was computed this round for direct comparison — flagged as a gap for a follow-up round, not fabricated here |
| Force-field convergence rate | **RDKit-favor** | chematic's 200-iteration-capped MMFF94-with-fallback shows only 14.0% converged; not directly compared to RDKit's own convergence rate this round (not measured on the RDKit side) — reported as a chematic-side finding, corroborating existing issues #185/#188 |
| Performance — force-field arms | **RDKit-favor, substantial** | ~11x slower wall-clock p50 on `uff_only`/`mmff94_with_uff_fallback` vs. RDKit's equivalents; RDKit's 10-conformer best-of-N is about as fast as chematic's single attempt |
| Performance — `no_ff`/legacy | **Roughly comparable** | chematic `no_ff` p50 9ms, legacy p50 2ms vs. RDKit raw p50 25ms — chematic faster here, but `no_ff`/legacy skip force-field work RDKit's `raw` arm doesn't |
| Reference-geometry accuracy (RMSD vs. experimental conformers) | **Insufficient evidence** | no reference conformers available this round — not fabricated |
| Torsion fingerprint deviation / conformer diversity / duplicate-conformer rate | **Insufficient evidence** | not computed this round (requires the reference-geometry subset above, or a dedicated ensemble-diversity study) |
| Overall "does chematic beat RDKit" | **Not claimed** | per this program's explicit rule — findings are class/metric-specific, several favor RDKit clearly (MMFF94 coverage, force-field convergence, raw speed), one is a known pre-existing gap (legacy path), several are genuinely comparable |

