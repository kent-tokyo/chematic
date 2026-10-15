# V3000 positive-zero z parsing — 2026-10-15

The production optimization is adopted in the attached source candidate. No toolkit-wide superiority, current-main or release claim follows from this run.

Decision: Adopt the six-decimal revision as a small unmerged source optimization. All original 2,064 snapshots/rewrites and 64 continuation variants match exactly; 844 Rust tests and 10 runner gates pass. First-profile parse is 1.010x [1.002, 1.014], 16/21 faster; NCI roundtrip 1.015x [1.008, 1.024], 15/21; controlled-chain parse 1.040x [1.030, 1.049] and roundtrip 1.031x [1.023, 1.049], both 20/21. Continuation parse/roundtrip intervals are also positive. These source-relative lanes do not meet the all-21-block strict criterion. No initial interval is wholly below 1, so no adaptive follow-up is selected. Diverse Open Babel lanes remain unranked for semantic differences; only controlled chains pass every comparator gate. The two rejected trials remain separate and are not pooled.

The V3000 coordinate reader now handles the exact finite positive-zero z lexemes `0`, `0.0000` and `0.000000` directly, avoiding the general float parser for those common 2D values. Every other spelling, including signed zero, still uses the original parser and finite-value check. x and y are evaluated first in the original order, preserving the first typed error, its axis, input text and source line. No public API, dependency, tolerance or unsafe code is added.

The immediate baseline contains all seven preceding accepted candidates, including the latest aromaticity component reuse. It excludes the rejected direct-line experiment. Both native release/LTO executables were built freshly before and after this change; source, binaries, corpus and Open Babel libraries are hashed. All original 2,064 records retain byte-identical parser snapshots and rewrites between the two source arms, and all 64 physical-continuation variants match the original chains in both arms.

Scoped Rust validation: 844 passed, no failures or ignored tests. This includes the full MOL library and four integration suites for 2D E/Z, wedge stereo, 3D conformers and square-planar I/O. Runner semantic gates: 10 passed. Release all-target Clippy, formatting and refactor checks passed. The frozen original coordinate parser is compared on 62,500 Cartesian lexeme/source-line cases and 10,000 deterministic f64 fallback cases. Exact tuple bits and complete typed-error representations match, including negative-zero bits and x/y/z failure precedence. The original numeric recipe, license, minimization, oracle source hash and normalized payload hash are recorded in the new zero-z fixture manifest. No new molecular records were added.

## Protocol

All fixed profiles are retained without filtering: 1,000 records each at offsets 0 and 5,000, and 64 linear chains of lengths 16 through 1,024 in steps of 16. This is an exposed development corpus on one macOS arm64 host, without a sealed holdout. Open Babel 3.2.1 native library/plugins and RDKit 2026.09.1 observer are pinned. Each original profile has 21 rotating three-arm fresh-process blocks and five full-batch repeats per mode. All builds and validation precede timing; profiles and diagnostics run sequentially. Only this run owned, completed build target was removed after retaining both measured binaries.

Ratios are comparator time / candidate time; above 1 means faster. Brackets are paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. Strict wins require semantic agreement, CI lower > 1 and all 21 blocks faster. Every record, block, mismatch, failure and initial/adaptive negative is retained. The writer is unchanged; writer timing movement is not attributed to a new writer algorithm.

Timers exclude startup, input loading/framing and plugin discovery. The native harness prepares a parsed batch and performs an untimed warm write before timing. Parse retains the complete diagnostic batch, prepared writing creates output strings, and roundtrip parses and writes fresh records. RSS is whole-process high-water usage including preparation and plugins; it is not per-operation allocation evidence. The general-parser bypass follows from the implementation; no per-operation allocation counter is measured.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 167.816 | 165.853 | 1.010 [1.002, 1.014]; 16/21 | unranked (semantic mismatch) |
| 0 | write | 21.907 | 21.767 | 1.007 [0.993, 1.018]; 13/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 190.339 | 188.220 | 1.006 [1.000, 1.015]; 15/21 | unranked (semantic mismatch) |
| 5000 | parse | 51.372 | 52.654 | 1.008 [0.987, 1.024]; 13/21 | unranked (semantic mismatch) |
| 5000 | write | 8.224 | 8.442 | 0.985 [0.974, 1.019]; 8/21 | unranked (semantic mismatch) |
| 5000 | roundtrip | 60.479 | 61.422 | 1.015 [1.008, 1.024]; 15/21 | unranked (semantic mismatch) |
| chains | parse | 38.854 | 33.683 | 1.040 [1.030, 1.049]; 20/21 | 123.444 [122.911, 124.777]; 21/21 **strict win** |
| chains | write | 17.131 | 15.811 | 1.009 [0.991, 1.023]; 11/21 | 6.301 [6.101, 6.405]; 21/21 **strict win** |
| chains | roundtrip | 56.578 | 49.879 | 1.031 [1.023, 1.049]; 20/21 | 87.510 [86.351, 97.128]; 21/21 **strict win** |

## Semantic exclusions

Open Babel disagreements remain unranked. No stereo marker is removed and no row is dropped to obtain a speed result. Diverse inputs include graph mismatches, explicit unknown-double-stereo loss and rejected rewrites; those lanes cannot support equivalent-output speed wins. Controlled-chain agreement does not establish broad chemical correctness or superiority.

| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |
|---|---:|---:|---:|---:|---:|
| 0 | 3 | 0 | 153 | 1 | 1 |
| 5000 | 137 | 0 | 0 | 137 | 137 |
| chains | 0 | 0 | 0 | 0 | 0 |

## Separate continuation diagnostic

The predeclared source-only diagnostic splits every numeric V30 atom and bond row of all controlled chains after its id into continued physical lines. All original and continued snapshots and rewrites agree between both source arms. It uses 21 alternating two-arm fresh-process blocks, all modes and five repeats; it has no Open Babel ranking or RSS measurement and never replaces an original observation.

| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |
|---|---:|---:|---|
| parse | 53.506 | 51.616 | 1.026 [1.021, 1.029]; 18/21 |
| write | 15.457 | 15.635 | 0.994 [0.986, 1.004]; 8/21 |
| roundtrip | 68.441 | 67.010 | 1.021 [1.015, 1.025]; 19/21 |

## Adaptive regression checks

Any initial baseline-relative interval wholly below 1 selects that profile for one separate 21-block source-only follow-up in all three modes. All initial results remain above; follow-ups are not pooled or substituted, and new follow-up negatives do not cause recursive retries.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {}
}
```

## Reproduction and status

[The evidence packet](2026-10-15-openbabel-v3000-zero-z6-artifacts/README.md) contains both measured source executables, every initial/follow-up result and semantic output, full build/test logs, source manifests, fixture recipe/provenance, the exact patch and a complete build/input archive. Every 10053 candidate archive hash and 10050 reconstructed baseline hash was verified. The patch also reconstructs the changed files in an independent tiny Git repository. The preceding [flat-geometry allocation experiment](2026-10-15-openbabel-v3000-flat-geometry.md) was rejected because none of its original-profile parse intervals established a gain. Its source, outcomes and logs remain archived separately. Both its point-allocation changes and shared geometry helper are excluded from this trial. The subsequent [two-spelling zero-z trial](2026-10-15-openbabel-v3000-zero-z.md) was also rejected: its fast spellings did not occur in any fixed input. Its changes are excluded from the seven-candidate baseline; the current candidate extends the spelling set. An audit found all z lexemes are `0.000000`; this separate revision adds that six-decimal spelling before its fixed measurement runs. All three experiments remain separate and are not pooled.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared refs were not refreshed; hashes identify version 1.1.4 source candidates, not current main. This run has no hosted CI, merge or publication. Previous source-candidate negatives in unrelated Python or native lanes are not resolved by an incremental comparison. The earlier rejected direct-line candidate remains archived and excluded.
