# V3000 flat-coordinate geometry allocation — 2026-10-15

The production optimization is rejected for adoption; only its experiment report and evidence are delivered. No toolkit-wide superiority, current-main or release claim follows from this run.

Decision: Reject for adoption: no original-profile parse or roundtrip interval establishes an improvement. Original parse medians are 0.996, 0.995 and 1.005 with intervals including 1; continuation parse is 0.984 [0.964, 1.007]. All outputs and tests match, but avoiding this allocation did not establish faster complete operations. Existing Open Babel chain wins are inherited and do not justify this source change.

The V3000 reader previously built a temporary `Vec<Point3>` for every geometry check, even when all z coordinates were near zero and the vector was immediately discarded. A small shared helper now resolves the original empty/flat prefix directly from the z channel; the V3000 reader materializes points only for the remaining cases. The shared classifier uses that same helper. The original strict `abs(z) < 1e-4` comparison, empty-input rank and non-finite fallback are unchanged. The z scan truncates at the xy channel length exactly as the previous zip did. Non-flat geometry, optional conformers, dimensional diagnostics and square-planar/wedge checks retain the old path. No public API, dependency or unsafe code is added.

The immediate baseline contains all seven preceding accepted candidates, including the latest aromaticity component reuse. It excludes the rejected direct-line experiment. Both native release/LTO executables were built freshly before and after this change; source, binaries, corpus and Open Babel libraries are hashed. All original 2,064 records retain byte-identical parser snapshots and rewrites between the two source arms, and all 64 physical-continuation variants match the original chains in both arms.

Scoped Rust validation: 844 passed, no failures or ignored tests. This includes the full MOL library and four integration suites for 2D E/Z, wedge stereo, 3D conformers and square-planar I/O. Runner semantic gates: 10 passed. Release all-target Clippy, formatting and refactor checks passed. The frozen original classifier and complete V3000 geometry function are compared against 1,911 numeric channel cases and all 10,000 existing molecules under 30 generated geometry/dimension conditions each. Exact conformer point bits, diagnostics, atom and bond fields match. The original recipe, license, minimization, oracle source hashes and normalized payload hash are recorded in the new flat-geometry fixture manifest. No new molecular records were added.

## Protocol

All fixed profiles are retained without filtering: 1,000 records each at offsets 0 and 5,000, and 64 linear chains of lengths 16 through 1,024 in steps of 16. This is an exposed development corpus on one macOS arm64 host, without a sealed holdout. Open Babel 3.2.1 native library/plugins and RDKit 2026.09.1 observer are pinned. Each original profile has 21 rotating three-arm fresh-process blocks and five full-batch repeats per mode. All builds and validation precede timing; profiles and diagnostics run sequentially. Only this run owned, completed build target was removed after retaining both measured binaries.

Ratios are comparator time / candidate time; above 1 means faster. Brackets are paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. Strict wins require semantic agreement, CI lower > 1 and all 21 blocks faster. Every record, block, mismatch, failure and initial/adaptive negative is retained. The writer is unchanged; writer timing movement is not attributed to a new writer algorithm.

Timers exclude startup, input loading/framing and plugin discovery. The native harness prepares a parsed batch and performs an untimed warm write before timing. Parse retains the complete diagnostic batch, prepared writing creates output strings, and roundtrip parses and writes fresh records. RSS is whole-process high-water usage including preparation and plugins; it is not per-operation allocation evidence. The removed temporary allocation follows from the implementation, not an allocation-counter measurement.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 164.196 | 164.027 | 0.996 [0.994, 1.000]; 7/21 | unranked (semantic mismatch) |
| 0 | write | 21.514 | 21.397 | 1.006 [0.993, 1.020]; 12/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 186.445 | 186.759 | 0.999 [0.991, 1.003]; 7/21 | unranked (semantic mismatch) |
| 5000 | parse | 47.417 | 47.612 | 0.995 [0.990, 1.003]; 7/21 | unranked (semantic mismatch) |
| 5000 | write | 7.700 | 7.735 | 0.995 [0.994, 1.000]; 7/21 | unranked (semantic mismatch) |
| 5000 | roundtrip | 55.221 | 55.326 | 1.001 [0.993, 1.004]; 11/21 | unranked (semantic mismatch) |
| chains | parse | 34.549 | 34.403 | 1.005 [0.999, 1.014]; 13/21 | 119.640 [118.411, 119.980]; 21/21 **strict win** |
| chains | write | 15.367 | 15.374 | 0.989 [0.979, 1.011]; 10/21 | 6.190 [6.118, 6.263]; 21/21 **strict win** |
| chains | roundtrip | 49.954 | 50.141 | 0.999 [0.990, 1.004]; 10/21 | 84.216 [83.630, 84.928]; 21/21 **strict win** |

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
| parse | 54.252 | 55.371 | 0.984 [0.964, 1.007]; 8/21 |
| write | 15.888 | 16.074 | 0.997 [0.971, 1.008]; 10/21 |
| roundtrip | 70.482 | 71.830 | 0.978 [0.969, 1.005]; 7/21 |

## Adaptive regression checks

Any initial baseline-relative interval wholly below 1 selects that profile for one separate 21-block source-only follow-up in all three modes. All initial results remain above; follow-ups are not pooled or substituted, and new follow-up negatives do not cause recursive retries.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {}
}
```

## Reproduction and status

[The evidence packet](2026-10-15-openbabel-v3000-flat-geometry-artifacts/README.md) contains both measured source executables, every initial/follow-up result and semantic output, full build/test logs, source manifests, fixture recipe/provenance, the exact patch and a complete build/input archive. Every 10053 candidate archive hash and 10050 reconstructed baseline hash was verified. The patch also reconstructs the changed files in an independent tiny Git repository. The initial Cargo invocation used the wrong working directory; its preserved error log precedes the successful baseline build and no timing used that failure.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared refs were not refreshed; hashes identify version 1.1.4 source candidates, not current main. This run has no hosted CI, merge or publication. Previous source-candidate negatives in unrelated Python or native lanes are not resolved by an incremental comparison. The earlier rejected direct-line candidate remains archived and excluded.
