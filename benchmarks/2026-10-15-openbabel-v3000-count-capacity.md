# V3000 bounded count-based capacity — 2026-10-15

The production optimization is rejected for adoption; only its experiment report and evidence are delivered. No toolkit-wide superiority, current-main or release claim follows from this run.

Decision: Reject exact-capacity reservations for adoption. Initial parser ratios improve 1.038x (all21), 1.064x (20/21) and 1.039x (20/21), but every original profile has an initial prepared-writer CI wholly below 1. Separate one-per-profile follow-ups retain clear writer regressions at NCI 0.967x [0.944,0.986] and chains 0.959x [0.930,0.986]. Original-profile roundtrip gains do not remove this prepared-operation cost; continuation/other follow-ups include wide or crossing intervals and are retained. All output and 847 native/10 runner gates pass. Only the report/evidence is delivered. A separate architecture revision will round reservation capacities toward the old Vec growth sizes while retaining the logical-row cap; no earlier observation will be pooled or replaced.

The V3000 reader reserves known atom/bond capacities in its builder, atom/bond source-ID maps, XY coordinates and Z-coordinate storage. Reservation occurs once at BEGIN ATOM after the final COUNTS row, so repeated COUNTS does not allocate or discard stores repeatedly. Every requested capacity is bounded by the actual number of logical V30 rows; a huge declared count in a short input does not trigger a huge allocation. Counts remain validation values, never a limit on parsed records. Parsing order, state transitions, typed-error precedence, coordinates, stereo/geometry and metadata are unchanged. No public API, dependency, tolerance or unsafe code is added.

The immediate baseline contains all nine preceding accepted candidates, including the six-decimal positive-zero parser and descriptor normalization storage reuse. All rejected direct-line, flat-geometry and two-spelling zero-z trials remain excluded. Both native release/LTO executables were built freshly before and after this change; source, binaries, corpus and Open Babel libraries are hashed. All original 2,064 records retain byte-identical parser snapshots and rewrites between the two source arms, and all 64 physical-continuation variants match the original chains in both arms.

Scoped Rust validation: 847 passed, no failures or ignored tests. This includes the full MOL library and four integration suites for 2D E/Z, wedge stereo, 3D conformers and square-planar I/O. Runner semantic gates: 10 passed. Release all-target Clippy, formatting and refactor checks passed. The frozen original COUNTS transition is compared on 27,000 count-lexeme, row-budget and source-line combinations, with a prior COUNTS row to check repeated declarations and partial mutations. Incomplete rows and block markers are also checked. A separate test repeats maximum-platform-integer declarations 16 times without allocation, then checks exact bounded capacities at BEGIN ATOM. Full typed-error Debug representations, expected counts and state discriminants agree. The original synthetic recipe, license, minimization, pinned source oracle and normalized payload hash are recorded in the new count-capacity fixture manifest. The first native-test build stopped on a private constructor signature change; its failure log and pre-repair source are retained. Constructor compatibility was restored and all relevant gates rerun successfully before timing. No new molecular records were added.

## Protocol

All fixed profiles are retained without filtering: 1,000 records each at offsets 0 and 5,000, and 64 linear chains of lengths 16 through 1,024 in steps of 16. This is an exposed development corpus on one macOS arm64 host, without a sealed holdout. Open Babel 3.2.1 native library/plugins and RDKit 2026.09.1 observer are pinned. Each original profile has 21 rotating three-arm fresh-process blocks and five full-batch repeats per mode. All builds and validation precede timing; profiles and diagnostics run sequentially. Only this run owned, completed build target was removed after retaining both measured binaries.

Ratios are comparator time / candidate time; above 1 means faster. Brackets are paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. Strict wins require semantic agreement, CI lower > 1 and all 21 blocks faster. Every record, block, mismatch, failure and initial/adaptive negative is retained. The writer is unchanged; writer timing movement is not attributed to a new writer algorithm.

Timers exclude startup, input loading/framing and plugin discovery. The native harness prepares a parsed batch and performs an untimed warm write before timing. Parse retains the complete diagnostic batch, prepared writing creates output strings, and roundtrip parses and writes fresh records. RSS is whole-process high-water usage including preparation and plugins; it is not per-operation allocation evidence. The avoided geometric vector growth follows from the implementation and capacity assertions; no operation-level allocation count is measured.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 166.948 | 162.522 | 1.038 [1.024, 1.054]; 21/21 **strict win** | unranked (semantic mismatch) |
| 0 | write | 21.723 | 22.068 | 0.977 [0.974, 0.999]; 6/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 188.919 | 184.721 | 1.021 [1.019, 1.039]; 20/21 | unranked (semantic mismatch) |
| 5000 | parse | 56.075 | 52.370 | 1.064 [1.051, 1.072]; 20/21 | unranked (semantic mismatch) |
| 5000 | write | 9.104 | 9.385 | 0.970 [0.943, 0.992]; 5/21 | unranked (semantic mismatch) |
| 5000 | roundtrip | 65.548 | 61.593 | 1.056 [1.047, 1.062]; 20/21 | unranked (semantic mismatch) |
| chains | parse | 38.932 | 37.598 | 1.039 [1.029, 1.041]; 20/21 | 127.722 [126.049, 128.214]; 21/21 **strict win** |
| chains | write | 17.974 | 18.506 | 0.970 [0.956, 0.983]; 3/21 | 5.980 [5.879, 6.173]; 21/21 **strict win** |
| chains | roundtrip | 57.334 | 56.983 | 1.011 [1.001, 1.021]; 17/21 | 87.049 [86.063, 88.836]; 21/21 **strict win** |

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
| parse | 60.407 | 58.430 | 1.012 [0.983, 1.031]; 13/21 |
| write | 17.409 | 17.784 | 0.961 [0.873, 0.993]; 5/21 |
| roundtrip | 83.576 | 77.889 | 0.986 [0.916, 1.009]; 7/21 |

## Adaptive regression checks

Any initial baseline-relative interval wholly below 1 selects that profile for one separate 21-block source-only follow-up in all three modes. All initial results remain above; follow-ups are not pooled or substituted, and new follow-up negatives do not cause recursive retries.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {
    "0": {
      "negative_initial_lanes": {
        "write": {
          "median": 0.977,
          "ci95": [
            0.974,
            0.999
          ],
          "winning_blocks": 6,
          "blocks": 21
        }
      },
      "corpus": "/private/tmp/chematic-openbabel-v3000-20261014/corpus-0.molblocks",
      "records": 1000,
      "gate": {
        "source_snapshot_sha256": {
          "baseline": "e2fce68d8c33616317f46397cf79ee877ef856460a35b84a622ddff1d8eced74",
          "candidate": "e2fce68d8c33616317f46397cf79ee877ef856460a35b84a622ddff1d8eced74"
        },
        "byte_identical_rewrites": true,
        "source_parse_snapshots_equal": true
      }
    },
    "5000": {
      "negative_initial_lanes": {
        "write": {
          "median": 0.97,
          "ci95": [
            0.943,
            0.992
          ],
          "winning_blocks": 5,
          "blocks": 21
        }
      },
      "corpus": "/private/tmp/chematic-openbabel-v3000-20261014/corpus-5000.molblocks",
      "records": 1000,
      "gate": {
        "source_snapshot_sha256": {
          "baseline": "7fa169eebb0d112936fd4e135d08702eb1150faa6791111f52792c51d504d67d",
          "candidate": "7fa169eebb0d112936fd4e135d08702eb1150faa6791111f52792c51d504d67d"
        },
        "byte_identical_rewrites": true,
        "source_parse_snapshots_equal": true
      }
    },
    "chains": {
      "negative_initial_lanes": {
        "write": {
          "median": 0.97,
          "ci95": [
            0.956,
            0.983
          ],
          "winning_blocks": 3,
          "blocks": 21
        }
      },
      "corpus": "/private/tmp/chematic-openbabel-v3000-20261014/corpus-chains.molblocks",
      "records": 64,
      "gate": {
        "source_snapshot_sha256": {
          "baseline": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550",
          "candidate": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550"
        },
        "byte_identical_rewrites": true,
        "source_parse_snapshots_equal": true
      }
    },
    "continued": {
      "negative_initial_lanes": {
        "write": {
          "median": 0.961,
          "ci95": [
            0.873,
            0.993
          ],
          "winning_blocks": 5,
          "blocks": 21
        }
      },
      "corpus": "/private/tmp/chematic-openbabel-count-capacity-20261015/corpus-continued-chains.molblocks",
      "records": 64,
      "gate": {
        "all_original_and_continued_source_snapshots_and_rewrites_byte_identical": true,
        "sha256": {
          "continued-baseline.snapshot": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550",
          "original-candidate.molblocks": "b67b1fb9cdb4c2e8977a34dd309a13ebc035050ba7e41473114c62f2c7d5ca95",
          "continued-candidate.snapshot": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550",
          "original-baseline.molblocks": "b67b1fb9cdb4c2e8977a34dd309a13ebc035050ba7e41473114c62f2c7d5ca95",
          "continued-candidate.molblocks": "b67b1fb9cdb4c2e8977a34dd309a13ebc035050ba7e41473114c62f2c7d5ca95",
          "original-baseline.snapshot": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550",
          "continued-baseline.molblocks": "b67b1fb9cdb4c2e8977a34dd309a13ebc035050ba7e41473114c62f2c7d5ca95",
          "original-candidate.snapshot": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550"
        },
        "records": 64
      }
    }
  }
}
```

### focused21-chains

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---:|
| parse | 1.039 [1.005, 1.080]; 15/21 |
| write | 0.959 [0.930, 0.986]; 6/21 |
| roundtrip | 0.996 [0.967, 1.016]; 8/21 |

### focused21-0

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---:|
| parse | 1.004 [0.821, 1.328]; 11/21 |
| write | 0.962 [0.680, 1.173]; 10/21 |
| roundtrip | 0.960 [0.898, 1.173]; 10/21 |

### focused21-5000

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---:|
| parse | 1.049 [1.030, 1.066]; 16/21 |
| write | 0.967 [0.944, 0.986]; 6/21 |
| roundtrip | 1.032 [1.021, 1.049]; 17/21 |

### focused21-continued

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---:|
| parse | 1.036 [0.906, 1.180]; 11/21 |
| write | 0.975 [0.931, 1.088]; 8/21 |
| roundtrip | 1.058 [0.916, 1.149]; 12/21 |

## Reproduction and status

[The evidence packet](2026-10-15-openbabel-v3000-count-capacity-artifacts/README.md) contains both measured source executables, every initial/follow-up result and semantic output, full build/test logs, source manifests, fixture recipe/provenance, the exact patch and a complete build/input archive. Every 10057 candidate archive hash and 10054 reconstructed baseline hash was verified. The patch also reconstructs the changed files in an independent tiny Git repository. All previous rejected experiments remain separate, archived and excluded. The fixed inputs and benchmark policy are unchanged; no observation from an earlier run is pooled with this run.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared refs were not refreshed; hashes identify version 1.1.4 source candidates, not current main. This run has no hosted CI, merge or publication. Previous source-candidate negatives in unrelated Python or native lanes are not resolved by an incremental comparison. The earlier rejected direct-line candidate remains archived and excluded.
