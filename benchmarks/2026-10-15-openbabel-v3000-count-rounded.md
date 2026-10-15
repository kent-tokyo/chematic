# V3000 bounded rounded capacity — 2026-10-15

The production optimization is rejected for adoption; only its experiment report and evidence are delivered. No toolkit-wide superiority, current-main or release claim follows from this run.

Decision: Reject rounded capacities alone: original profile parser/roundtrip gains are retained, but prepared-writer regressions remain in at least one prescribed separate follow-up. Rounded capacities did not eliminate the prepared-operation cost. All output and 848 native/10 runner gates pass. Archive all initial/continuation/follow-up results separately. The next implementation adds a direct writer optimization for the common zero atom-map field; it will be tested and measured as a separate full candidate, never pooled with this run.

The V3000 reader reserves known atom/bond capacities in its builder, atom/bond source-ID maps, XY coordinates and Z-coordinate storage. Reservation occurs once at BEGIN ATOM after the final COUNTS row, so repeated COUNTS does not allocate or discard stores repeatedly. Requested capacities approach the original geometric Vec growth size and remain bounded by the actual number of logical V30 rows; a huge declared count in a short input does not trigger a huge allocation. Counts remain validation values, never a limit on parsed records. Parsing order, state transitions, typed-error precedence, coordinates, stereo/geometry and metadata are unchanged. No public API, dependency, tolerance or unsafe code is added.

The immediate baseline contains all nine preceding accepted candidates, including the six-decimal positive-zero parser and descriptor normalization storage reuse. All rejected direct-line, flat-geometry and two-spelling zero-z trials remain excluded. The baseline native release/LTO executable was built before both capacity trials; its exact binary and all baseline source hashes are verified and reused. The candidate is built freshly after this revision; source, binaries, corpus and Open Babel libraries are hashed. All original 2,064 records retain byte-identical parser snapshots and rewrites between the two source arms, and all 64 physical-continuation variants match the original chains in both arms.

Scoped Rust validation: 848 passed, no failures or ignored tests. This includes the full MOL library and four integration suites for 2D E/Z, wedge stereo, 3D conformers and square-planar I/O. Runner semantic gates: 10 passed. Release all-target Clippy, formatting and refactor checks passed. The frozen original COUNTS transition is compared on 27,000 count-lexeme, row-budget and source-line combinations, with a prior COUNTS row to check repeated declarations and partial mutations. Incomplete rows and block markers are also checked. A separate test repeats maximum-platform-integer declarations 16 times without allocation, then checks exact bounded capacities at BEGIN ATOM. Full typed-error Debug representations, expected counts and state discriminants agree. The original synthetic recipe, license, minimization, pinned source oracle and normalized payload hash are recorded in the new count-capacity fixture manifest. The exact-capacity trial is archived separately and rejected for reproduced prepared-writer regressions. This revision rounds reservations toward the original Vec growth capacities, still bounded by actual logical rows, and reruns every validation and fixed timing profile. A new test compares capacity hints with original Vec growth for counts 0 through 1,024. The exact same verified nine-candidate baseline executable is reused, with its build log and origin hash; new timing observations are not pooled. No new molecular records were added.

## Protocol

All fixed profiles are retained without filtering: 1,000 records each at offsets 0 and 5,000, and 64 linear chains of lengths 16 through 1,024 in steps of 16. This is an exposed development corpus on one macOS arm64 host, without a sealed holdout. Open Babel 3.2.1 native library/plugins and RDKit 2026.09.1 observer are pinned. Each original profile has 21 rotating three-arm fresh-process blocks and five full-batch repeats per mode. All builds and validation precede timing; profiles and diagnostics run sequentially. Only this run owned, completed build target was removed after retaining both measured binaries.

Ratios are comparator time / candidate time; above 1 means faster. Brackets are paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. Strict wins require semantic agreement, CI lower > 1 and all 21 blocks faster. Every record, block, mismatch, failure and initial/adaptive negative is retained. The writer is unchanged; writer timing movement is not attributed to a new writer algorithm.

Timers exclude startup, input loading/framing and plugin discovery. The native harness prepares a parsed batch and performs an untimed warm write before timing. Parse retains the complete diagnostic batch, prepared writing creates output strings, and roundtrip parses and writes fresh records. RSS is whole-process high-water usage including preparation and plugins; it is not per-operation allocation evidence. The avoided geometric vector growth follows from the implementation and capacity assertions; no operation-level allocation count is measured.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 162.490 | 159.117 | 1.021 [1.020, 1.029]; 20/21 | unranked (semantic mismatch) |
| 0 | write | 21.211 | 21.638 | 0.976 [0.971, 0.983]; 2/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 184.553 | 181.334 | 1.017 [1.009, 1.021]; 18/21 | unranked (semantic mismatch) |
| 5000 | parse | 47.087 | 44.312 | 1.065 [1.062, 1.072]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | write | 7.779 | 7.943 | 0.980 [0.972, 0.990]; 2/21 | unranked (semantic mismatch) |
| 5000 | roundtrip | 54.813 | 52.154 | 1.051 [1.047, 1.052]; 21/21 **strict win** | unranked (semantic mismatch) |
| chains | parse | 33.287 | 32.447 | 1.025 [1.019, 1.034]; 20/21 | 126.763 [126.269, 127.012]; 21/21 **strict win** |
| chains | write | 15.500 | 15.710 | 0.978 [0.963, 0.992]; 4/21 | 6.059 [5.997, 6.078]; 21/21 **strict win** |
| chains | roundtrip | 48.669 | 48.220 | 1.010 [1.003, 1.017]; 18/21 | 87.321 [86.932, 87.984]; 21/21 **strict win** |

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
| parse | 52.277 | 51.738 | 1.015 [0.996, 1.041]; 14/21 |
| write | 15.863 | 15.925 | 0.976 [0.955, 1.014]; 7/21 |
| roundtrip | 68.334 | 67.173 | 1.015 [0.997, 1.023]; 14/21 |

## Adaptive regression checks

Any initial baseline-relative interval wholly below 1 selects that profile for one separate 21-block source-only follow-up in all three modes. All initial results remain above; follow-ups are not pooled or substituted, and new follow-up negatives do not cause recursive retries.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {
    "0": {
      "negative_initial_lanes": {
        "write": {
          "median": 0.976,
          "ci95": [
            0.971,
            0.983
          ],
          "winning_blocks": 2,
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
          "median": 0.98,
          "ci95": [
            0.972,
            0.99
          ],
          "winning_blocks": 2,
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
          "median": 0.978,
          "ci95": [
            0.963,
            0.992
          ],
          "winning_blocks": 4,
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
    }
  }
}
```

### focused21-chains

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---:|
| parse | 1.003 [0.945, 1.056]; 11/21 |
| write | 0.899 [0.864, 0.967]; 6/21 |
| roundtrip | 0.961 [0.919, 0.993]; 5/21 |

### focused21-0

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---:|
| parse | 0.990 [0.848, 1.170]; 10/21 |
| write | 1.066 [0.863, 1.134]; 12/21 |
| roundtrip | 1.002 [0.926, 1.126]; 11/21 |

### focused21-5000

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---:|
| parse | 1.052 [1.042, 1.061]; 19/21 |
| write | 0.971 [0.953, 0.991]; 3/21 |
| roundtrip | 1.037 [1.033, 1.041]; 20/21 |

## Reproduction and status

[The evidence packet](2026-10-15-openbabel-v3000-count-rounded-artifacts/README.md) contains both measured source executables, every initial/follow-up result and semantic output, full build/test logs, source manifests, fixture recipe/provenance, the exact patch and a complete build/input archive. Every 10057 candidate archive hash and 10054 reconstructed baseline hash was verified. The patch also reconstructs the changed files in an independent tiny Git repository. All previous rejected experiments remain separate, archived and excluded. The fixed inputs and benchmark policy are unchanged; no observation from an earlier run is pooled with this run.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared refs were not refreshed; hashes identify version 1.1.4 source candidates, not current main. This run has no hosted CI, merge or publication. Previous source-candidate negatives in unrelated Python or native lanes are not resolved by an incremental comparison. The earlier rejected direct-line candidate remains archived and excluded.
