# Borrowed V3000 logical rows — 2026-10-16

Local adoption gate: **FAIL**. Complete Open Babel-wide superiority remains unproved. Native source measurements on one shared macOS arm64 host are separate from published artifacts.

Both primary parse intervals improve (1.020 [1.018, 1.023] and 1.041 [1.032, 1.049]), but the fixed retention gate fails: SDF writing has medians 0.986 and 0.983, and continued V3000 parsing has median 0.989, each below the 0.99 floor. The SDF writer algorithm is unchanged; these measurements do not isolate the cause. This candidate must remain unpromoted even if correctness CI succeeds.

## Implementation and preservation

Baseline: clean source main `2251a9f5d2e706087efaf55cb0f8ab1ff20b8016`, source version 1.1.5. Both Rust arms are freshly built from this base and the one-file production delta; earlier unpromoted sparse-ID and descriptor trials are excluded. The baseline includes requested-only writer losses and the dative-carbonyl E/Z correction.

V3000 reading previously built a complete vector of logical-line wrappers even when all payloads were already borrowed from the raw input. The candidate detects that no exact V30 payload ends in a continuation hyphen, counts the same logical rows, and feeds their borrowed payloads directly to the unchanged parser. Any continuation takes the original collection/joining path. Physical-line validation, error precedence, Unicode tokenization, exact allocation hints, chemical interpretation, full diagnostics, public signatures and dependencies are unchanged. No unsafe code. Physical lines and the full input remain resident; this is not a streaming-file API or memory bound.

A frozen complete original reader from the named base compares full returned graphs, atom/bond/tag/R-group/stereo-neighbor and group-membership data, coordinate and conformer float bits, metadata, every diagnostic family, and rewrites. It checks all 2,064 pinned V3000 records twice, long/short/empty/Unicode/continued rows, exact malformed-input error variants and fields, unusual headers, huge counts, allocation capacities, 3D stereo and enhanced collections. The original reader function tokens are independently checked after restoring only its name/visibility. Existing rich/query/opaque/3D tests remain active.

Full MOL checks: 1,271 passed, 0 failed, 0 ignored across 55 suites. Ten runner semantic gates, release all-target MOL Clippy, formatting and refactor budget pass. Native correctness uses opt3/LTO=false/16 codegen units, separate from the full-LTO/one-codegen-unit timing binaries. Complete preflight snapshots and rewrite bytes match across all 4,128 profile inputs. These inputs overlap in chemical structures.

## Frozen protocol

Two unfiltered 1,000-record SDF slices and two unfiltered 1,000-record V3000 slices at offsets 0 and 5,000, plus 64 controlled chains of 16–1,024 atoms and their continuation-only variants. Every atom row in the latter is split after its x coordinate; all seven declared observer fields match the original 64 chains. The two chain profiles represent the same chemistry with different physical formatting. All inputs are exposed development data, not a sealed holdout. Existing corpus provenance/licenses and third-party notices apply.

Twenty-one rotating three-arm fresh-process blocks per profile and three native batch repetitions per lane. Startup, input loading and initial plugin discovery are outside timers. SDF framing is inside parsing; V3000 delimiter splitting is outside. Parse retains the full diagnostic batch. Prepared write follows an untimed warm write; roundtrip parses fresh inputs and writes inside the timer. Rust SDF reuses its String; V3000 and C++ use allocating native return strings. Whole-process RSS includes loading, preparation, warmup, allocator context and plugins; it is not operation-only or streaming memory. External host activity is uncontrolled and load averages are retained.

Rust 1.97.0, release opt3/full LTO/one codegen unit/one build job, incremental=false, strip=none. Open Babel 3.2.1 C++ -O3 harnesses are reused by exact binary/library hashes; no historical timings are reused. RDKit 2026.09.1 is an independent observer, not a timed competitor. The 281 package-local native-file entries are pinned before/after each profile; aliases can repeat contents. System/transitive libraries outside the named package directories and native data files are outside this seal.

Ratios are reference time / candidate time. Paired median intervals use the unchanged deterministic 10,000-resample bootstrap and three-decimal reporting. Strict Open Babel wins require semantic prerequisites, interval lower>1 and 21/21 faster blocks. Internal canonical graph equality gates parsing; all seven declared rewrite fields gate write/roundtrip. Rejected/mismatching records remain in their original positions and unranked without filtering or repair.

The local plan was fixed before timing: exact source snapshots, rewrites and full reports/errors; V3000 parse interval lower>1 on both diverse slices; all 18 source-relative medians >=0.99; disclose every negative interval. No adaptive timing reruns, pooled epochs, substituted samples or relaxed gates. An unchanged lane moving does not establish its cause.

## Complete matrix

Times are median milliseconds per complete batch, native elapsed divided by three repetitions. Ratios use paired observations.

| Profile | Mode | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 12.681 | 12.523 | 1.013 [1.006, 1.015]; 17/21 | unranked: semantic prerequisite fails |
| sdf-0 | write | 3.929 | 3.990 | 0.986 [0.976, 0.989]; 3/21 | 6.090 [6.009, 6.248]; 21/21 **strict win** |
| sdf-0 | roundtrip | 16.802 | 16.737 | 1.004 [0.997, 1.010]; 12/21 | 7.170 [7.146, 7.215]; 21/21 **strict win** |
| sdf-5000 | parse | 4.919 | 4.827 | 1.012 [1.002, 1.023]; 15/21 | 8.994 [8.940, 9.052]; 21/21 **strict win** |
| sdf-5000 | write | 1.737 | 1.742 | 0.983 [0.970, 0.994]; 5/21 | 8.298 [8.147, 8.507]; 21/21 **strict win** |
| sdf-5000 | roundtrip | 6.678 | 6.668 | 0.993 [0.988, 1.015]; 7/21 | 8.747 [8.695, 8.923]; 21/21 **strict win** |
| v3000-0 | parse | 13.071 | 12.804 | 1.020 [1.018, 1.023]; 18/21 | unranked: semantic prerequisite fails |
| v3000-0 | write | 2.710 | 2.722 | 0.995 [0.990, 0.996]; 4/21 | unranked: semantic prerequisite fails |
| v3000-0 | roundtrip | 15.969 | 15.709 | 1.017 [1.010, 1.021]; 16/21 | unranked: semantic prerequisite fails |
| v3000-5000 | parse | 5.128 | 4.947 | 1.041 [1.032, 1.049]; 21/21 **strict win** | unranked: semantic prerequisite fails |
| v3000-5000 | write | 0.695 | 0.697 | 0.997 [0.992, 1.008]; 9/21 | unranked: semantic prerequisite fails |
| v3000-5000 | roundtrip | 5.833 | 5.673 | 1.033 [1.031, 1.039]; 20/21 | unranked: semantic prerequisite fails |
| v3000-chains | parse | 5.649 | 5.549 | 1.017 [1.009, 1.029]; 19/21 | 146.420 [145.448, 147.928]; 21/21 **strict win** |
| v3000-chains | write | 1.226 | 1.234 | 0.997 [0.990, 1.004]; 10/21 | 15.594 [15.424, 15.724]; 21/21 **strict win** |
| v3000-chains | roundtrip | 6.870 | 6.777 | 1.011 [1.006, 1.023]; 18/21 | 123.139 [122.222, 124.095]; 21/21 **strict win** |
| v3000-continued | parse | 7.736 | 7.809 | 0.989 [0.986, 0.998]; 6/21 | 105.363 [105.026, 106.091]; 21/21 **strict win** |
| v3000-continued | write | 1.218 | 1.215 | 1.001 [0.997, 1.007]; 12/21 | unranked: semantic prerequisite fails |
| v3000-continued | roundtrip | 9.025 | 9.038 | 1.000 [0.992, 1.003]; 10/21 | unranked: semantic prerequisite fails |

## Full decision and every negative interval

```json
{
  "local_adoption_gate_passed": false,
  "all_source_snapshots_and_rewrites_exact": true,
  "primary": {
    "v3000-0": {
      "speedup": {
        "median": 1.02,
        "ci95": [
          1.018,
          1.023
        ],
        "winning_blocks": 18,
        "blocks": 21
      },
      "semantic_prerequisite": true,
      "strict_speed_win": false
    },
    "v3000-5000": {
      "speedup": {
        "median": 1.041,
        "ci95": [
          1.032,
          1.049
        ],
        "winning_blocks": 21,
        "blocks": 21
      },
      "semantic_prerequisite": true,
      "strict_speed_win": true
    }
  },
  "primary_ci_pass": true,
  "all_18_source_median_floor": 0.99,
  "below_floor": {
    "sdf-0": {
      "write": {
        "median": 0.986,
        "ci95": [
          0.976,
          0.989
        ],
        "winning_blocks": 3,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "write": {
        "median": 0.983,
        "ci95": [
          0.97,
          0.994
        ],
        "winning_blocks": 5,
        "blocks": 21
      }
    },
    "v3000-0": {},
    "v3000-5000": {},
    "v3000-chains": {},
    "v3000-continued": {
      "parse": {
        "median": 0.989,
        "ci95": [
          0.986,
          0.998
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    }
  },
  "unranked_source_lanes": [],
  "all_negative_source_intervals": {
    "sdf-0": {
      "write": {
        "median": 0.986,
        "ci95": [
          0.976,
          0.989
        ],
        "winning_blocks": 3,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "write": {
        "median": 0.983,
        "ci95": [
          0.97,
          0.994
        ],
        "winning_blocks": 5,
        "blocks": 21
      }
    },
    "v3000-0": {
      "write": {
        "median": 0.995,
        "ci95": [
          0.99,
          0.996
        ],
        "winning_blocks": 4,
        "blocks": 21
      }
    },
    "v3000-5000": {},
    "v3000-chains": {},
    "v3000-continued": {
      "parse": {
        "median": 0.989,
        "ci95": [
          0.986,
          0.998
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    }
  },
  "openbabel_strict_wins": {
    "sdf-0": {
      "parse": false,
      "write": true,
      "roundtrip": true
    },
    "sdf-5000": {
      "parse": true,
      "write": true,
      "roundtrip": true
    },
    "v3000-0": {
      "parse": false,
      "write": false,
      "roundtrip": false
    },
    "v3000-5000": {
      "parse": false,
      "write": false,
      "roundtrip": false
    },
    "v3000-chains": {
      "parse": true,
      "write": true,
      "roundtrip": true
    },
    "v3000-continued": {
      "parse": true,
      "write": false,
      "roundtrip": false
    }
  },
  "no_timing_samples_reused": true,
  "universal_superiority": false
}
```

## Whole-process peak RSS

| Profile | Mode | Baseline MiB | Candidate MiB | Open Babel MiB |
|---|---|---:|---:|---:|
| sdf-0 | parse | 15.19 | 15.12 | 52.72 |
| sdf-0 | write | 10.66 | 10.67 | 39.69 |
| sdf-0 | roundtrip | 15.23 | 15.19 | 52.83 |
| sdf-5000 | parse | 10.95 | 10.95 | 41.31 |
| sdf-5000 | write | 7.86 | 7.83 | 33.33 |
| sdf-5000 | roundtrip | 10.94 | 10.97 | 41.42 |
| v3000-0 | parse | 16.55 | 16.48 | 51.94 |
| v3000-0 | write | 11.03 | 10.94 | 38.95 |
| v3000-0 | roundtrip | 16.59 | 16.52 | 52.02 |
| v3000-5000 | parse | 11.02 | 10.95 | 40.89 |
| v3000-5000 | write | 7.70 | 7.66 | 32.92 |
| v3000-5000 | roundtrip | 10.97 | 10.98 | 40.92 |
| v3000-chains | parse | 13.91 | 13.72 | 50.75 |
| v3000-chains | write | 9.50 | 9.23 | 39.36 |
| v3000-chains | roundtrip | 13.92 | 13.73 | 51.47 |
| v3000-continued | parse | 14.23 | 14.31 | 51.03 |
| v3000-continued | write | 9.86 | 9.81 | 39.64 |
| v3000-continued | roundtrip | 14.27 | 14.38 | 51.78 |

## Retention and limits

The immutable file manifest also includes pre-existing Python bytecode caches
copied with the baseline tree. Those files are unused by both native Rust
executables and the external RDKit observer environment. The
`manifest-entry-classification.json` inventory separates them from Git-tracked
code/build inputs; cache hashes are retained as incidental files, not claimed
as reconstructable source. All native source/build inputs match the committed
files. No timing observation or adoption decision changed.

Existing tiny test fixtures excluded by sparse checkout were restored from unchanged Git contents before testing. The timing wrapper reserves 64MiB and can retry only the same in-memory owned output write after ENOSPC; it never repeats a worker or timer.

The [2026-10-16-openbabel-v3000-borrowed-rows-artifacts](2026-10-16-openbabel-v3000-borrowed-rows-artifacts/README.md) packet retains complete raw matrices and every mismatch/negative interval, decisions, source/compiler/corpus/library/binary hashes, frozen plan, reconstructable patch, full-report oracle proof, build/test logs and launch/verification drivers. Full binaries, generated profile inputs and complete source snapshot/rewrite streams remain locally under `/private/tmp/chematic-v3000-record-stream-evidence-20261016` with versioned hashes. Source reconstruction does not imply bit-reproducible compiler output.

Additional hosts, sealed holdout chemistry, richer query/SGROUP throughput, cold CLI startup, Python/WASM throughput, published packages, other formats and FF/3D performance remain open. Correctness checks do not prove toolkit-wide superiority, and source measurements are separate from main integration, hosted CI and publication.
