# Borrowed complete SDF records — 2026-10-16

This candidate passes its predeclared local adoption gate. The full decision and all negative intervals are retained below. Complete speed superiority over Open Babel is not established. This is native source evidence from one shared macOS arm64 host, not a published-package comparison.

## Implementation and preservation

The existing SDF file reader copies each physical line into an owned record buffer. The candidate scans complete newline-terminated lines in the bytes already exposed by BufRead. When a complete delimiter-terminated record fits that buffer and all existing limits allow it, the unchanged molecule/SD-field parser borrows that slice directly. The source is consumed by exactly the record and delimiter length after parsing. Empty records remain iterative. Trailing incomplete records, missing final newlines and every potential resource-limit case use the original read-until path.

The line reader and its resource-limit validator are byte-identical to the pinned baseline. The shared record parser only substitutes a byte slice and diagnostics boolean for access through self. The [preservation proof](#raw-artifacts-retired) checks these boundaries. Chemical interpretation, diagnostics, SD properties, error precedence, failed-record continuation and logical input positions are retained. Interrupted reads use the same retry behavior. There is no new public API, unsafe code or per-record heap wrapper.

The lockfile's existing memchr 2.8.3 becomes a direct dependency for newline scanning; no new package version is introduced. All source hashes, including Cargo manifests/lockfile, are pinned. The default diagnostics path and lightweight SDMolSupplier-compatible path are checked separately. This optimization targets records that fit the caller's available buffer; it does not promise that every streaming input avoids a record copy.

## Correctness and local checks

- The final source passes all 1,251 MOL tests across 53 suites, zero failed/ignored, including all six buffer tests. These correctness tests use release opt-level 3, LTO disabled and 16 codegen units in a separate target. The measured source executables use full LTO and one codegen unit; their complete pre-timing snapshots/rewrites are separately compared exactly. Test and timing build boundaries are recorded rather than conflated.
- Differential tests match the original iterator path on 5,796 boundary/resource-limit cases and 3,774 deterministic I/O-fault/interruption cases, including consumed bytes, byte/record counters, termination, exact typed errors and full successful records.
- The two existing unfiltered 1,000-record SDF corpora contribute 20,000 record comparisons across five chunk sizes and both diagnostic modes. Actual BufReader capacities contribute another 10,000 comparisons on the first 1,000 records; 3,340 complete spans use the borrowing path and other spans exercise fallback/refills. These are repeated boundary checks, not 30,000 independent molecules.
- Before timing, both source executables have byte-identical complete internal snapshots and rewrites for all 4,064 input records across the five profiles. SDF and MOL profiles overlap in chemistry; this is an API/format coverage count, not 4,064 independent structures.
- Ten unchanged benchmark-runner gates, release all-target MOL Clippy, workspace formatting and refactor-budget checks pass. An initial test compile compared neighbor iterators directly; the corrected test compares their ordered values. Clippy also rejected an untimed custom enum draft for its large variant; the final helper uses the standard Option/Result return form, with no per-record heap wrapper or lint suppression. Earlier draft checks/failures remain separately retained and no timing used them.

The [validation record](#raw-artifacts-retired), [exact source equality](#raw-artifacts-retired), [source manifest](#raw-artifacts-retired) and [reconstructable patch](#raw-artifacts-retired) preserve those claims. The patch was independently reconstructed from its named Git base and every changed-file hash verified. Local checks are separate from hosted CI, main integration and publication.

## Frozen protocol

Baseline: clean main revision b9f1e3cdd959a20f1cfa4c635237418a550ffa00, source version 1.1.5. Baseline and candidate native executables are freshly built with Rust 1.97.0, release opt-level 3, full LTO, one codegen unit, one build job and incremental builds disabled. The previous sparse-ID and aromatic-normalization prototypes are excluded. The C++ harnesses are unchanged, freshly compiled with Apple Clang -O3 -std=c++17 against Open Babel 3.2.1. RDKit 2026.09.1 is an independent observer, not a timed comparator. All installed Open Babel package libraries/plugins and RDKit package native files are pinned; system and other transitive libraries outside those package directories are excluded.

The two SDF profiles and two diverse V3000 profiles retain all 1,000 rows at offsets 0 and 5,000 of the existing frozen data. The fifth profile retains all 64 controlled linear chains of 16–1,024 atoms. These are exposed development corpora, not an independent holdout. Existing [corpus provenance](../validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json) and [third-party notices](../THIRD_PARTY_NOTICES.md) apply.

Each profile uses 21 sequential rotating three-arm fresh-process blocks and three native batch repetitions per lane. Input loading, process startup and initial plugin discovery are excluded. Parse retains the entire diagnostic batch; SDF framing is inside its parse timer, whereas native MOL blocks are split outside the V3000 timer. Prepared writing follows one untimed warm write. Roundtrip reads a fresh batch and writes it within one timer. The existing Rust SDF writer reuses a String; the native Open Babel writer allocates its API return strings. Whole-process high-water RSS includes loading, preparation, warmup and plugins; it is not a streaming-memory guarantee. Own builds/tests/cache cleanup finish before formal timing. Other host workloads are uncontrolled; every block records load averages.

Ratios above 1 favor the candidate. Paired median-ratio intervals use the unchanged runner's 10,000 deterministic bootstrap resamples and three-decimal reporting. Strict Open Babel wins require semantic prerequisites, CI lower >1 and 21/21 faster blocks. Parse comparator eligibility checks canonical internal graphs; rewrite eligibility checks all seven declared observer fields. Semantic mismatches remain unranked without filtering or repairs. Writer code is unchanged; writer timing movement is not an isolated writer-algorithm improvement and its cause is not established. Source equality and source-relative speed, equivalent-output comparator wins and toolkit-wide superiority are distinct conclusions.

The [predeclared plan](#raw-artifacts-retired) requires exact source snapshots/rewrites, positive parse and roundtrip confidence lower bounds on both SDF slices, and every one of all 15 SDF/V3000 source-relative lane medians >=0.99. Every negative interval is disclosed, including unchanged writer operations. There are no adaptive timing reruns, pooled samples, substitutions or relaxed gates.

## Complete timing matrix

| Profile | Operation | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 43.760 | 43.096 | 1.019 [1.017, 1.025]; 20/21 | unranked: semantic mismatch |
| sdf-0 | write | 13.621 | 13.316 | 1.022 [1.011, 1.030]; 17/21 | 6.109 [6.072, 6.196]; 21/21 **strict win** |
| sdf-0 | roundtrip | 58.389 | 57.133 | 1.024 [1.020, 1.030]; 21/21 **strict win** | 7.217 [7.155, 7.251]; 21/21 **strict win** |
| sdf-5000 | parse | 17.373 | 16.875 | 1.030 [1.023, 1.040]; 18/21 | 9.105 [9.075, 9.262]; 21/21 **strict win** |
| sdf-5000 | write | 6.052 | 5.737 | 1.045 [1.037, 1.057]; 19/21 | 8.731 [8.636, 8.785]; 21/21 **strict win** |
| sdf-5000 | roundtrip | 23.685 | 22.804 | 1.028 [1.021, 1.043]; 20/21 | 8.866 [8.734, 9.016]; 21/21 **strict win** |
| v3000-0 | parse | 38.927 | 39.023 | 0.995 [0.993, 1.003]; 8/21 | unranked: semantic mismatch |
| v3000-0 | write | 8.203 | 8.205 | 0.998 [0.990, 1.000]; 9/21 | unranked: semantic mismatch |
| v3000-0 | roundtrip | 47.620 | 47.883 | 0.995 [0.990, 0.999]; 6/21 | unranked: semantic mismatch |
| v3000-5000 | parse | 18.454 | 17.804 | 1.009 [0.977, 1.024]; 12/21 | unranked: semantic mismatch |
| v3000-5000 | write | 2.661 | 2.479 | 1.009 [0.998, 1.034]; 14/21 | unranked: semantic mismatch |
| v3000-5000 | roundtrip | 21.191 | 20.867 | 1.000 [0.984, 1.024]; 11/21 | unranked: semantic mismatch |
| v3000-chains | parse | 17.041 | 17.036 | 1.000 [0.991, 1.005]; 11/21 | 144.159 [143.276, 144.844]; 21/21 **strict win** |
| v3000-chains | write | 3.820 | 3.777 | 1.006 [0.992, 1.020]; 14/21 | 15.174 [15.051, 15.291]; 21/21 **strict win** |
| v3000-chains | roundtrip | 20.824 | 20.890 | 0.996 [0.989, 1.003]; 9/21 | 119.669 [119.021, 120.167]; 21/21 **strict win** |

## Source adoption decision

```json
{
  "exact_source_snapshots_and_rewrites": true,
  "sdf_primary": {
    "sdf-0": {
      "parse": {
        "median": 1.019,
        "ci95": [
          1.017,
          1.025
        ],
        "winning_blocks": 20,
        "blocks": 21
      },
      "roundtrip": {
        "median": 1.024,
        "ci95": [
          1.02,
          1.03
        ],
        "winning_blocks": 21,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "parse": {
        "median": 1.03,
        "ci95": [
          1.023,
          1.04
        ],
        "winning_blocks": 18,
        "blocks": 21
      },
      "roundtrip": {
        "median": 1.028,
        "ci95": [
          1.021,
          1.043
        ],
        "winning_blocks": 20,
        "blocks": 21
      }
    }
  },
  "sdf_primary_ci_pass": true,
  "all_15_source_lane_median_floor": 0.99,
  "below_floor": {
    "sdf-0": {},
    "sdf-5000": {},
    "v3000-0": {},
    "v3000-5000": {},
    "v3000-chains": {}
  },
  "all_negative_intervals": {
    "sdf-0": {},
    "sdf-5000": {},
    "v3000-0": {
      "roundtrip": {
        "median": 0.995,
        "ci95": [
          0.99,
          0.999
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    },
    "v3000-5000": {},
    "v3000-chains": {}
  },
  "local_adoption_gate_passed": true,
  "no_adaptive_timing_reruns": true
}
```

## Whole-process peak RSS

| Profile | Operation | Baseline MiB | Candidate MiB | Open Babel MiB |
|---|---|---:|---:|---:|
| sdf-0 | parse | 15.16 | 15.22 | 52.69 |
| sdf-0 | write | 10.61 | 10.80 | 39.67 |
| sdf-0 | roundtrip | 15.16 | 15.25 | 52.88 |
| sdf-5000 | parse | 10.95 | 11.02 | 41.33 |
| sdf-5000 | write | 7.84 | 7.94 | 33.30 |
| sdf-5000 | roundtrip | 10.95 | 11.05 | 41.41 |
| v3000-0 | parse | 16.59 | 16.48 | 51.94 |
| v3000-0 | write | 11.05 | 10.92 | 38.97 |
| v3000-0 | roundtrip | 16.62 | 16.52 | 52.03 |
| v3000-5000 | parse | 11.12 | 11.05 | 40.92 |
| v3000-5000 | write | 7.88 | 7.80 | 32.95 |
| v3000-5000 | roundtrip | 11.12 | 11.11 | 40.97 |
| v3000-chains | parse | 13.97 | 13.92 | 50.70 |
| v3000-chains | write | 9.47 | 9.38 | 39.36 |
| v3000-chains | roundtrip | 14.05 | 13.97 | 51.48 |

## Reproduction and limits

Check out the named baseline in an isolated worktree. Decompress/apply the candidate patch for a separate candidate. Build both using CARGO_BUILD_JOBS=1 CARGO_INCREMENTAL=0 cargo build --release --locked -p chematic-mol --example openbabel_sdf_throughput --example openbabel_v3000_throughput; preserve baseline executables before rebuilding. Compile the unchanged C++ harnesses with the recorded pkg-config Open Babel flags. Decompress the five named checked-in corpora. Run the unchanged scripts/bench_openbabel_sdf_throughput.py and scripts/bench_openbabel_v3000_throughput.py with --baseline, --candidate, --openbabel, --corpus, --records <1000-or-64>, --blocks 21 --repeats 3 --source-manifest <manifest.json> --output <report.json> in the pinned observer environment.

The retained launch scripts check complete source/library/binary/corpus hashes before and after each profile. Their 64 MiB reserve only retries an owned checkpoint write with the same in-memory observation on ENOSPC; workers, timers and statistics are unchanged. Darwin time -l needs OS-statistics access. Adjust local paths when reproducing.

Raw complete [SDF 0](#raw-artifacts-retired), [SDF 5000](#raw-artifacts-retired), [MOL 0](#raw-artifacts-retired), [MOL 5000](#raw-artifacts-retired) and [chains](#raw-artifacts-retired) matrices contain every block, mismatch, hash and RSS observation. Build/test/format/refactor/Clippy/timing logs and the exact launch scripts are retained here. Full executables and snapshot/rewrite streams remain locally in /private/tmp/chematic-sdf-buffer-scan-evidence-20261016; their [hashes](#raw-artifacts-retired) are versioned without duplicating binaries and generated outputs into Git.

Other formats, diverse held-out chemistry, rich query/SGROUP performance, buffers smaller than a record, filesystem/CLI/Python overhead, cold startup, multithreading, additional hosts and published artifacts remain separate comparisons. Existing Open Babel semantic disagreements are not resolved by a faster SDF record boundary. This packet does not validate later main commits or published packages.


## Raw artifacts retired

The checked-in wheel, source-tarball, and intermediate-result directory for this historical experiment was removed on 2026-10-23. The summary and conclusions above remain, and the raw files can be recovered from Git history before the benchmark cleanup.
