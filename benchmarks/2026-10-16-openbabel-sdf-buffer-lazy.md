# Lazy owned storage for borrowed SDF records — 2026-10-16

This candidate fails its predeclared local adoption gate. The full decision and all negative intervals are retained below. Complete speed superiority over Open Babel is not established. This is native source evidence from one shared macOS arm64 host, not a published-package comparison.

## Integration boundary

Main advanced during the first comparison. PR #924 changes the geometry classifier shared with SDF and MOL2, so the full five-profile comparison is repeated against the new source baseline before adoption. The complete [first packet](2026-10-16-openbabel-sdf-buffered.md) remains frozen. This packet uses freshly built baseline/candidate Rust executables on the new base, keeps all thresholds and inputs, and does not pool or replace earlier observations. The [first integration candidate](2026-10-16-openbabel-sdf-buffered-main.md) retains parse/roundtrip gains but regresses prepared SDF writing by 1.5–2.3%, failing the same floor. This third, distinct source candidate defers the unused owned-buffer allocation until fallback. All earlier matrices remain separate and no timing is reused.

## Implementation and preservation

The existing SDF file reader copies each physical line into an owned record buffer. The candidate scans complete newline-terminated lines in the bytes already exposed by BufRead. When a complete delimiter-terminated record fits that buffer and all existing limits allow it, the unchanged molecule/SD-field parser borrows that slice directly. The source is consumed by exactly the record and delimiter length after parsing. Complete borrowed records now allocate no owned framing buffer. This removes an unused allocation, without assuming that allocation caused the preceding writer regression. Empty records remain iterative. Trailing incomplete records, missing final newlines and every potential resource-limit case use the original read-until path.

The constructor starts with Vec::new instead of reserving 2,048 bytes. The original line fallback reserves that storage only when needed, after clearing the reused buffer. Removing that single deferred reservation from the fallback makes the remaining line reader and resource-limit validator byte-identical to the pinned baseline. The oracle independently freezes the original line reader and restores its original initial reservation; it does not share the changed fallback storage path. The shared record parser only substitutes a byte slice and diagnostics boolean for access through self. The [preservation proof](#raw-artifacts-retired) checks these boundaries. Chemical interpretation, diagnostics, SD properties, error precedence, failed-record continuation and logical input positions are retained. Interrupted reads use the same retry behavior. There is no new public API, unsafe code or per-record heap wrapper.

The lockfile's existing memchr 2.8.3 becomes a direct dependency for newline scanning; no new package version is introduced. All source hashes, including Cargo manifests/lockfile, are pinned. The default diagnostics path and lightweight SDMolSupplier-compatible path are checked separately. This optimization targets records that fit the caller's available buffer; it does not promise that every streaming input avoids a record copy.

## Correctness and local checks

- The final source passes all 1,254 MOL tests across 54 suites, zero failed/ignored, including all six buffer tests. These correctness tests use release opt-level 3, LTO disabled and 16 codegen units in a separate target. The measured source executables use full LTO and one codegen unit; their complete pre-timing snapshots/rewrites are separately compared exactly. Test and timing build boundaries are recorded rather than conflated.
- Differential tests match the original iterator path on 5,796 boundary/resource-limit cases and 3,774 deterministic I/O-fault/interruption cases, including consumed bytes, byte/record counters, termination, exact typed errors and full successful records.
- The two existing unfiltered 1,000-record SDF corpora contribute 20,000 record comparisons across five chunk sizes and both diagnostic modes. Actual BufReader capacities contribute another 10,000 comparisons on the first 1,000 records; 3,340 complete spans use the borrowing path and other spans exercise fallback/refills. These are repeated boundary checks, not 30,000 independent molecules.
- Before timing, both source executables have byte-identical complete internal snapshots and rewrites for all 4,064 input records across the five profiles. SDF and MOL profiles overlap in chemistry; this is an API/format coverage count, not 4,064 independent structures.
- Ten unchanged benchmark-runner gates, release all-target MOL Clippy, workspace formatting and refactor-budget checks pass. An initial test compile compared neighbor iterators directly; the corrected test compares their ordered values. Clippy also rejected an untimed custom enum draft for its large variant; the final helper uses the standard Option/Result return form, with no per-record heap wrapper or lint suppression. Earlier draft checks/failures remain separately retained in the first packet and no timing used them.

The [validation record](#raw-artifacts-retired), [exact source equality](#raw-artifacts-retired), [source manifest](#raw-artifacts-retired) and [reconstructable patch](#raw-artifacts-retired) preserve those claims. The patch was independently reconstructed from its named Git base and every changed-file hash verified. Local checks are separate from hosted CI, main integration and publication.

## Frozen protocol

Baseline: clean main revision 025fa2f63da324869aa144a7520626ae4ebbfd22, source version 1.1.5. The clean integration-baseline executables are reused by exact binary hash from the preceding packet, where they were freshly built from the named commit. The candidate is freshly built. Both use Rust 1.97.0, release opt-level 3, full LTO, one codegen unit, one build job and incremental builds disabled. No timing samples are reused. The previous sparse-ID and aromatic-normalization prototypes are excluded. The unchanged C++ harness executables are reused by exact binary/library hash from the first packet, where they were freshly compiled with Apple Clang -O3 -std=c++17 against Open Babel 3.2.1. No timing samples are reused. RDKit 2026.09.1 is an independent observer, not a timed comparator. All installed Open Babel package libraries/plugins and RDKit package native files are pinned; system and other transitive libraries outside those package directories are excluded.

The two SDF profiles and two diverse V3000 profiles retain all 1,000 rows at offsets 0 and 5,000 of the existing frozen data. The fifth profile retains all 64 controlled linear chains of 16–1,024 atoms. These are exposed development corpora, not an independent holdout. Existing [corpus provenance](../validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json) and [third-party notices](../THIRD_PARTY_NOTICES.md) apply.

Each profile uses 21 sequential rotating three-arm fresh-process blocks and three native batch repetitions per lane. Input loading, process startup and initial plugin discovery are excluded. Parse retains the entire diagnostic batch; SDF framing is inside its parse timer, whereas native MOL blocks are split outside the V3000 timer. Prepared writing follows one untimed warm write. Roundtrip reads a fresh batch and writes it within one timer. The existing Rust SDF writer reuses a String; the native Open Babel writer allocates its API return strings. Whole-process high-water RSS includes loading, preparation, warmup and plugins; it is not a streaming-memory guarantee. Own builds/tests/cache cleanup finish before formal timing. Other host workloads are uncontrolled; every block records load averages.

Ratios above 1 favor the candidate. Paired median-ratio intervals use the unchanged runner's 10,000 deterministic bootstrap resamples and three-decimal reporting. Strict Open Babel wins require semantic prerequisites, CI lower >1 and 21/21 faster blocks. Parse comparator eligibility checks canonical internal graphs; rewrite eligibility checks all seven declared observer fields. Semantic mismatches remain unranked without filtering or repairs. Writer code is unchanged; writer timing movement is not an isolated writer-algorithm improvement and its cause is not established. Source equality and source-relative speed, equivalent-output comparator wins and toolkit-wide superiority are distinct conclusions.

The [predeclared plan](#raw-artifacts-retired) requires exact source snapshots/rewrites, positive parse and roundtrip confidence lower bounds on both SDF slices, and every one of all 15 SDF/V3000 source-relative lane medians >=0.99. Every negative interval is disclosed, including unchanged writer operations. There are no adaptive timing reruns, pooled samples, substitutions or relaxed gates.

## Complete timing matrix

| Profile | Operation | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 38.650 | 37.750 | 1.020 [1.016, 1.026]; 21/21 **strict win** | unranked: semantic mismatch |
| sdf-0 | write | 11.641 | 11.838 | 0.983 [0.974, 0.992]; 6/21 | 6.018 [5.940, 6.069]; 21/21 **strict win** |
| sdf-0 | roundtrip | 50.799 | 50.275 | 1.013 [1.012, 1.021]; 20/21 | 7.149 [7.100, 7.184]; 21/21 **strict win** |
| sdf-5000 | parse | 14.534 | 14.189 | 1.026 [1.024, 1.032]; 21/21 **strict win** | 8.937 [8.882, 9.019]; 21/21 **strict win** |
| sdf-5000 | write | 4.892 | 4.988 | 0.980 [0.969, 0.985]; 4/21 | 8.334 [8.168, 8.417]; 21/21 **strict win** |
| sdf-5000 | roundtrip | 19.681 | 19.331 | 1.017 [1.012, 1.026]; 19/21 | 8.724 [8.697, 8.784]; 21/21 **strict win** |
| v3000-0 | parse | 39.263 | 39.516 | 0.998 [0.991, 1.013]; 10/21 | unranked: semantic mismatch |
| v3000-0 | write | 8.310 | 8.383 | 0.995 [0.986, 1.014]; 9/21 | unranked: semantic mismatch |
| v3000-0 | roundtrip | 48.467 | 48.377 | 1.001 [0.989, 1.012]; 11/21 | unranked: semantic mismatch |
| v3000-5000 | parse | 19.357 | 19.674 | 0.999 [0.970, 1.006]; 10/21 | unranked: semantic mismatch |
| v3000-5000 | write | 2.698 | 2.754 | 0.993 [0.974, 1.002]; 7/21 | unranked: semantic mismatch |
| v3000-5000 | roundtrip | 22.896 | 22.742 | 0.994 [0.989, 1.011]; 9/21 | unranked: semantic mismatch |
| v3000-chains | parse | 17.304 | 17.226 | 1.000 [0.988, 1.013]; 11/21 | 145.065 [142.727, 145.909]; 21/21 **strict win** |
| v3000-chains | write | 3.801 | 3.862 | 0.994 [0.980, 1.014]; 9/21 | 15.460 [15.201, 15.581]; 21/21 **strict win** |
| v3000-chains | roundtrip | 21.201 | 21.220 | 0.998 [0.994, 1.002]; 9/21 | 120.682 [118.954, 121.810]; 21/21 **strict win** |

## Source adoption decision

```json
{
  "exact_source_snapshots_and_rewrites": true,
  "sdf_primary": {
    "sdf-0": {
      "parse": {
        "median": 1.02,
        "ci95": [
          1.016,
          1.026
        ],
        "winning_blocks": 21,
        "blocks": 21
      },
      "roundtrip": {
        "median": 1.013,
        "ci95": [
          1.012,
          1.021
        ],
        "winning_blocks": 20,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "parse": {
        "median": 1.026,
        "ci95": [
          1.024,
          1.032
        ],
        "winning_blocks": 21,
        "blocks": 21
      },
      "roundtrip": {
        "median": 1.017,
        "ci95": [
          1.012,
          1.026
        ],
        "winning_blocks": 19,
        "blocks": 21
      }
    }
  },
  "sdf_primary_ci_pass": true,
  "all_15_source_lane_median_floor": 0.99,
  "below_floor": {
    "sdf-0": {
      "write": {
        "median": 0.983,
        "ci95": [
          0.974,
          0.992
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "write": {
        "median": 0.98,
        "ci95": [
          0.969,
          0.985
        ],
        "winning_blocks": 4,
        "blocks": 21
      }
    },
    "v3000-0": {},
    "v3000-5000": {},
    "v3000-chains": {}
  },
  "all_negative_intervals": {
    "sdf-0": {
      "write": {
        "median": 0.983,
        "ci95": [
          0.974,
          0.992
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "write": {
        "median": 0.98,
        "ci95": [
          0.969,
          0.985
        ],
        "winning_blocks": 4,
        "blocks": 21
      }
    },
    "v3000-0": {},
    "v3000-5000": {},
    "v3000-chains": {}
  },
  "local_adoption_gate_passed": false,
  "no_adaptive_timing_reruns": true
}
```

## Whole-process peak RSS

| Profile | Operation | Baseline MiB | Candidate MiB | Open Babel MiB |
|---|---|---:|---:|---:|
| sdf-0 | parse | 15.23 | 15.17 | 52.75 |
| sdf-0 | write | 10.69 | 10.62 | 39.69 |
| sdf-0 | roundtrip | 15.27 | 15.19 | 52.92 |
| sdf-5000 | parse | 11.03 | 10.97 | 41.34 |
| sdf-5000 | write | 7.83 | 7.84 | 33.31 |
| sdf-5000 | roundtrip | 11.02 | 11.00 | 41.42 |
| v3000-0 | parse | 16.52 | 16.59 | 51.94 |
| v3000-0 | write | 11.03 | 11.08 | 38.97 |
| v3000-0 | roundtrip | 16.55 | 16.61 | 52.05 |
| v3000-5000 | parse | 11.14 | 11.14 | 40.97 |
| v3000-5000 | write | 7.86 | 7.78 | 32.95 |
| v3000-5000 | roundtrip | 11.16 | 11.11 | 41.00 |
| v3000-chains | parse | 13.98 | 13.88 | 50.73 |
| v3000-chains | write | 9.44 | 9.53 | 39.38 |
| v3000-chains | roundtrip | 14.02 | 13.92 | 51.47 |

## Reproduction and limits

Check out the named baseline in an isolated worktree. Decompress/apply the candidate patch for a separate candidate. Build both using CARGO_BUILD_JOBS=1 CARGO_INCREMENTAL=0 cargo build --release --locked -p chematic-mol --example openbabel_sdf_throughput --example openbabel_v3000_throughput; preserve baseline executables before rebuilding. Compile the unchanged C++ harnesses with the recorded pkg-config Open Babel flags. Decompress the five named checked-in corpora. Run the unchanged scripts/bench_openbabel_sdf_throughput.py and scripts/bench_openbabel_v3000_throughput.py with --baseline, --candidate, --openbabel, --corpus, --records <1000-or-64>, --blocks 21 --repeats 3 --source-manifest <manifest.json> --output <report.json> in the pinned observer environment.

The retained launch scripts check complete source/library/binary/corpus hashes before and after each profile. Their 64 MiB reserve only retries an owned checkpoint write with the same in-memory observation on ENOSPC; workers, timers and statistics are unchanged. Darwin time -l needs OS-statistics access. Adjust local paths when reproducing.

Raw complete [SDF 0](#raw-artifacts-retired), [SDF 5000](#raw-artifacts-retired), [MOL 0](#raw-artifacts-retired), [MOL 5000](#raw-artifacts-retired) and [chains](#raw-artifacts-retired) matrices contain every block, mismatch, hash and RSS observation. Build/test/format/refactor/Clippy/timing logs and the exact launch scripts are retained here. Full executables and snapshot/rewrite streams remain locally in /private/tmp/chematic-sdf-buffer-lazy-evidence-20261016; their [hashes](#raw-artifacts-retired) are versioned without duplicating binaries and generated outputs into Git.

Other formats, diverse held-out chemistry, rich query/SGROUP performance, buffers smaller than a record, filesystem/CLI/Python overhead, cold startup, multithreading, additional hosts and published artifacts remain separate comparisons. Existing Open Babel semantic disagreements are not resolved by a faster SDF record boundary. This packet does not validate later main commits or published packages.


## Raw artifacts retired

The checked-in wheel, source-tarball, and intermediate-result directory for this historical experiment was removed on 2026-10-23. The summary and conclusions above remain, and the raw files can be recovered from Git history before the benchmark cleanup.
