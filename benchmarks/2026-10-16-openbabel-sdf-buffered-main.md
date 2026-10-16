# Borrowed SDF records after main geometry update — 2026-10-16

This candidate fails its predeclared local adoption gate. The full decision and all negative intervals are retained below. Complete speed superiority over Open Babel is not established. This is native source evidence from one shared macOS arm64 host, not a published-package comparison.

## Integration boundary

Main advanced during the first comparison. PR #924 changes the geometry classifier shared with SDF and MOL2, so the full five-profile comparison is repeated against the new source baseline before adoption. The complete [first packet](2026-10-16-openbabel-sdf-buffered.md) remains frozen. This packet uses freshly built baseline/candidate Rust executables on the new base, keeps all thresholds and inputs, and does not pool or replace earlier observations. The SDF optimization itself is unchanged.

## Implementation and preservation

The existing SDF file reader copies each physical line into an owned record buffer. The candidate scans complete newline-terminated lines in the bytes already exposed by BufRead. When a complete delimiter-terminated record fits that buffer and all existing limits allow it, the unchanged molecule/SD-field parser borrows that slice directly. The source is consumed by exactly the record and delimiter length after parsing. Empty records remain iterative. Trailing incomplete records, missing final newlines and every potential resource-limit case use the original read-until path.

The line reader and its resource-limit validator are byte-identical to the pinned baseline. The shared record parser only substitutes a byte slice and diagnostics boolean for access through self. The [preservation proof](2026-10-16-openbabel-sdf-buffered-main-artifacts/fallback-preservation.json) checks these boundaries. Chemical interpretation, diagnostics, SD properties, error precedence, failed-record continuation and logical input positions are retained. Interrupted reads use the same retry behavior. There is no new public API, unsafe code or per-record heap wrapper.

The lockfile's existing memchr 2.8.3 becomes a direct dependency for newline scanning; no new package version is introduced. All source hashes, including Cargo manifests/lockfile, are pinned. The default diagnostics path and lightweight SDMolSupplier-compatible path are checked separately. This optimization targets records that fit the caller's available buffer; it does not promise that every streaming input avoids a record copy.

## Correctness and local checks

- The final source passes all 1,254 MOL tests across 54 suites, zero failed/ignored, including all six buffer tests. These correctness tests use release opt-level 3, LTO disabled and 16 codegen units in a separate target. The measured source executables use full LTO and one codegen unit; their complete pre-timing snapshots/rewrites are separately compared exactly. Test and timing build boundaries are recorded rather than conflated.
- Differential tests match the original iterator path on 5,796 boundary/resource-limit cases and 3,774 deterministic I/O-fault/interruption cases, including consumed bytes, byte/record counters, termination, exact typed errors and full successful records.
- The two existing unfiltered 1,000-record SDF corpora contribute 20,000 record comparisons across five chunk sizes and both diagnostic modes. Actual BufReader capacities contribute another 10,000 comparisons on the first 1,000 records; 3,340 complete spans use the borrowing path and other spans exercise fallback/refills. These are repeated boundary checks, not 30,000 independent molecules.
- Before timing, both source executables have byte-identical complete internal snapshots and rewrites for all 4,064 input records across the five profiles. SDF and MOL profiles overlap in chemistry; this is an API/format coverage count, not 4,064 independent structures.
- Ten unchanged benchmark-runner gates, release all-target MOL Clippy, workspace formatting and refactor-budget checks pass. An initial test compile compared neighbor iterators directly; the corrected test compares their ordered values. Clippy also rejected an untimed custom enum draft for its large variant; the final helper uses the standard Option/Result return form, with no per-record heap wrapper or lint suppression. Earlier draft checks/failures remain separately retained in the first packet and no timing used them.

The [validation record](2026-10-16-openbabel-sdf-buffered-main-artifacts/validation.json), [exact source equality](2026-10-16-openbabel-sdf-buffered-main-artifacts/preflight-source-equality.json), [source manifest](2026-10-16-openbabel-sdf-buffered-main-artifacts/source-manifest.json) and [reconstructable patch](2026-10-16-openbabel-sdf-buffered-main-artifacts/candidate.patch.gz) preserve those claims. The patch was independently reconstructed from its named Git base and every changed-file hash verified. Local checks are separate from hosted CI, main integration and publication.

## Frozen protocol

Baseline: clean main revision 025fa2f63da324869aa144a7520626ae4ebbfd22, source version 1.1.5. Baseline and candidate native executables are freshly built with Rust 1.97.0, release opt-level 3, full LTO, one codegen unit, one build job and incremental builds disabled. The previous sparse-ID and aromatic-normalization prototypes are excluded. The unchanged C++ harness executables are reused by exact binary/library hash from the first packet, where they were freshly compiled with Apple Clang -O3 -std=c++17 against Open Babel 3.2.1. No timing samples are reused. RDKit 2026.09.1 is an independent observer, not a timed comparator. All installed Open Babel package libraries/plugins and RDKit package native files are pinned; system and other transitive libraries outside those package directories are excluded.

The two SDF profiles and two diverse V3000 profiles retain all 1,000 rows at offsets 0 and 5,000 of the existing frozen data. The fifth profile retains all 64 controlled linear chains of 16–1,024 atoms. These are exposed development corpora, not an independent holdout. Existing [corpus provenance](../validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json) and [third-party notices](../THIRD_PARTY_NOTICES.md) apply.

Each profile uses 21 sequential rotating three-arm fresh-process blocks and three native batch repetitions per lane. Input loading, process startup and initial plugin discovery are excluded. Parse retains the entire diagnostic batch; SDF framing is inside its parse timer, whereas native MOL blocks are split outside the V3000 timer. Prepared writing follows one untimed warm write. Roundtrip reads a fresh batch and writes it within one timer. The existing Rust SDF writer reuses a String; the native Open Babel writer allocates its API return strings. Whole-process high-water RSS includes loading, preparation, warmup and plugins; it is not a streaming-memory guarantee. Own builds/tests/cache cleanup finish before formal timing. Other host workloads are uncontrolled; every block records load averages.

Ratios above 1 favor the candidate. Paired median-ratio intervals use the unchanged runner's 10,000 deterministic bootstrap resamples and three-decimal reporting. Strict Open Babel wins require semantic prerequisites, CI lower >1 and 21/21 faster blocks. Parse comparator eligibility checks canonical internal graphs; rewrite eligibility checks all seven declared observer fields. Semantic mismatches remain unranked without filtering or repairs. Writer code is unchanged; writer timing movement is not an isolated writer-algorithm improvement and its cause is not established. Source equality and source-relative speed, equivalent-output comparator wins and toolkit-wide superiority are distinct conclusions.

The [predeclared plan](2026-10-16-openbabel-sdf-buffered-main-artifacts/measurement-plan.json) requires exact source snapshots/rewrites, positive parse and roundtrip confidence lower bounds on both SDF slices, and every one of all 15 SDF/V3000 source-relative lane medians >=0.99. Every negative interval is disclosed, including unchanged writer operations. There are no adaptive timing reruns, pooled samples, substitutions or relaxed gates.

## Complete timing matrix

| Profile | Operation | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 37.847 | 37.112 | 1.025 [1.018, 1.028]; 21/21 **strict win** | unranked: semantic mismatch |
| sdf-0 | write | 11.409 | 11.654 | 0.985 [0.977, 0.987]; 1/21 | 6.084 [5.920, 6.108]; 21/21 **strict win** |
| sdf-0 | roundtrip | 49.914 | 49.297 | 1.013 [1.010, 1.017]; 18/21 | 7.109 [7.083, 7.121]; 21/21 **strict win** |
| sdf-5000 | parse | 14.297 | 13.907 | 1.028 [1.026, 1.031]; 21/21 **strict win** | 8.920 [8.887, 8.991]; 21/21 **strict win** |
| sdf-5000 | write | 4.817 | 4.935 | 0.977 [0.969, 0.982]; 0/21 | 8.268 [8.157, 8.379]; 21/21 **strict win** |
| sdf-5000 | roundtrip | 19.351 | 19.027 | 1.018 [1.011, 1.022]; 18/21 | 8.678 [8.653, 8.735]; 21/21 **strict win** |
| v3000-0 | parse | 39.655 | 39.688 | 0.997 [0.991, 1.004]; 9/21 | unranked: semantic mismatch |
| v3000-0 | write | 8.380 | 8.361 | 1.002 [0.996, 1.012]; 12/21 | unranked: semantic mismatch |
| v3000-0 | roundtrip | 48.645 | 48.806 | 0.999 [0.995, 1.002]; 9/21 | unranked: semantic mismatch |
| v3000-5000 | parse | 15.444 | 15.509 | 0.996 [0.992, 0.999]; 6/21 | unranked: semantic mismatch |
| v3000-5000 | write | 2.161 | 2.163 | 1.000 [0.993, 1.004]; 11/21 | unranked: semantic mismatch |
| v3000-5000 | roundtrip | 17.660 | 17.706 | 0.996 [0.992, 0.999]; 6/21 | unranked: semantic mismatch |
| v3000-chains | parse | 21.178 | 20.702 | 1.006 [0.992, 1.044]; 11/21 | 143.475 [142.288, 146.787]; 21/21 **strict win** |
| v3000-chains | write | 4.606 | 4.613 | 1.006 [0.990, 1.020]; 12/21 | 15.372 [15.166, 15.812]; 21/21 **strict win** |
| v3000-chains | roundtrip | 25.288 | 25.680 | 0.998 [0.991, 1.006]; 10/21 | 120.925 [119.193, 122.219]; 21/21 **strict win** |

## Source adoption decision

```json
{
  "exact_source_snapshots_and_rewrites": true,
  "sdf_primary": {
    "sdf-0": {
      "parse": {
        "median": 1.025,
        "ci95": [
          1.018,
          1.028
        ],
        "winning_blocks": 21,
        "blocks": 21
      },
      "roundtrip": {
        "median": 1.013,
        "ci95": [
          1.01,
          1.017
        ],
        "winning_blocks": 18,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "parse": {
        "median": 1.028,
        "ci95": [
          1.026,
          1.031
        ],
        "winning_blocks": 21,
        "blocks": 21
      },
      "roundtrip": {
        "median": 1.018,
        "ci95": [
          1.011,
          1.022
        ],
        "winning_blocks": 18,
        "blocks": 21
      }
    }
  },
  "sdf_primary_ci_pass": true,
  "all_15_source_lane_median_floor": 0.99,
  "below_floor": {
    "sdf-0": {
      "write": {
        "median": 0.985,
        "ci95": [
          0.977,
          0.987
        ],
        "winning_blocks": 1,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "write": {
        "median": 0.977,
        "ci95": [
          0.969,
          0.982
        ],
        "winning_blocks": 0,
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
        "median": 0.985,
        "ci95": [
          0.977,
          0.987
        ],
        "winning_blocks": 1,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "write": {
        "median": 0.977,
        "ci95": [
          0.969,
          0.982
        ],
        "winning_blocks": 0,
        "blocks": 21
      }
    },
    "v3000-0": {},
    "v3000-5000": {
      "parse": {
        "median": 0.996,
        "ci95": [
          0.992,
          0.999
        ],
        "winning_blocks": 6,
        "blocks": 21
      },
      "roundtrip": {
        "median": 0.996,
        "ci95": [
          0.992,
          0.999
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    },
    "v3000-chains": {}
  },
  "local_adoption_gate_passed": false,
  "no_adaptive_timing_reruns": true
}
```

## Whole-process peak RSS

| Profile | Operation | Baseline MiB | Candidate MiB | Open Babel MiB |
|---|---|---:|---:|---:|
| sdf-0 | parse | 15.20 | 15.19 | 52.80 |
| sdf-0 | write | 10.67 | 10.69 | 39.69 |
| sdf-0 | roundtrip | 15.27 | 15.27 | 52.84 |
| sdf-5000 | parse | 11.03 | 11.05 | 41.33 |
| sdf-5000 | write | 7.84 | 7.94 | 33.33 |
| sdf-5000 | roundtrip | 11.05 | 11.06 | 41.41 |
| v3000-0 | parse | 16.56 | 16.58 | 51.92 |
| v3000-0 | write | 11.03 | 11.06 | 38.97 |
| v3000-0 | roundtrip | 16.55 | 16.59 | 52.03 |
| v3000-5000 | parse | 11.14 | 11.06 | 40.92 |
| v3000-5000 | write | 7.84 | 7.77 | 32.92 |
| v3000-5000 | roundtrip | 11.14 | 11.09 | 40.94 |
| v3000-chains | parse | 13.94 | 13.84 | 50.73 |
| v3000-chains | write | 9.44 | 9.55 | 39.36 |
| v3000-chains | roundtrip | 14.02 | 13.95 | 51.52 |

## Reproduction and limits

Check out the named baseline in an isolated worktree. Decompress/apply the candidate patch for a separate candidate. Build both using CARGO_BUILD_JOBS=1 CARGO_INCREMENTAL=0 cargo build --release --locked -p chematic-mol --example openbabel_sdf_throughput --example openbabel_v3000_throughput; preserve baseline executables before rebuilding. Compile the unchanged C++ harnesses with the recorded pkg-config Open Babel flags. Decompress the five named checked-in corpora. Run the unchanged scripts/bench_openbabel_sdf_throughput.py and scripts/bench_openbabel_v3000_throughput.py with --baseline, --candidate, --openbabel, --corpus, --records <1000-or-64>, --blocks 21 --repeats 3 --source-manifest <manifest.json> --output <report.json> in the pinned observer environment.

The retained launch scripts check complete source/library/binary/corpus hashes before and after each profile. Their 64 MiB reserve only retries an owned checkpoint write with the same in-memory observation on ENOSPC; workers, timers and statistics are unchanged. Darwin time -l needs OS-statistics access. Adjust local paths when reproducing.

Raw complete [SDF 0](2026-10-16-openbabel-sdf-buffered-main-artifacts/paired-sdf-0.json.gz), [SDF 5000](2026-10-16-openbabel-sdf-buffered-main-artifacts/paired-sdf-5000.json.gz), [MOL 0](2026-10-16-openbabel-sdf-buffered-main-artifacts/paired-v3000-0.json.gz), [MOL 5000](2026-10-16-openbabel-sdf-buffered-main-artifacts/paired-v3000-5000.json.gz) and [chains](2026-10-16-openbabel-sdf-buffered-main-artifacts/paired-v3000-chains.json.gz) matrices contain every block, mismatch, hash and RSS observation. Build/test/format/refactor/Clippy/timing logs and the exact launch scripts are retained here. Full executables and snapshot/rewrite streams remain locally in /private/tmp/chematic-sdf-buffer-main-evidence-20261016; their [hashes](2026-10-16-openbabel-sdf-buffered-main-artifacts/local-semantic-streams.json) are versioned without duplicating binaries and generated outputs into Git.

Other formats, diverse held-out chemistry, rich query/SGROUP performance, buffers smaller than a record, filesystem/CLI/Python overhead, cold startup, multithreading, additional hosts and published artifacts remain separate comparisons. Existing Open Babel semantic disagreements are not resolved by a faster SDF record boundary. This packet does not validate later main commits or published packages.
