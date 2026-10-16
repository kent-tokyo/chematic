# Borrowed SDF reads and requested-only stereo-loss collection — 2026-10-16

This candidate passes its predeclared local adoption gate. The full decision and all negative intervals are retained below. Complete speed superiority over Open Babel is not established. This is native source evidence from one shared macOS arm64 host, not a published-package comparison.

## Integration boundary

Main advanced during the first comparison. PR #924 changes the geometry classifier shared with SDF and MOL2, so the full five-profile comparison is repeated against the new source baseline before adoption. The complete [first packet](2026-10-16-openbabel-sdf-buffered.md) remains frozen. This packet uses freshly built baseline/candidate Rust executables on the new base, keeps all thresholds and inputs, and does not pool or replace earlier observations. The [first integration candidate](2026-10-16-openbabel-sdf-buffered-main.md) retains parse/roundtrip gains but regresses prepared SDF writing by 1.5–2.3%, failing the same floor. The [lazy-buffer candidate](2026-10-16-openbabel-sdf-buffer-lazy.md) also retains a 1.7–2.0% prepared-write retreat. The [direct-field candidate](2026-10-16-openbabel-sdf-direct-fields.md) still misses the floor for first-slice writing and diverse V3000 parsing. This fifth, distinct source candidate additionally skips unused stereo-loss collection for plain output; reporting APIs explicitly retain it. All earlier matrices remain separate and no timing is reused.

## Implementation and preservation

The existing SDF file reader copies each physical line into an owned record buffer. The candidate scans complete newline-terminated lines in the bytes already exposed by BufRead. When a complete delimiter-terminated record fits that buffer and all existing limits allow it, the unchanged molecule/SD-field parser borrows that slice directly. The source is consumed by exactly the record and delimiter length after parsing. Complete borrowed records now allocate no owned framing buffer. This removes an unused allocation, without assuming that allocation caused the preceding writer regression. Empty records remain iterative. Trailing incomplete records, missing final newlines and every potential resource-limit case use the original read-until path.

The constructor starts with Vec::new instead of reserving 2,048 bytes. The original line fallback reserves that storage only when needed, after clearing the reused buffer. Removing that single deferred reservation from the fallback makes the remaining line reader and resource-limit validator byte-identical to the pinned baseline. The oracle independently freezes the original line reader and restores its original initial reservation; it does not share the changed fallback storage path. The shared record parser only substitutes a byte slice and diagnostics boolean for access through self. The [preservation proof](2026-10-16-openbabel-sdf-loss-work-artifacts/fallback-preservation.json) checks these boundaries. Chemical interpretation, diagnostics, SD properties, error precedence, failed-record continuation and logical input positions are retained. Interrupted reads use the same retry behavior. There is no new public API, unsafe code or per-record heap wrapper.

The lockfile's existing memchr 2.8.3 becomes a direct dependency for newline scanning; no new package version is introduced. All source hashes, including Cargo manifests/lockfile, are pinned. The default diagnostics path and lightweight SDMolSupplier-compatible path are checked separately. This optimization targets records that fit the caller's available buffer; it does not promise that every streaming input avoids a record copy.

The fourth delta appends the same literal header, key, value and line endings directly. Property iteration order, underscore/private filtering, aggregate capacity reservation and delimiter remain unchanged. A frozen original formatter compares 9,217 property-map/prefix cases, including empty, Unicode, embedded control/newline characters and a large value. These are byte-preservation cases, not independent chemical-validity claims.

The fifth delta leaves stereo depiction and output construction unchanged, but plain output callers no longer clone and scan loss vectors that they discard. write_mol_with_stereo_report and write_laid_out_mol_into request the original collection explicitly. A frozen original V2000 writer compares plain bytes, reporting bytes and every loss field across all 2,000 pinned SDF records, including the public reporting API. The underlying loss collector is unchanged. This does not skip stereo perception or suppress a returned diagnostic.

## Correctness and local checks

- The final source passes all 1,256 MOL tests across 54 suites, zero failed/ignored, including all six buffer tests. These correctness tests use release opt-level 3, LTO disabled and 16 codegen units in a separate target. The measured source executables use full LTO and one codegen unit; their complete pre-timing snapshots/rewrites are separately compared exactly. Test and timing build boundaries are recorded rather than conflated.
- Differential tests match the original iterator path on 5,796 boundary/resource-limit cases and 3,774 deterministic I/O-fault/interruption cases, including consumed bytes, byte/record counters, termination, exact typed errors and full successful records.
- The two existing unfiltered 1,000-record SDF corpora contribute 20,000 record comparisons across five chunk sizes and both diagnostic modes. Actual BufReader capacities contribute another 10,000 comparisons on the first 1,000 records; 3,340 complete spans use the borrowing path and other spans exercise fallback/refills. These are repeated boundary checks, not 30,000 independent molecules.
- Before timing, both source executables have byte-identical complete internal snapshots and rewrites for all 4,064 input records across the five profiles. SDF and MOL profiles overlap in chemistry; this is an API/format coverage count, not 4,064 independent structures.
- Ten unchanged benchmark-runner gates, release all-target MOL Clippy, workspace formatting and refactor-budget checks pass. An initial test compile compared neighbor iterators directly; the corrected test compares their ordered values. Clippy also rejected an untimed custom enum draft for its large variant; the final helper uses the standard Option/Result return form, with no per-record heap wrapper or lint suppression. Earlier draft checks/failures remain separately retained in the first packet and no timing used them.

The [validation record](2026-10-16-openbabel-sdf-loss-work-artifacts/validation.json), [exact source equality](2026-10-16-openbabel-sdf-loss-work-artifacts/preflight-source-equality.json), [source manifest](2026-10-16-openbabel-sdf-loss-work-artifacts/source-manifest.json) and [reconstructable patch](2026-10-16-openbabel-sdf-loss-work-artifacts/candidate.patch.gz) preserve those claims. The patch was independently reconstructed from its named Git base and every changed-file hash verified. Local checks are separate from hosted CI, main integration and publication.

## Frozen protocol

Baseline: clean main revision 025fa2f63da324869aa144a7520626ae4ebbfd22, source version 1.1.5. The clean integration-baseline executables are reused by exact binary hash from the preceding packet, where they were freshly built from the named commit. The candidate is freshly built. Both use Rust 1.97.0, release opt-level 3, full LTO, one codegen unit, one build job and incremental builds disabled. No timing samples are reused. The previous sparse-ID and aromatic-normalization prototypes are excluded. The unchanged C++ harness executables are reused by exact binary/library hash from the first packet, where they were freshly compiled with Apple Clang -O3 -std=c++17 against Open Babel 3.2.1. No timing samples are reused. RDKit 2026.09.1 is an independent observer, not a timed comparator. All installed Open Babel package libraries/plugins and RDKit package native files are pinned; system and other transitive libraries outside those package directories are excluded.

The two SDF profiles and two diverse V3000 profiles retain all 1,000 rows at offsets 0 and 5,000 of the existing frozen data. The fifth profile retains all 64 controlled linear chains of 16–1,024 atoms. These are exposed development corpora, not an independent holdout. Existing [corpus provenance](../validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json) and [third-party notices](../THIRD_PARTY_NOTICES.md) apply.

Each profile uses 21 sequential rotating three-arm fresh-process blocks and three native batch repetitions per lane. Input loading, process startup and initial plugin discovery are excluded. Parse retains the entire diagnostic batch; SDF framing is inside its parse timer, whereas native MOL blocks are split outside the V3000 timer. Prepared writing follows one untimed warm write. Roundtrip reads a fresh batch and writes it within one timer. The existing Rust SDF writer reuses a String; the native Open Babel writer allocates its API return strings. Whole-process high-water RSS includes loading, preparation, warmup and plugins; it is not a streaming-memory guarantee. Own builds/tests/cache cleanup finish before formal timing. Other host workloads are uncontrolled; every block records load averages.

Ratios above 1 favor the candidate. Paired median-ratio intervals use the unchanged runner's 10,000 deterministic bootstrap resamples and three-decimal reporting. Strict Open Babel wins require semantic prerequisites, CI lower >1 and 21/21 faster blocks. Parse comparator eligibility checks canonical internal graphs; rewrite eligibility checks all seven declared observer fields. Semantic mismatches remain unranked without filtering or repairs. SD-field text appending and requested-only loss collection change in the writer; MOL numeric, stereo depiction and graph serialization remain unchanged. The matrix measures the combined reader/writer source, not an isolated causal decomposition of each delta. Source equality and source-relative speed, equivalent-output comparator wins and toolkit-wide superiority are distinct conclusions.

The [predeclared plan](2026-10-16-openbabel-sdf-loss-work-artifacts/measurement-plan.json) requires exact source snapshots/rewrites, positive parse, prepared-write and roundtrip confidence lower bounds on both SDF slices, and every one of all 15 SDF/V3000 source-relative lane medians >=0.99. The prepared-write improvement condition is added before this candidate measurement; every previous adoption condition remains. Every negative interval is disclosed. There are no adaptive timing reruns, pooled samples, substitutions or relaxed gates.

## Complete timing matrix

| Profile | Operation | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 38.587 | 37.928 | 1.019 [1.011, 1.023]; 19/21 | unranked: semantic mismatch |
| sdf-0 | write | 11.679 | 11.540 | 1.011 [1.004, 1.027]; 16/21 | 6.181 [6.082, 6.224]; 21/21 **strict win** |
| sdf-0 | roundtrip | 50.831 | 49.828 | 1.022 [1.018, 1.028]; 21/21 **strict win** | 7.169 [7.128, 7.182]; 21/21 **strict win** |
| sdf-5000 | parse | 14.306 | 13.946 | 1.029 [1.026, 1.031]; 21/21 **strict win** | 8.920 [8.828, 8.973]; 21/21 **strict win** |
| sdf-5000 | write | 4.811 | 4.719 | 1.017 [1.013, 1.029]; 19/21 | 8.657 [8.594, 8.695]; 21/21 **strict win** |
| sdf-5000 | roundtrip | 19.449 | 18.902 | 1.031 [1.024, 1.035]; 21/21 **strict win** | 8.773 [8.717, 8.802]; 21/21 **strict win** |
| v3000-0 | parse | 39.094 | 39.333 | 0.993 [0.988, 0.996]; 5/21 | unranked: semantic mismatch |
| v3000-0 | write | 8.248 | 8.229 | 1.006 [0.994, 1.011]; 12/21 | unranked: semantic mismatch |
| v3000-0 | roundtrip | 47.703 | 48.040 | 0.995 [0.993, 0.997]; 3/21 | unranked: semantic mismatch |
| v3000-5000 | parse | 15.132 | 15.167 | 0.999 [0.992, 1.002]; 9/21 | unranked: semantic mismatch |
| v3000-5000 | write | 2.123 | 2.125 | 0.999 [0.992, 1.004]; 9/21 | unranked: semantic mismatch |
| v3000-5000 | roundtrip | 17.343 | 17.344 | 0.998 [0.996, 1.002]; 7/21 | unranked: semantic mismatch |
| v3000-chains | parse | 16.950 | 17.049 | 0.990 [0.986, 1.003]; 8/21 | 143.925 [143.001, 145.426]; 21/21 **strict win** |
| v3000-chains | write | 3.767 | 3.792 | 0.999 [0.987, 1.016]; 10/21 | 15.170 [15.019, 15.349]; 21/21 **strict win** |
| v3000-chains | roundtrip | 20.779 | 20.932 | 0.993 [0.989, 0.996]; 4/21 | 119.595 [118.564, 120.510]; 21/21 **strict win** |

## Source adoption decision

```json
{
  "exact_source_snapshots_and_rewrites": true,
  "sdf_primary": {
    "sdf-0": {
      "parse": {
        "median": 1.019,
        "ci95": [
          1.011,
          1.023
        ],
        "winning_blocks": 19,
        "blocks": 21
      },
      "write": {
        "median": 1.011,
        "ci95": [
          1.004,
          1.027
        ],
        "winning_blocks": 16,
        "blocks": 21
      },
      "roundtrip": {
        "median": 1.022,
        "ci95": [
          1.018,
          1.028
        ],
        "winning_blocks": 21,
        "blocks": 21
      }
    },
    "sdf-5000": {
      "parse": {
        "median": 1.029,
        "ci95": [
          1.026,
          1.031
        ],
        "winning_blocks": 21,
        "blocks": 21
      },
      "write": {
        "median": 1.017,
        "ci95": [
          1.013,
          1.029
        ],
        "winning_blocks": 19,
        "blocks": 21
      },
      "roundtrip": {
        "median": 1.031,
        "ci95": [
          1.024,
          1.035
        ],
        "winning_blocks": 21,
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
      "parse": {
        "median": 0.993,
        "ci95": [
          0.988,
          0.996
        ],
        "winning_blocks": 5,
        "blocks": 21
      },
      "roundtrip": {
        "median": 0.995,
        "ci95": [
          0.993,
          0.997
        ],
        "winning_blocks": 3,
        "blocks": 21
      }
    },
    "v3000-5000": {},
    "v3000-chains": {
      "roundtrip": {
        "median": 0.993,
        "ci95": [
          0.989,
          0.996
        ],
        "winning_blocks": 4,
        "blocks": 21
      }
    }
  },
  "local_adoption_gate_passed": true,
  "no_adaptive_timing_reruns": true
}
```

## Whole-process peak RSS

| Profile | Operation | Baseline MiB | Candidate MiB | Open Babel MiB |
|---|---|---:|---:|---:|
| sdf-0 | parse | 15.22 | 15.05 | 52.72 |
| sdf-0 | write | 10.69 | 10.56 | 39.69 |
| sdf-0 | roundtrip | 15.28 | 15.08 | 52.86 |
| sdf-5000 | parse | 11.03 | 10.81 | 41.28 |
| sdf-5000 | write | 7.81 | 7.73 | 33.33 |
| sdf-5000 | roundtrip | 11.03 | 10.86 | 41.38 |
| v3000-0 | parse | 16.53 | 16.47 | 51.92 |
| v3000-0 | write | 10.98 | 10.97 | 38.97 |
| v3000-0 | roundtrip | 16.55 | 16.53 | 52.06 |
| v3000-5000 | parse | 11.11 | 11.08 | 40.88 |
| v3000-5000 | write | 7.86 | 7.80 | 32.92 |
| v3000-5000 | roundtrip | 11.14 | 11.09 | 40.94 |
| v3000-chains | parse | 13.97 | 13.94 | 50.72 |
| v3000-chains | write | 9.44 | 9.44 | 39.33 |
| v3000-chains | roundtrip | 14.02 | 13.91 | 51.45 |

The independent [packet audit](2026-10-16-openbabel-sdf-loss-work-artifacts/packet-audit.json) and [audit script](2026-10-16-openbabel-sdf-loss-work-artifacts/verify_packets.py) reconstruct all five patches and verify all 25 complete matrices against their separate frozen decisions. The first packet passes only at its older baseline; the fifth passes at the named integration baseline. No prior observation is replaced.

## Latest-main correctness integration

The measured source is committed as a1e8c8e87c99a92966547cdf384953063d42ff09
on the named 025fa2f6 baseline. Later main ecda5c53 adds shared 2D/3D E/Z
and stereo-enumeration changes. Integration commit
247b9f8614ebc33d353f9088b3826195c1b34b69 retains all seven optimization-file
hashes, and separately passes 1,258 MOL tests plus 271 perception and 394
SMILES library tests (1,923 passed, zero failed, two existing ignored).
All-target Clippy for those three crates, formatting, refactor budget and
the 481-record index also pass. See the
[integration validation](2026-10-16-openbabel-sdf-loss-work-artifacts/integration-validation.json)
and [post-integration audit](2026-10-16-openbabel-sdf-loss-work-artifacts/integration-packet-audit.json).

These later checks establish correctness integration. The timings above
remain bound to the measured checkpoint; they are not new speed measurements
of the later main merge or published packages. The audit's optional
--require-exact-source flag requires the measured checkpoint rather than a
later checkout.

## Reproduction and limits

Check out the named baseline in an isolated worktree. Decompress/apply the candidate patch for a separate candidate. Build both using CARGO_BUILD_JOBS=1 CARGO_INCREMENTAL=0 cargo build --release --locked -p chematic-mol --example openbabel_sdf_throughput --example openbabel_v3000_throughput; preserve baseline executables before rebuilding. Compile the unchanged C++ harnesses with the recorded pkg-config Open Babel flags. Decompress the five named checked-in corpora. Run the unchanged scripts/bench_openbabel_sdf_throughput.py and scripts/bench_openbabel_v3000_throughput.py with --baseline, --candidate, --openbabel, --corpus, --records <1000-or-64>, --blocks 21 --repeats 3 --source-manifest <manifest.json> --output <report.json> in the pinned observer environment.

The retained launch scripts check complete source/library/binary/corpus hashes before and after each profile. Their 64 MiB reserve only retries an owned checkpoint write with the same in-memory observation on ENOSPC; workers, timers and statistics are unchanged. Darwin time -l needs OS-statistics access. Adjust local paths when reproducing.

Raw complete [SDF 0](2026-10-16-openbabel-sdf-loss-work-artifacts/paired-sdf-0.json.gz), [SDF 5000](2026-10-16-openbabel-sdf-loss-work-artifacts/paired-sdf-5000.json.gz), [MOL 0](2026-10-16-openbabel-sdf-loss-work-artifacts/paired-v3000-0.json.gz), [MOL 5000](2026-10-16-openbabel-sdf-loss-work-artifacts/paired-v3000-5000.json.gz) and [chains](2026-10-16-openbabel-sdf-loss-work-artifacts/paired-v3000-chains.json.gz) matrices contain every block, mismatch, hash and RSS observation. Build/test/format/refactor/Clippy/timing logs and the exact launch scripts are retained here. Full executables and snapshot/rewrite streams remain locally in /private/tmp/chematic-sdf-loss-work-evidence-20261016; their [hashes](2026-10-16-openbabel-sdf-loss-work-artifacts/local-semantic-streams.json) are versioned without duplicating binaries and generated outputs into Git.

Other formats, diverse held-out chemistry, rich query/SGROUP performance, buffers smaller than a record, filesystem/CLI/Python overhead, cold startup, multithreading, additional hosts and published artifacts remain separate comparisons. Existing Open Babel semantic disagreements are not resolved by a faster SDF record boundary. This packet does not validate later main commits or published packages.
