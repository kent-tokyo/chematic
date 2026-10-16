# Requested-only V3000 stereo-loss collection — 2026-10-16

Local adoption gate: **PASS**. Complete Open Babel-wide superiority remains unproved. These are native source measurements on one shared macOS arm64 host; they are not published-package results.

## Implementation and preservation

Baseline: clean source main `57d318c03b5bd22e46b8482e72ba47d190d77701`, source version 1.1.5. Both Rust arms were freshly built from this base and the archived one-file production delta. Earlier unpromoted sparse-ID and descriptor candidates are excluded. The baseline already includes the prior SDF reader/writer improvements.

Plain V3000 output previously called the reporting writer and discarded its `MolStereoLoss`. The candidate keeps the same depiction and serialization, but collects losses only when `write_mol_v3000_with_stereo_report` requests them. Plain output skips unused vector clones, atom scans, sorting and deduplication. The loss collector itself is unchanged. No stereo perception, returned diagnostic, public signature, dependency or unsafe-code change. The 3D and checked writers are unchanged.

The output-construction prefix is byte-identical to the original. A separately frozen original reporting function compares ordinary bytes, reporting bytes and **all** loss fields on 2,064 V3000 records and 2,000 V2000-read graphs, cold and repeated calls. These profiles overlap in structures. Directed cases cover empty graphs, missing/all-zero coordinates, tetrahedral/EZ/square-planar tags, enhanced stereo groups and metadata text; the public reporting API must retain a nonempty square-planar loss list. Existing tests also retain numeric/nonfinite, opaque-property and 3D writer boundaries.

Full MOL validation: 1264 passed, 0 failed, 0 ignored across 55 suites. Ten runner gates, release all-target MOL Clippy, formatting and refactor-budget checks pass. Correctness tests use opt3/LTO=false/16 codegen units; measured executables use full LTO/one codegen unit. Before timing, complete internal snapshots and exact rewrites match on all 4,064 profile input records. These are input/API counts, not 4,064 independent chemical structures.

## Frozen protocol

Two unfiltered 1,000-record SDF slices, two unfiltered 1,000-record V3000 slices at offsets 0 and 5,000, and all 64 controlled chains of 16–1,024 atoms. All are existing exposed development corpora, not a sealed holdout. Existing corpus provenance and third-party notices apply.

Twenty-one rotating three-arm fresh-process blocks per profile, three native batch repetitions per lane. Startup, input loading and initial plugin discovery are outside native timers. SDF framing is inside parsing; V3000 delimiter splitting precedes timing. Parse retains the complete diagnostic batch. Prepared writing follows one untimed warm write. Roundtrip parses a fresh batch and writes within one timer. Rust SDF reuses its String; V3000 and Open Babel use their allocating native string APIs. Native return objects are retained according to the unchanged harness. Whole-process peak RSS includes loading, preparation, warmup and plugins and is not a streaming-memory guarantee. External host activity is uncontrolled and load averages are recorded.

Rust 1.97.0, release opt3/full LTO/one codegen unit/one build job, incremental=false, strip=none. Open Babel 3.2.1 C++ -O3 harnesses are reused by exact binary/library hashes from the prior packet; no timing samples are reused. RDKit 2026.09.1 is an independent observer. The 281 package-local native-file entries are pinned before/after each profile; aliases may repeat content. System/transitive libraries outside those package directories and native data files are outside this seal.

Ratios are reference time / candidate time; >1 favors the candidate. Intervals use the unchanged 10,000-resample paired median bootstrap and three-decimal reporting. Strict Open Babel wins require semantic prerequisites, CI lower>1 and 21/21 faster blocks. Internal canonical graph equality gates parse; all seven declared rewrite fields gate write/roundtrip. Rejected/mismatching rows stay in their original positions and remain unranked, without selection or repair. Output-preserving source speed and equivalent-output comparator wins are separate conclusions.

Local adoption was fixed before timing: exact source snapshots/rewrites and reported losses; prepared V3000 write CI lower>1 on both diverse slices; all 15 source-relative medians >=0.99; disclose every negative interval. No adaptive timing reruns, pooled samples, substituted observations or relaxed gates. Movement in an unchanged lane is not evidence of a new algorithm or its cause.

## Complete matrix

Times are median milliseconds per complete batch (native elapsed divided by the three repetitions). Ratios use paired observations, rather than ratios of independent medians.

| Profile | Mode | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 12.589 | 12.619 | 0.999 [0.997, 1.001]; 10/21 | unranked: semantic prerequisite fails |
| sdf-0 | write | 3.932 | 3.938 | 1.003 [0.993, 1.009]; 13/21 | 6.021 [5.980, 6.077]; 21/21 **strict win** |
| sdf-0 | roundtrip | 16.717 | 16.769 | 0.997 [0.993, 1.005]; 10/21 | 7.090 [7.057, 7.142]; 21/21 **strict win** |
| sdf-5000 | parse | 4.652 | 4.657 | 0.999 [0.994, 1.007]; 10/21 | 8.872 [8.828, 8.914]; 21/21 **strict win** |
| sdf-5000 | write | 1.635 | 1.621 | 1.009 [0.998, 1.025]; 13/21 | 8.406 [8.289, 8.473]; 21/21 **strict win** |
| sdf-5000 | roundtrip | 6.374 | 6.375 | 1.003 [0.993, 1.005]; 11/21 | 8.622 [8.559, 8.674]; 21/21 **strict win** |
| v3000-0 | parse | 13.262 | 13.110 | 1.004 [0.998, 1.011]; 13/21 | unranked: semantic prerequisite fails |
| v3000-0 | write | 2.760 | 2.723 | 1.015 [1.008, 1.021]; 17/21 | unranked: semantic prerequisite fails |
| v3000-0 | roundtrip | 16.243 | 16.112 | 1.005 [0.999, 1.008]; 14/21 | unranked: semantic prerequisite fails |
| v3000-5000 | parse | 5.038 | 5.094 | 0.993 [0.990, 0.998]; 4/21 | unranked: semantic prerequisite fails |
| v3000-5000 | write | 0.707 | 0.695 | 1.015 [1.011, 1.025]; 17/21 | unranked: semantic prerequisite fails |
| v3000-5000 | roundtrip | 5.780 | 5.769 | 1.003 [0.996, 1.005]; 13/21 | unranked: semantic prerequisite fails |
| v3000-chains | parse | 5.657 | 5.623 | 1.006 [1.002, 1.010]; 18/21 | 144.794 [144.198, 145.223]; 21/21 **strict win** |
| v3000-chains | write | 1.276 | 1.215 | 1.047 [1.035, 1.057]; 20/21 | 15.697 [15.504, 15.803]; 21/21 **strict win** |
| v3000-chains | roundtrip | 6.945 | 6.854 | 1.013 [1.012, 1.016]; 19/21 | 121.273 [120.597, 121.877]; 21/21 **strict win** |

## Full adoption decision and every negative interval

```json
{
  "local_adoption_gate_passed": true,
  "all_source_snapshots_and_rewrites_exact": true,
  "primary": {
    "v3000-0": {
      "speedup": {
        "median": 1.015,
        "ci95": [
          1.008,
          1.021
        ],
        "winning_blocks": 17,
        "blocks": 21
      },
      "semantic_prerequisite": true,
      "strict_speed_win": false
    },
    "v3000-5000": {
      "speedup": {
        "median": 1.015,
        "ci95": [
          1.011,
          1.025
        ],
        "winning_blocks": 17,
        "blocks": 21
      },
      "semantic_prerequisite": true,
      "strict_speed_win": false
    }
  },
  "primary_ci_pass": true,
  "all_15_source_median_floor": 0.99,
  "below_floor": {
    "sdf-0": {},
    "sdf-5000": {},
    "v3000-0": {},
    "v3000-5000": {},
    "v3000-chains": {}
  },
  "unranked_source_lanes": [],
  "all_negative_source_intervals": {
    "sdf-0": {},
    "sdf-5000": {},
    "v3000-0": {},
    "v3000-5000": {
      "parse": {
        "median": 0.993,
        "ci95": [
          0.99,
          0.998
        ],
        "winning_blocks": 4,
        "blocks": 21
      }
    },
    "v3000-chains": {}
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
    }
  },
  "no_timing_samples_reused": true,
  "universal_superiority": false
}
```

## Retention and limits

Initial native checks stopped first on an existing SDF test fixture excluded from the sparse checkout, then on ENOSPC during test compilation. The original fixture bytes were restored from Git and the complete successful MOL suite was run before any formal timing. Both failures and successful logs remain. No failed build supplied timing samples.

Only completed owned build/test caches were removed after source/binary/log preservation and process/open-file checks. Unchanged historical benchmark checkout copies were removed through sparse checkout; Git originals and separate raw evidence remain. No broad temporary-directory cleanup or shared active cache deletion was used. The timing wrapper reserves 64MiB and may retry only the same in-memory output write after ENOSPC; it never repeats workers, timers or statistical samples.

The [2026-10-16-openbabel-v3000-loss-work-artifacts](2026-10-16-openbabel-v3000-loss-work-artifacts/README.md) packet retains the complete raw matrices, decisions, source/compiler/corpus/library/binary manifests, frozen plan, patch, preservation proofs, logs and exact launch/verification drivers. Full binaries and full generated snapshots/rewrites remain locally under `/private/tmp/chematic-v3000-loss-work-evidence-20261016` with versioned hashes. Patch reconstruction proves source bytes, not bit-reproducible compiler output.

Additional hosts, held-out chemistry, richer query/SGROUP performance, cold CLI startup, Python/WASM bindings, published artifacts, other formats, force fields and complete Open Babel-wide speed superiority remain outside this packet. Main integration/hosted CI/publication are separate from the named measured source.
