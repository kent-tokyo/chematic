# Interrupted balanced throughput: outlined stereocenter classification

Decision: **INCOMPLETE AND FAILED: outlined candidate not adopted**. Baseline `35b932f911851f51c5cdb89fa371f1c2c4a1b11a`; measured source `a774a05c55d038aa4abf3122d09eaa87b604ffaa`. One macOS host, native source builds, steady-state throughput. This does not establish universal superiority, cold-process latency, physical-conformer quality or published-wheel performance.

## Interruption and missing data

The original run stopped with transient ENOSPC while opening the sixth block result for `v3000-continued`. The original JSON remains readable with five blocks. Six completed native calls in the failed block have no persisted timing values; they are not reconstructed. Of 1,584 planned native rows, 1,470 are retained, six failed persistence and 108 were never measured. Nine profiles have complete 24-block epochs; the three continued-profile lanes remain unranked. This experiment is ineligible for adoption even apart from the failed SDF floors.

After checking process/open-file use, unused owned build caches were removed with all copied formal binaries, source, inputs, logs and data preserved. The exact allocator/process causing the transient allocation shortage was not identified. Only the four profiles never reached were then measured for the first time with the original timing-loop bytes and fixed orders, while all six preexisting raw files were SHA-checked unchanged. The completed profiles and interrupted profile were never rerun. No observations within a profile were pooled across this interruption. The original full-plan assembler is retained; a separate interrupted-report assembler checks the available coverage and forces the overall decision false.

## Source and correctness

The production change is the necessary degree>=3 precondition in `stereo_centers`, plus `inline(never)` on that classifier. Every shape emitted by the encoded baseline classifier needs three or four heavy neighbors, whose count cannot exceed graph degree. The guard preserves global ring-cache preparation, valid-center branches and final shape tests. Numeric formatting is unchanged. No API, dependency or production unsafe additions.

The full baseline classifier is frozen byte-exact after renaming: 14,072 pinned parsed graphs, 6,144 five-atom/H variants and directed boundaries match. All 4,384 fresh native parse snapshots and warmed rewrites are byte-identical across source arms and match the preceding C baseline streams. Fresh SDF observer assessment also reproduces all accepted rows and mismatch counts; semantic failures remain unranked. These establish preservation of existing behavior within the stated tests, not chemical correctness of all possible inputs.

Fresh native validation: 3306 passed, 0 failed, 3 existing ignored in 110 suites (perception, MOL, chemistry, fingerprints). All-target release Clippy, formatting, unchanged refactor ratchet and 24 runner/accounting tests pass. The nine-mode instrumentation check is excluded from formal timings. Three preparation failures stopped before the first formal sample (an import error, an inappropriate comparison to Open Babel's time-bearing program header, and a shadowed observer variable). All original drivers, logs, seals and unmeasured outputs are retained. A standalone ten-profile preflight then passed before any timing epoch. Open Babel SDF semantic prerequisites are recomputed independently from fresh output; its program-header timestamp is outside the unchanged declared molecular fields. Neither source byte-exact conditions nor performance thresholds/orders were changed.

Actual generated code calls a separate classifier from the 3D reader. The new reader owner has 3628 instructions and the separate classifier 274; the prior inlined baseline owner had 3,881 instructions. Instruction counts explain source realization, not total performance or the cause of previous writer penalties. Full disassembly is retained.

## Prospective protocol and limits

The formatter-only and compound attempts in #984 and the degree-only C attempt in [its full report](2026-10-17-openbabel-degree-steady.md) failed their fixed gates. C retains its SDF write penalties and all negative observations. Its same-binary controls showed a nonunit cyclic write interval, while all balanced-control intervals included one. Those controls do not establish causation and do not erase C's failed result. This is a distinct outlined source candidate with a newly fixed protocol, not an unchanged-candidate rerun or pooled analysis.

Before any D timing: 24 blocks/profile; all six three-arm permutations used four times for the two SDF profiles, and twelve AB/twelve BA source orders for the remaining eight profiles, with one fixed shuffle per profile. Planned-order prefixes are audited for every saved row; complete-profile position balance is audited. The interrupted five-block prefix is not ranked. Source adoption requires both predeclared branched-3D parse CI lower bounds>1 and every one of 30 unrounded paired source medians>=0.99, plus exact outputs. No threshold exceptions, filtering, performance-informed reruns, retiming or reuse of previous elapsed samples.

The same ten complete corpus hashes total 4,384 records: two 1,000-row SDF slices, corresponding V3000 slices, 64 linear chains, continued variants and four arbitrary-coordinate 3D chain/F-branched variants spanning 16..1024 carbons. These synthetic coordinates are deterministic nonplanar inputs, not physical conformers or a representative holdout. Source-only timing was planned for all ten, but only nine epochs completed; fresh external Open Babel timings cover only the two SDF slices. There are no new external speed rankings for the other eight profiles, and old speed ratios are never multiplied.

Each mode receives >=500ms warmup, then >=100ms timed work and >=3 complete batches. Actual repeats can differ by engine; raw elapsed, repeats, processed records, output units and budgets remain in every row. Normalization is raw elapsed/actual repeats, rounded down by <1ns per batch. All 1,470 persisted rows are independently audited; six unavailable values are not reconstructed and 108 planned rows were not measured. Parse/roundtrip create fresh objects and retain the full batch; prepared write reuses its batch. SDF Rust reuses a String; V3000 and Open Babel allocate. Input loading and startup/plugin discovery are excluded; SDF framing is inside timing and V3000 framing outside. Whole-process RSS includes preparation/warmup/batch retention, is not isolated-operation or streaming memory, and is not ranked. Host load is recorded but external activity and CPU frequency are uncontrolled.

The baseline reuses exact C sealed binaries, actual empty-target production compiler logs and source proof at the then-current main snapshot `35b932f911851f51c5cdb89fa371f1c2c4a1b11a`; no C timing samples are used. Baseline source and every copied binary still match that seal. The candidate uses a newly empty target with opt3/full LTO/one CGU/one build job/incremental=false/strip=none. Both production libraries are recompiled; seven common native probe sources are byte-identical in the two arms. Compiler versions are checked against C; package-local Open Babel files, sources, protocol, inputs and binaries are sealed before/after profiles. System/package-external runtime data are outside the seal.

Open Babel `Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27`, native C++ -O3; observer RDKit `2025.09.3`. External parse ranking requires accepted equivalent canonical graphs. Write/roundtrip additionally require declared atom/bond fields, coordinates to 1e-4, names/properties and explicit unknown stereo. Strict external wins require CI lower>1 and 24/24 faster blocks. Paired bootstrap uses 10,000 resamples, seed0 source/seed1 Open Babel; intervals display three decimals, adoption uses unrounded source medians.

## Complete coverage matrix, including the incomplete profile

Milliseconds per whole unfiltered batch, not process duration. Paired medians need not equal ratios of marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Fresh Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 12.264 | 12.254 | 1.000 [0.996, 1.005]; 13/24 | unranked: semantic prerequisite fails |
| sdf-0 | write | 3.708 | 3.794 | 0.980 [0.970, 0.989]; 3/24 | 6.187 [6.118, 6.233]; 24/24 strict win |
| sdf-0 | roundtrip | 16.335 | 16.314 | 0.998 [0.993, 1.004]; 11/24 | 7.161 [7.134, 7.206]; 24/24 strict win |
| sdf-5000 | parse | 4.461 | 4.481 | 0.991 [0.986, 1.002]; 8/24 | 9.198 [9.145, 9.236]; 24/24 strict win |
| sdf-5000 | write | 1.525 | 1.574 | 0.969 [0.961, 0.978]; 1/24 | 8.599 [8.502, 8.700]; 24/24 strict win |
| sdf-5000 | roundtrip | 6.137 | 6.186 | 0.992 [0.984, 0.997]; 6/24 | 8.888 [8.825, 8.961]; 24/24 strict win |
| v3000-0 | parse | 12.937 | 12.839 | 1.009 [1.004, 1.010]; 20/24 | not measured in this candidate |
| v3000-0 | write | 2.700 | 2.697 | 1.003 [0.990, 1.008]; 15/24 | not measured in this candidate |
| v3000-0 | roundtrip | 15.946 | 15.857 | 1.004 [1.002, 1.009]; 20/24 | not measured in this candidate |
| v3000-5000 | parse | 4.866 | 4.910 | 0.992 [0.987, 1.001]; 9/24 | not measured in this candidate |
| v3000-5000 | write | 0.692 | 0.692 | 0.998 [0.996, 1.003]; 8/24 | not measured in this candidate |
| v3000-5000 | roundtrip | 5.629 | 5.614 | 1.003 [0.992, 1.007]; 14/24 | not measured in this candidate |
| v3000-chains | parse | 5.511 | 5.541 | 0.995 [0.990, 1.005]; 10/24 | not measured in this candidate |
| v3000-chains | write | 1.233 | 1.245 | 0.993 [0.986, 0.998]; 6/24 | not measured in this candidate |
| v3000-chains | roundtrip | 6.782 | 6.808 | 0.998 [0.991, 1.006]; 10/24 | not measured in this candidate |
| v3000-continued | parse | — | — | incomplete 5/24; no ratio or CI | not measured |
| v3000-continued | write | — | — | incomplete 5/24; no ratio or CI | not measured |
| v3000-continued | roundtrip | — | — | incomplete 5/24; no ratio or CI | not measured |
| v3000_3d-spiral | parse | 6.620 | 6.584 | 1.004 [1.002, 1.008]; 17/24 | not measured in this candidate |
| v3000_3d-spiral | write | 2.482 | 2.468 | 1.008 [0.996, 1.019]; 15/24 | not measured in this candidate |
| v3000_3d-spiral | roundtrip | 9.166 | 9.071 | 1.012 [1.006, 1.015]; 21/24 | not measured in this candidate |
| v3000_3d-branched-spiral | parse | 83.039 | 9.927 | 8.364 [8.281, 8.422]; 24/24 | not measured in this candidate |
| v3000_3d-branched-spiral | write | 2.496 | 2.449 | 1.016 [1.003, 1.025]; 17/24 | not measured in this candidate |
| v3000_3d-branched-spiral | roundtrip | 85.475 | 12.418 | 6.890 [6.840, 6.938]; 24/24 | not measured in this candidate |
| v3000_3d-offset | parse | 6.733 | 6.724 | 1.004 [0.996, 1.009]; 14/24 | not measured in this candidate |
| v3000_3d-offset | write | 2.472 | 2.412 | 1.024 [1.016, 1.034]; 20/24 | not measured in this candidate |
| v3000_3d-offset | roundtrip | 9.277 | 9.150 | 1.013 [1.006, 1.018]; 18/24 | not measured in this candidate |
| v3000_3d-branched-offset | parse | 83.270 | 10.012 | 8.345 [8.257, 8.366]; 24/24 | not measured in this candidate |
| v3000_3d-branched-offset | write | 2.478 | 2.420 | 1.030 [1.010, 1.047]; 21/24 | not measured in this candidate |
| v3000_3d-branched-offset | roundtrip | 85.972 | 12.531 | 6.850 [6.815, 6.883]; 24/24 | not measured in this candidate |

## Decision and all negative observations

```json
{
  "experiment_complete": false,
  "local_adoption_gate_passed": false,
  "reasons": [
    "Two completed SDF source write medians below fixed0.99 floor",
    "v3000-continued stopped during save; three source lanes lack a complete24-block epoch"
  ],
  "primary": {
    "v3000_3d-branched-spiral/parse": {
      "median": 8.364,
      "ci95": [
        8.281,
        8.422
      ],
      "winning_blocks": 24,
      "blocks": 24
    },
    "v3000_3d-branched-offset/parse": {
      "median": 8.345,
      "ci95": [
        8.257,
        8.366
      ],
      "winning_blocks": 24,
      "blocks": 24
    }
  },
  "source_ratios": {
    "sdf-0/parse": {
      "median": 1.0,
      "ci95": [
        0.996,
        1.005
      ],
      "winning_blocks": 13,
      "blocks": 24
    },
    "sdf-0/write": {
      "median": 0.98,
      "ci95": [
        0.97,
        0.989
      ],
      "winning_blocks": 3,
      "blocks": 24
    },
    "sdf-0/roundtrip": {
      "median": 0.998,
      "ci95": [
        0.993,
        1.004
      ],
      "winning_blocks": 11,
      "blocks": 24
    },
    "sdf-5000/parse": {
      "median": 0.991,
      "ci95": [
        0.986,
        1.002
      ],
      "winning_blocks": 8,
      "blocks": 24
    },
    "sdf-5000/write": {
      "median": 0.969,
      "ci95": [
        0.961,
        0.978
      ],
      "winning_blocks": 1,
      "blocks": 24
    },
    "sdf-5000/roundtrip": {
      "median": 0.992,
      "ci95": [
        0.984,
        0.997
      ],
      "winning_blocks": 6,
      "blocks": 24
    },
    "v3000-0/parse": {
      "median": 1.009,
      "ci95": [
        1.004,
        1.01
      ],
      "winning_blocks": 20,
      "blocks": 24
    },
    "v3000-0/write": {
      "median": 1.003,
      "ci95": [
        0.99,
        1.008
      ],
      "winning_blocks": 15,
      "blocks": 24
    },
    "v3000-0/roundtrip": {
      "median": 1.004,
      "ci95": [
        1.002,
        1.009
      ],
      "winning_blocks": 20,
      "blocks": 24
    },
    "v3000-5000/parse": {
      "median": 0.992,
      "ci95": [
        0.987,
        1.001
      ],
      "winning_blocks": 9,
      "blocks": 24
    },
    "v3000-5000/write": {
      "median": 0.998,
      "ci95": [
        0.996,
        1.003
      ],
      "winning_blocks": 8,
      "blocks": 24
    },
    "v3000-5000/roundtrip": {
      "median": 1.003,
      "ci95": [
        0.992,
        1.007
      ],
      "winning_blocks": 14,
      "blocks": 24
    },
    "v3000-chains/parse": {
      "median": 0.995,
      "ci95": [
        0.99,
        1.005
      ],
      "winning_blocks": 10,
      "blocks": 24
    },
    "v3000-chains/write": {
      "median": 0.993,
      "ci95": [
        0.986,
        0.998
      ],
      "winning_blocks": 6,
      "blocks": 24
    },
    "v3000-chains/roundtrip": {
      "median": 0.998,
      "ci95": [
        0.991,
        1.006
      ],
      "winning_blocks": 10,
      "blocks": 24
    },
    "v3000_3d-spiral/parse": {
      "median": 1.004,
      "ci95": [
        1.002,
        1.008
      ],
      "winning_blocks": 17,
      "blocks": 24
    },
    "v3000_3d-spiral/write": {
      "median": 1.008,
      "ci95": [
        0.996,
        1.019
      ],
      "winning_blocks": 15,
      "blocks": 24
    },
    "v3000_3d-spiral/roundtrip": {
      "median": 1.012,
      "ci95": [
        1.006,
        1.015
      ],
      "winning_blocks": 21,
      "blocks": 24
    },
    "v3000_3d-branched-spiral/parse": {
      "median": 8.364,
      "ci95": [
        8.281,
        8.422
      ],
      "winning_blocks": 24,
      "blocks": 24
    },
    "v3000_3d-branched-spiral/write": {
      "median": 1.016,
      "ci95": [
        1.003,
        1.025
      ],
      "winning_blocks": 17,
      "blocks": 24
    },
    "v3000_3d-branched-spiral/roundtrip": {
      "median": 6.89,
      "ci95": [
        6.84,
        6.938
      ],
      "winning_blocks": 24,
      "blocks": 24
    },
    "v3000_3d-offset/parse": {
      "median": 1.004,
      "ci95": [
        0.996,
        1.009
      ],
      "winning_blocks": 14,
      "blocks": 24
    },
    "v3000_3d-offset/write": {
      "median": 1.024,
      "ci95": [
        1.016,
        1.034
      ],
      "winning_blocks": 20,
      "blocks": 24
    },
    "v3000_3d-offset/roundtrip": {
      "median": 1.013,
      "ci95": [
        1.006,
        1.018
      ],
      "winning_blocks": 18,
      "blocks": 24
    },
    "v3000_3d-branched-offset/parse": {
      "median": 8.345,
      "ci95": [
        8.257,
        8.366
      ],
      "winning_blocks": 24,
      "blocks": 24
    },
    "v3000_3d-branched-offset/write": {
      "median": 1.03,
      "ci95": [
        1.01,
        1.047
      ],
      "winning_blocks": 21,
      "blocks": 24
    },
    "v3000_3d-branched-offset/roundtrip": {
      "median": 6.85,
      "ci95": [
        6.815,
        6.883
      ],
      "winning_blocks": 24,
      "blocks": 24
    }
  },
  "negative_intervals": [
    {
      "lane": "sdf-0/write",
      "speedup": {
        "median": 0.98,
        "ci95": [
          0.97,
          0.989
        ],
        "winning_blocks": 3,
        "blocks": 24
      }
    },
    {
      "lane": "sdf-5000/write",
      "speedup": {
        "median": 0.969,
        "ci95": [
          0.961,
          0.978
        ],
        "winning_blocks": 1,
        "blocks": 24
      }
    },
    {
      "lane": "sdf-5000/roundtrip",
      "speedup": {
        "median": 0.992,
        "ci95": [
          0.984,
          0.997
        ],
        "winning_blocks": 6,
        "blocks": 24
      }
    },
    {
      "lane": "v3000-chains/write",
      "speedup": {
        "median": 0.993,
        "ci95": [
          0.986,
          0.998
        ],
        "winning_blocks": 6,
        "blocks": 24
      }
    }
  ],
  "source_floor_failures": [
    {
      "lane": "sdf-0/write",
      "speedup": {
        "median": 0.98,
        "ci95": [
          0.97,
          0.989
        ],
        "winning_blocks": 3,
        "blocks": 24
      },
      "unrounded_median": 0.979881088853413
    },
    {
      "lane": "sdf-5000/write",
      "speedup": {
        "median": 0.969,
        "ci95": [
          0.961,
          0.978
        ],
        "winning_blocks": 1,
        "blocks": 24
      },
      "unrounded_median": 0.9694023385136397
    }
  ],
  "missing_lanes": [
    "v3000-continued/parse",
    "v3000-continued/write",
    "v3000-continued/roundtrip"
  ],
  "completed_profiles": [
    "sdf-0",
    "sdf-5000",
    "v3000-0",
    "v3000-5000",
    "v3000-chains",
    "v3000_3d-spiral",
    "v3000_3d-branched-spiral",
    "v3000_3d-offset",
    "v3000_3d-branched-offset"
  ],
  "incomplete_profiles": [
    {
      "profile": "v3000-continued",
      "persisted_blocks": 5,
      "planned_blocks": 24,
      "reason": "ENOSPC during persistence of sixth block; original five blocks retained byte-exact",
      "unpersisted_native_rows": 6
    }
  ],
  "source_valid_epochs": 9,
  "planned_native_rows": 1584,
  "persisted_native_rows": 1470,
  "unpersisted_native_rows": 6,
  "never_measured_native_rows": 108,
  "fresh_SDF_Open_Babel": {
    "sdf-0/parse": {
      "semantic_prerequisite": false,
      "speedup": null,
      "strict_speed_win": false
    },
    "sdf-0/write": {
      "semantic_prerequisite": true,
      "speedup": {
        "median": 6.187,
        "ci95": [
          6.118,
          6.233
        ],
        "winning_blocks": 24,
        "blocks": 24
      },
      "strict_speed_win": true
    },
    "sdf-0/roundtrip": {
      "semantic_prerequisite": true,
      "speedup": {
        "median": 7.161,
        "ci95": [
          7.134,
          7.206
        ],
        "winning_blocks": 24,
        "blocks": 24
      },
      "strict_speed_win": true
    },
    "sdf-5000/parse": {
      "semantic_prerequisite": true,
      "speedup": {
        "median": 9.198,
        "ci95": [
          9.145,
          9.236
        ],
        "winning_blocks": 24,
        "blocks": 24
      },
      "strict_speed_win": true
    },
    "sdf-5000/write": {
      "semantic_prerequisite": true,
      "speedup": {
        "median": 8.599,
        "ci95": [
          8.502,
          8.7
        ],
        "winning_blocks": 24,
        "blocks": 24
      },
      "strict_speed_win": true
    },
    "sdf-5000/roundtrip": {
      "semantic_prerequisite": true,
      "speedup": {
        "median": 8.888,
        "ci95": [
          8.825,
          8.961
        ],
        "winning_blocks": 24,
        "blocks": 24
      },
      "strict_speed_win": true
    }
  },
  "all_4384_preflight_source_snapshots_and_rewrites_exact": true,
  "interrupted_profile_not_rerun": true,
  "four_remaining_profiles_measured_for_first_time_with_identical_original_timing_loop": true,
  "universal_superiority": false
}
```

## Evidence and reproduction

[Hashed evidence packet](2026-10-17-openbabel-degree-outline-artifacts/README.md) contains every persisted raw native sample, counts/budgets, input and output streams, exact source patch/common instrumentation, reused baseline provenance, fresh candidate code, native test/compiler logs, full disassembly and the independently recomputed audit. Full local evidence: `/private/tmp/chematic-degree-outline-evidence-20261017`. Rebuild baseline main with the retained common probe patch and the measured candidate in independent targets with recorded flags; run `build_native.py`, `instrument.py`, `seal.py`, `formal.py` after adapting original paths. The recorded interrupted run additionally uses `resume_remaining.py` for only the four previously unmeasured profiles and `prepare_interrupted_report.py` / `assemble_interrupted.py` to report the incomplete result; the original full-plan `assemble.py` cannot validate this incomplete run. The fresh process after native validation must be free of own competing builds/profiling. Retain all failures and do not pool previous elapsed data. Hosted CI completion, integration revisions and main merge are separate from this measured source.
