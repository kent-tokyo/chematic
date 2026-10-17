# Open Babel comparison: degree-only stereocenter filtering with steady native batches

Decision: **FAIL: degree-only candidate not adopted**. Both branched-3D parse CI lower bounds must exceed 1; every one of 30 unrounded source median ratios must be >=0.99; source outcomes must be exact. This is source-build evidence on one macOS host, not universal or published-package superiority.

## Change and correctness

Baseline `35b932f911851f51c5cdb89fa371f1c2c4a1b11a`; measured candidate `f7e3c0dfb00ec797f9aabac07b5ef6ac747be00d`. The production change is only the degree>=3 precondition in `stereo_centers`. Every old return shape requires three or four heavy neighbors, including constrained ring nitrogen, and heavy-neighbor count cannot exceed graph degree. Rejecting lower-degree atoms before Vec allocation, H counting and rooted branch traversal preserves every valid-center path and the final shape check. Global `find_sssr` cache preparation remains before the loop. Numerical formatting is byte-identical to the baseline source. No public API, dependency or production unsafe additions.

The complete old classifier is byte-exact after only renaming its declaration; 14,072 parsed graphs, all 6,144 five-vertex/H-override variants and directed boundaries match it. All 4,384 native source snapshots and warmed rewrites match byte for byte. Native validation: 3306 passed, 0 failed, 3 existing ignored across 110 suites covering perception, MOL, chemistry and fingerprints. Release all-target Clippy, formatting and the unchanged refactor ratchet pass; all 24 runner/accounting tests pass. An additional nine-mode native instrumentation preflight verifies raw counts, elapsed time and budgets in all three engines; those observations and the baseline CPU diagnostic are excluded from every formal sample.

## New fixed throughput protocol

The rejected formatter and compound candidates remain in PR #984, with their complete observations. This candidate removes numeric changes and isolates the necessary classifier condition. Before its timing began, the source baseline was updated to latest main including #981/#986. The original prospective plan is retained. No earlier elapsed samples are reused, pooled or used as new observations. The two parser primary lanes and the 0.99 retention threshold remain fixed.

Ten unchanged corpus profiles/hashes: two unfiltered 1,000-row diverse SDF slices, corresponding V3000 slices, 64 linear chains, continued versions, two synthetic 3D chain profiles and two F-branched variants; total 4,384 records. Synthetic chains span 16..1024 carbons with arbitrary deterministic finite nonplanar coordinates. They are not physical conformers, embedding-quality evidence or a representative holdout.

Each profile has 21 rotating three-arm fresh-process blocks, with parse, prepared write and roundtrip. All engines first run the measured mode for at least 500ms, then time at least 100ms and at least three whole batches. Faster engines may process more complete batches. Every sample records raw elapsed, actual iterations, actual records and units; ns per batch is raw elapsed divided by actual iterations, rounded down by less than 1ns. Paired ratios use those per-batch values. All 1,890 sample counts, budgets and normalizations are independently rechecked. This is steady-state throughput, not cold-process or startup latency. Absolute times and whole-process RSS are not directly comparable to prior fixed-three-repeat experiments.

Parse and roundtrip create fresh result objects on every iteration; prepared write reuses its batch. Rust SDF reuses a String; V3000 and Open Babel allocate. Rewrites are checked after the same write warmup/time budget, including all unranked rows. Input loading and native startup/plugin discovery are excluded; V3000 framing is outside timing and SDF framing inside. Whole-process RSS includes preparation, warmup and batch retention; it is not isolated-operation/streaming memory and is not ranked. External host load remains uncontrolled and is logged at each block start.

Both initially empty independent targets use opt3/full LTO/one CGU/one build job/incremental=false/strip=none. All baseline files match Git except seven explicitly listed probe-only instrumentation files, whose bytes match candidate files exactly. Both real production compilations are verified in verbose logs; the actual 3D owner containing the classifier has different machine code. The instrumentation patch and exact sources are retained. Comparator Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27, native C++ probe -O3; observer RDKit 2025.09.3. Sources, baseline files, inputs, binaries, protocol and package-local Open Babel files are sealed before/after profiles; system/package-external runtime data are outside the seal.

Paired bootstrap: 10,000 resamples, seed0 for source, seed1 for Open Babel, intervals reported to three decimals. Raw source medians determine retention. External parse ranking requires accepted equivalent canonical graphs; write/roundtrip also require declared atom/bond fields, coordinates to 1e-4, names/properties and explicit unknown stereo. Rejections/mismatches are retained and unranked. Strict external wins additionally require 21/21 faster blocks. Earlier source ratios are never multiplied into these measurements.

## Complete matrix

Milliseconds per complete batch from actual repetition counts, not total native process duration. Paired medians need not equal ratios of displayed marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 12.246 | 12.250 | 1.001 [0.998, 1.003]; 11/21 | unranked: semantic prerequisite fails |
| sdf-0 | write | 3.705 | 3.783 | 0.980 [0.976, 0.994]; 4/21 | 6.194 [6.158, 6.253]; 21/21 strict win |
| sdf-0 | roundtrip | 16.288 | 16.374 | 0.993 [0.991, 1.002]; 8/21 | 7.147 [7.120, 7.197]; 21/21 strict win |
| sdf-5000 | parse | 4.465 | 4.495 | 1.003 [0.985, 1.007]; 11/21 | 9.162 [9.109, 9.253]; 21/21 strict win |
| sdf-5000 | write | 1.531 | 1.581 | 0.970 [0.957, 0.983]; 2/21 | 8.596 [8.428, 8.675]; 21/21 strict win |
| sdf-5000 | roundtrip | 6.173 | 6.204 | 0.993 [0.983, 1.000]; 7/21 | 8.823 [8.767, 8.860]; 21/21 strict win |
| v3000-0 | parse | 12.917 | 12.900 | 1.001 [0.998, 1.005]; 11/21 | unranked: semantic prerequisite fails |
| v3000-0 | write | 2.691 | 2.670 | 1.004 [1.000, 1.008]; 14/21 | unranked: semantic prerequisite fails |
| v3000-0 | roundtrip | 15.921 | 15.908 | 1.001 [1.000, 1.003]; 13/21 | unranked: semantic prerequisite fails |
| v3000-5000 | parse | 4.879 | 4.914 | 0.998 [0.990, 1.006]; 9/21 | unranked: semantic prerequisite fails |
| v3000-5000 | write | 0.690 | 0.690 | 1.001 [0.999, 1.005]; 12/21 | unranked: semantic prerequisite fails |
| v3000-5000 | roundtrip | 5.583 | 5.640 | 0.992 [0.984, 1.002]; 7/21 | unranked: semantic prerequisite fails |
| v3000-chains | parse | 5.519 | 5.553 | 0.992 [0.986, 1.001]; 8/21 | 150.108 [149.326, 150.821]; 21/21 strict win |
| v3000-chains | write | 1.229 | 1.232 | 0.997 [0.988, 1.002]; 7/21 | 15.321 [15.266, 15.410]; 21/21 strict win |
| v3000-chains | roundtrip | 6.779 | 6.813 | 0.992 [0.983, 1.001]; 7/21 | 124.175 [123.456, 125.232]; 21/21 strict win |
| v3000-continued | parse | 7.649 | 7.719 | 0.991 [0.981, 1.002]; 8/21 | 108.030 [107.143, 109.566]; 21/21 strict win |
| v3000-continued | write | 1.235 | 1.229 | 1.012 [1.002, 1.021]; 15/21 | unranked: semantic prerequisite fails |
| v3000-continued | roundtrip | 8.967 | 8.947 | 1.000 [0.992, 1.007]; 12/21 | unranked: semantic prerequisite fails |
| v3000_3d-spiral | parse | 6.606 | 6.592 | 1.003 [0.997, 1.007]; 12/21 | 126.054 [125.533, 127.033]; 21/21 strict win |
| v3000_3d-spiral | write | 2.482 | 2.444 | 1.017 [1.007, 1.029]; 20/21 | 8.600 [8.571, 8.705]; 21/21 strict win |
| v3000_3d-spiral | roundtrip | 9.177 | 9.070 | 1.011 [1.007, 1.015]; 19/21 | 93.837 [92.895, 93.988]; 21/21 strict win |
| v3000_3d-branched-spiral | parse | 83.273 | 9.907 | 8.424 [8.367, 8.489]; 21/21 | unranked: semantic prerequisite fails |
| v3000_3d-branched-spiral | write | 2.493 | 2.472 | 1.018 [1.006, 1.030]; 17/21 | 8.580 [8.481, 8.740]; 21/21 strict win |
| v3000_3d-branched-spiral | roundtrip | 85.521 | 12.462 | 6.890 [6.854, 6.941]; 21/21 | 115.107 [114.448, 116.056]; 21/21 strict win |
| v3000_3d-offset | parse | 6.823 | 6.830 | 1.001 [0.993, 1.007]; 11/21 | 121.790 [120.782, 123.461]; 21/21 strict win |
| v3000_3d-offset | write | 2.483 | 2.452 | 1.010 [1.003, 1.020]; 16/21 | 8.902 [8.882, 8.928]; 21/21 strict win |
| v3000_3d-offset | roundtrip | 9.337 | 9.306 | 1.006 [1.002, 1.008]; 17/21 | 91.877 [91.427, 92.272]; 21/21 strict win |
| v3000_3d-branched-offset | parse | 83.268 | 10.033 | 8.298 [8.279, 8.347]; 21/21 | 141.134 [140.171, 141.730]; 21/21 strict win |
| v3000_3d-branched-offset | write | 2.500 | 2.432 | 1.020 [1.003, 1.044]; 18/21 | 8.903 [8.769, 8.974]; 21/21 strict win |
| v3000_3d-branched-offset | roundtrip | 86.336 | 12.584 | 6.855 [6.804, 6.887]; 21/21 | 114.425 [113.467, 114.840]; 21/21 strict win |

## Decision and every negative result

```json
{
  "local_adoption_gate_passed": false,
  "primary": {
    "v3000_3d-branched-spiral/parse": {
      "median": 8.424,
      "ci95": [
        8.367,
        8.489
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-branched-offset/parse": {
      "median": 8.298,
      "ci95": [
        8.279,
        8.347
      ],
      "winning_blocks": 21,
      "blocks": 21
    }
  },
  "source_ratios": {
    "sdf-0/parse": {
      "median": 1.001,
      "ci95": [
        0.998,
        1.003
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "sdf-0/write": {
      "median": 0.98,
      "ci95": [
        0.976,
        0.994
      ],
      "winning_blocks": 4,
      "blocks": 21
    },
    "sdf-0/roundtrip": {
      "median": 0.993,
      "ci95": [
        0.991,
        1.002
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "sdf-5000/parse": {
      "median": 1.003,
      "ci95": [
        0.985,
        1.007
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "sdf-5000/write": {
      "median": 0.97,
      "ci95": [
        0.957,
        0.983
      ],
      "winning_blocks": 2,
      "blocks": 21
    },
    "sdf-5000/roundtrip": {
      "median": 0.993,
      "ci95": [
        0.983,
        1.0
      ],
      "winning_blocks": 7,
      "blocks": 21
    },
    "v3000-0/parse": {
      "median": 1.001,
      "ci95": [
        0.998,
        1.005
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "v3000-0/write": {
      "median": 1.004,
      "ci95": [
        1.0,
        1.008
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "v3000-0/roundtrip": {
      "median": 1.001,
      "ci95": [
        1.0,
        1.003
      ],
      "winning_blocks": 13,
      "blocks": 21
    },
    "v3000-5000/parse": {
      "median": 0.998,
      "ci95": [
        0.99,
        1.006
      ],
      "winning_blocks": 9,
      "blocks": 21
    },
    "v3000-5000/write": {
      "median": 1.001,
      "ci95": [
        0.999,
        1.005
      ],
      "winning_blocks": 12,
      "blocks": 21
    },
    "v3000-5000/roundtrip": {
      "median": 0.992,
      "ci95": [
        0.984,
        1.002
      ],
      "winning_blocks": 7,
      "blocks": 21
    },
    "v3000-chains/parse": {
      "median": 0.992,
      "ci95": [
        0.986,
        1.001
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "v3000-chains/write": {
      "median": 0.997,
      "ci95": [
        0.988,
        1.002
      ],
      "winning_blocks": 7,
      "blocks": 21
    },
    "v3000-chains/roundtrip": {
      "median": 0.992,
      "ci95": [
        0.983,
        1.001
      ],
      "winning_blocks": 7,
      "blocks": 21
    },
    "v3000-continued/parse": {
      "median": 0.991,
      "ci95": [
        0.981,
        1.002
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "v3000-continued/write": {
      "median": 1.012,
      "ci95": [
        1.002,
        1.021
      ],
      "winning_blocks": 15,
      "blocks": 21
    },
    "v3000-continued/roundtrip": {
      "median": 1.0,
      "ci95": [
        0.992,
        1.007
      ],
      "winning_blocks": 12,
      "blocks": 21
    },
    "v3000_3d-spiral/parse": {
      "median": 1.003,
      "ci95": [
        0.997,
        1.007
      ],
      "winning_blocks": 12,
      "blocks": 21
    },
    "v3000_3d-spiral/write": {
      "median": 1.017,
      "ci95": [
        1.007,
        1.029
      ],
      "winning_blocks": 20,
      "blocks": 21
    },
    "v3000_3d-spiral/roundtrip": {
      "median": 1.011,
      "ci95": [
        1.007,
        1.015
      ],
      "winning_blocks": 19,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/parse": {
      "median": 8.424,
      "ci95": [
        8.367,
        8.489
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/write": {
      "median": 1.018,
      "ci95": [
        1.006,
        1.03
      ],
      "winning_blocks": 17,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/roundtrip": {
      "median": 6.89,
      "ci95": [
        6.854,
        6.941
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-offset/parse": {
      "median": 1.001,
      "ci95": [
        0.993,
        1.007
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "v3000_3d-offset/write": {
      "median": 1.01,
      "ci95": [
        1.003,
        1.02
      ],
      "winning_blocks": 16,
      "blocks": 21
    },
    "v3000_3d-offset/roundtrip": {
      "median": 1.006,
      "ci95": [
        1.002,
        1.008
      ],
      "winning_blocks": 17,
      "blocks": 21
    },
    "v3000_3d-branched-offset/parse": {
      "median": 8.298,
      "ci95": [
        8.279,
        8.347
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-branched-offset/write": {
      "median": 1.02,
      "ci95": [
        1.003,
        1.044
      ],
      "winning_blocks": 18,
      "blocks": 21
    },
    "v3000_3d-branched-offset/roundtrip": {
      "median": 6.855,
      "ci95": [
        6.804,
        6.887
      ],
      "winning_blocks": 21,
      "blocks": 21
    }
  },
  "negative_intervals": [
    {
      "lane": "sdf-0/write",
      "speedup": {
        "median": 0.98,
        "ci95": [
          0.976,
          0.994
        ],
        "winning_blocks": 4,
        "blocks": 21
      }
    },
    {
      "lane": "sdf-5000/write",
      "speedup": {
        "median": 0.97,
        "ci95": [
          0.957,
          0.983
        ],
        "winning_blocks": 2,
        "blocks": 21
      }
    }
  ],
  "source_floor_failures": [
    {
      "lane": "sdf-0/write",
      "speedup": {
        "median": 0.98,
        "ci95": [
          0.976,
          0.994
        ],
        "winning_blocks": 4,
        "blocks": 21
      },
      "unrounded_median": 0.9802092730095815
    },
    {
      "lane": "sdf-5000/write",
      "speedup": {
        "median": 0.97,
        "ci95": [
          0.957,
          0.983
        ],
        "winning_blocks": 2,
        "blocks": 21
      },
      "unrounded_median": 0.9700774957939933
    }
  ],
  "all_source_snapshots_and_rewrites_exact": true,
  "source_valid_epochs": 10,
  "universal_superiority": false
}
```

## Same-binary diagnostic, excluded from adoption

After all ten formal candidate epochs completed, both Rust labels execute the exact same `baseline-sdf` path and bytes. The same Open Babel third arm, unfiltered SDF-0 corpus, modes, warmup and time budgets are retained. The two designs were fixed in `control-plan.json` before any control samples: 21 cyclic rotations and 24 blocks spanning all six permutations, each used four times with one fixed seed. No control is repeated and no control observations enter the candidate ratios, retention gates or external rankings.

A nonunit ratio here indicates possible order/environment effects for identical executable bytes. A ratio consistent with one does not prove the candidate slowdown arose from code layout. The failed candidate remains unadopted regardless of these diagnostic observations.

| Design | Mode | First label / second label median | 95% CI | Second label wins |
|---|---|---:|---|---:|
| cyclic | parse | 1.001 | 0.997 .. 1.004 | 13/21 |
| cyclic | write | 1.016 | 1.002 .. 1.032 | 16/21 |
| cyclic | roundtrip | 1.004 | 1.000 .. 1.019 | 14/21 |
| balanced | parse | 1.001 | 0.994 .. 1.004 | 12/24 |
| balanced | write | 0.996 | 0.991 .. 1.010 | 9/24 |
| balanced | roundtrip | 0.996 | 0.995 .. 1.002 | 9/24 |

## Evidence and reproduction

[Complete packet](2026-10-17-openbabel-degree-steady-artifacts/README.md) includes raw samples/counts/budgets, complete input and outcome streams, exact common probe patch, source/binary/native seals, compiler logs, disassembly, frozen-classifier proof, native test logs and recomputed statistical audit. Rebuild the two named source revisions in independent roots, apply the same retained common probe patch to the baseline, run `build.py`, `native.py`, `instrument.py`, `profile.py`, `seal.py`, then `formal.py`; adapt the original absolute paths. `assemble.py` verifies the raw accounting and all statistics. Full local evidence: `/private/tmp/chematic-degree-only-steady-evidence-20261017`. A later integration revision requires correctness validation and is not the timed candidate.
