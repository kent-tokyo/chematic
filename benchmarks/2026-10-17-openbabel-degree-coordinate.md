# Open Babel comparison: coordinate pairs and impossible-center degree filtering

Decision: **FAIL: compound candidate not adopted**. Both diverse SDF write lanes and both branched-3D parse lanes must improve with paired CI lower bounds above 1. Every one of the 30 unrounded source median ratios must be >=0.99 and all source outputs must match. This is one macOS source experiment, not universal or published-package superiority.

## Change and correctness

Baseline `34f141266a81264f42990f98acbe25c91d9074af`; measured compound candidate `086a3da876de812ac01d6784911beb42e67e9c80`. Both contain the identical independent ring-partner repair #986 over main `f20304f4`. The formatter converts decimal digits in pairs with a 200-byte ASCII table and combines common padding into the same stack buffer. Rounding, negative zero and general-format fallbacks are retained.

The classifier now rejects graph degree <3 before allocating heavy-neighbor vectors, calculating H counts or traversing branch signatures. Every former return shape requires three or four heavy neighbors, including constrained ring nitrogen, and heavy-neighbor count cannot exceed graph degree. The original global `find_sssr` preparation still runs first; all valid-center paths and the final output-shape check are unchanged. No public API, dependency or unsafe additions. The complete frozen old classifier matches the base source exactly after declaration renaming and is compared on 14,072 pinned graphs, all 6,144 five-vertex/H-override graphs and directed boundaries. The frozen numeric formatter and Rust standard formatting remain independent byte oracles for every fractional pair, 100,000 deterministic random f64 bits, padding/sign/rounding boundaries and the pre-existing 400,000 generated coordinate cases.

All 4,384 source native snapshots and rewritten records are byte-exact. Native validation: 3307 passed, 0 failed, 3 existing ignored across 110 suites covering perception, MOL, chemistry and fingerprints. Release all-target Clippy, formatting, unchanged refactor budget and ten runner tests pass. A separate baseline branched-3D parse CPU trace retains classifier/branch-traversal ancestry; its worker elapsed value is excluded from all formal samples.

## Separate candidate and fixed protocol

The [formatter-only experiment](2026-10-17-openbabel-coordinate-pairs.md) improved SDF writes but failed its retention gate at branched-spiral 3D parse, median 0.983815. This compound candidate adds a necessary precondition to remove that workload's known impossible-center traversal. It has a new fixed plan and a separate timing epoch. No old elapsed observations are reused, pooled or filtered, and no retention threshold is relaxed. The previous candidate and its initial pre-measurement main-regression failure remain fully recorded in that earlier packet.

Ten unchanged input profiles/hashes: two unfiltered 1,000-row diverse SDF slices; corresponding V3000 slices; 64 linear chains; continued lines; two synthetic 3D chain profiles and two fluorine-branched versions. Synthetic chains span 16..1024 carbons and use arbitrary deterministic finite nonplanar coordinates. They are not physical conformers, embedding-quality evidence or a representative holdout. Four predeclared primary lanes and all 30 retention lanes are recorded below, including every negative interval and unranked arm.

Twenty-one rotating three-arm fresh-process blocks per profile, three native batch repeats, parse/prepared write/roundtrip. Native startup, input loading and plugin discovery are excluded. V3000 framing is outside the timer and SDF framing inside. Diagnostic readers retain full batches. Prepared writing has an untimed warm write; Rust SDF reuses a String while V3000 and Open Babel allocate. Whole-process high-water RSS includes preparation, retained batches and warmup, not isolated-operation or streaming memory. External host load is uncontrolled and logged at each block start.

The exact independent baseline binaries, source tree and verbose build logs are reused from the sealed formatter-only experiment; every source file still matches Git `34f141266a81264f42990f98acbe25c91d9074af` and every binary hash matches that prior seal. The candidate uses a new initially empty independent target. Rust compiler identity is verified unchanged; opt3/full LTO/one CGU/one build job/incremental=false/strip=none. Both production libraries are verified in the verbose logs. Actual coordinate-formatter bodies and the 3D perception owner containing the classifier have different machine code in the arms. `baseline-reuse.json` records the original baseline build command and source root; all timing blocks are new. Native comparator: Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27, C++ probe -O3; independent observer RDKit 2025.09.3. Seals cover tracked source, baseline files, corpora, binaries and package-local Open Babel files before/after profiles. System/package-external runtime data remain outside the seal.

Ratios are reference/candidate time, deterministic paired 10,000-resample bootstrap, seed 0 for source and 1 for Open Babel, reported to three decimals. Retention uses unrounded medians. Parse ranking requires accepted equivalent canonical graphs; write/roundtrip ranking also requires declared atom/bond fields, coordinates to 1e-4, names/properties and explicit unknown stereo. Mismatches/rejections remain unranked. Strict external wins additionally require 21/21 faster blocks. Ratios from the predecessor or older experiments are never multiplied into these measurements.

## Complete matrix

Milliseconds per full batch: native elapsed divided by three. Paired medians need not equal ratios of displayed marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 12.592 | 12.557 | 1.003 [0.995, 1.009]; 13/21 | unranked: semantic prerequisite fails |
| sdf-0 | write | 3.875 | 3.856 | 1.002 [0.993, 1.016]; 11/21 | 6.192 [6.118, 6.231]; 21/21 strict win |
| sdf-0 | roundtrip | 16.575 | 16.594 | 0.998 [0.995, 1.003]; 10/21 | 7.123 [7.084, 7.187]; 21/21 strict win |
| sdf-5000 | parse | 4.970 | 4.860 | 1.011 [0.995, 1.018]; 14/21 | 9.343 [9.096, 9.439]; 21/21 strict win |
| sdf-5000 | write | 1.661 | 1.674 | 0.985 [0.943, 1.022]; 8/21 | 8.919 [8.220, 9.462]; 21/21 strict win |
| sdf-5000 | roundtrip | 6.630 | 6.590 | 0.997 [0.986, 1.007]; 9/21 | 9.029 [8.826, 9.259]; 21/21 strict win |
| v3000-0 | parse | 13.268 | 13.342 | 0.996 [0.981, 1.003]; 8/21 | unranked: semantic prerequisite fails |
| v3000-0 | write | 2.821 | 2.825 | 1.002 [0.987, 1.014]; 12/21 | unranked: semantic prerequisite fails |
| v3000-0 | roundtrip | 16.500 | 16.323 | 1.005 [0.997, 1.024]; 14/21 | unranked: semantic prerequisite fails |
| v3000-5000 | parse | 5.139 | 5.149 | 1.003 [0.998, 1.008]; 13/21 | unranked: semantic prerequisite fails |
| v3000-5000 | write | 0.717 | 0.695 | 1.034 [1.027, 1.042]; 20/21 | unranked: semantic prerequisite fails |
| v3000-5000 | roundtrip | 5.887 | 5.859 | 1.008 [1.002, 1.020]; 15/21 | unranked: semantic prerequisite fails |
| v3000-chains | parse | 5.738 | 5.701 | 1.013 [0.999, 1.020]; 14/21 | 143.675 [142.970, 144.257]; 21/21 strict win |
| v3000-chains | write | 1.256 | 1.249 | 1.001 [0.989, 1.023]; 11/21 | 15.212 [15.103, 15.324]; 21/21 strict win |
| v3000-chains | roundtrip | 6.978 | 6.966 | 1.002 [0.993, 1.007]; 12/21 | 120.580 [119.065, 121.629]; 21/21 strict win |
| v3000-continued | parse | 7.991 | 7.920 | 1.012 [1.000, 1.023]; 14/21 | 104.837 [104.227, 105.789]; 21/21 strict win |
| v3000-continued | write | 1.256 | 1.253 | 1.008 [0.982, 1.032]; 11/21 | unranked: semantic prerequisite fails |
| v3000-continued | roundtrip | 9.239 | 9.163 | 1.012 [0.999, 1.016]; 13/21 | unranked: semantic prerequisite fails |
| v3000_3d-spiral | parse | 6.814 | 6.764 | 1.007 [1.006, 1.009]; 17/21 | 121.113 [120.703, 121.950]; 21/21 strict win |
| v3000_3d-spiral | write | 2.502 | 2.484 | 1.006 [0.998, 1.015]; 13/21 | 8.485 [8.406, 8.563]; 21/21 strict win |
| v3000_3d-spiral | roundtrip | 9.315 | 9.230 | 1.005 [1.000, 1.012]; 14/21 | 90.630 [90.167, 91.051]; 21/21 strict win |
| v3000_3d-branched-spiral | parse | 82.590 | 9.974 | 8.268 [8.194, 8.341]; 21/21 | unranked: semantic prerequisite fails |
| v3000_3d-branched-spiral | write | 2.490 | 2.496 | 0.992 [0.987, 1.005]; 8/21 | 8.401 [8.366, 8.464]; 21/21 strict win |
| v3000_3d-branched-spiral | roundtrip | 85.406 | 12.510 | 6.807 [6.774, 6.844]; 21/21 | 113.354 [113.204, 113.656]; 21/21 strict win |
| v3000_3d-offset | parse | 6.990 | 6.893 | 1.011 [1.005, 1.014]; 17/21 | 118.551 [117.244, 118.809]; 21/21 strict win |
| v3000_3d-offset | write | 2.490 | 2.464 | 1.007 [0.996, 1.021]; 14/21 | 8.790 [8.753, 8.829]; 21/21 strict win |
| v3000_3d-offset | roundtrip | 9.441 | 9.350 | 1.011 [1.005, 1.018]; 17/21 | 89.663 [89.486, 90.039]; 21/21 strict win |
| v3000_3d-branched-offset | parse | 82.500 | 10.131 | 8.154 [8.119, 8.191]; 21/21 | 138.244 [137.445, 139.149]; 21/21 strict win |
| v3000_3d-branched-offset | write | 2.505 | 2.510 | 0.999 [0.996, 1.011]; 10/21 | 8.670 [8.580, 8.728]; 21/21 strict win |
| v3000_3d-branched-offset | roundtrip | 85.156 | 12.605 | 6.771 [6.723, 6.791]; 21/21 | 112.649 [112.345, 113.007]; 21/21 strict win |

## Decision and every negative result

```json
{
  "local_adoption_gate_passed": false,
  "primary": {
    "sdf-0/write": {
      "median": 1.002,
      "ci95": [
        0.993,
        1.016
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "sdf-5000/write": {
      "median": 0.985,
      "ci95": [
        0.943,
        1.022
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/parse": {
      "median": 8.268,
      "ci95": [
        8.194,
        8.341
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-branched-offset/parse": {
      "median": 8.154,
      "ci95": [
        8.119,
        8.191
      ],
      "winning_blocks": 21,
      "blocks": 21
    }
  },
  "source_ratios": {
    "sdf-0/parse": {
      "median": 1.003,
      "ci95": [
        0.995,
        1.009
      ],
      "winning_blocks": 13,
      "blocks": 21
    },
    "sdf-0/write": {
      "median": 1.002,
      "ci95": [
        0.993,
        1.016
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "sdf-0/roundtrip": {
      "median": 0.998,
      "ci95": [
        0.995,
        1.003
      ],
      "winning_blocks": 10,
      "blocks": 21
    },
    "sdf-5000/parse": {
      "median": 1.011,
      "ci95": [
        0.995,
        1.018
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "sdf-5000/write": {
      "median": 0.985,
      "ci95": [
        0.943,
        1.022
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "sdf-5000/roundtrip": {
      "median": 0.997,
      "ci95": [
        0.986,
        1.007
      ],
      "winning_blocks": 9,
      "blocks": 21
    },
    "v3000-0/parse": {
      "median": 0.996,
      "ci95": [
        0.981,
        1.003
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "v3000-0/write": {
      "median": 1.002,
      "ci95": [
        0.987,
        1.014
      ],
      "winning_blocks": 12,
      "blocks": 21
    },
    "v3000-0/roundtrip": {
      "median": 1.005,
      "ci95": [
        0.997,
        1.024
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "v3000-5000/parse": {
      "median": 1.003,
      "ci95": [
        0.998,
        1.008
      ],
      "winning_blocks": 13,
      "blocks": 21
    },
    "v3000-5000/write": {
      "median": 1.034,
      "ci95": [
        1.027,
        1.042
      ],
      "winning_blocks": 20,
      "blocks": 21
    },
    "v3000-5000/roundtrip": {
      "median": 1.008,
      "ci95": [
        1.002,
        1.02
      ],
      "winning_blocks": 15,
      "blocks": 21
    },
    "v3000-chains/parse": {
      "median": 1.013,
      "ci95": [
        0.999,
        1.02
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "v3000-chains/write": {
      "median": 1.001,
      "ci95": [
        0.989,
        1.023
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "v3000-chains/roundtrip": {
      "median": 1.002,
      "ci95": [
        0.993,
        1.007
      ],
      "winning_blocks": 12,
      "blocks": 21
    },
    "v3000-continued/parse": {
      "median": 1.012,
      "ci95": [
        1.0,
        1.023
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "v3000-continued/write": {
      "median": 1.008,
      "ci95": [
        0.982,
        1.032
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "v3000-continued/roundtrip": {
      "median": 1.012,
      "ci95": [
        0.999,
        1.016
      ],
      "winning_blocks": 13,
      "blocks": 21
    },
    "v3000_3d-spiral/parse": {
      "median": 1.007,
      "ci95": [
        1.006,
        1.009
      ],
      "winning_blocks": 17,
      "blocks": 21
    },
    "v3000_3d-spiral/write": {
      "median": 1.006,
      "ci95": [
        0.998,
        1.015
      ],
      "winning_blocks": 13,
      "blocks": 21
    },
    "v3000_3d-spiral/roundtrip": {
      "median": 1.005,
      "ci95": [
        1.0,
        1.012
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/parse": {
      "median": 8.268,
      "ci95": [
        8.194,
        8.341
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/write": {
      "median": 0.992,
      "ci95": [
        0.987,
        1.005
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/roundtrip": {
      "median": 6.807,
      "ci95": [
        6.774,
        6.844
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-offset/parse": {
      "median": 1.011,
      "ci95": [
        1.005,
        1.014
      ],
      "winning_blocks": 17,
      "blocks": 21
    },
    "v3000_3d-offset/write": {
      "median": 1.007,
      "ci95": [
        0.996,
        1.021
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "v3000_3d-offset/roundtrip": {
      "median": 1.011,
      "ci95": [
        1.005,
        1.018
      ],
      "winning_blocks": 17,
      "blocks": 21
    },
    "v3000_3d-branched-offset/parse": {
      "median": 8.154,
      "ci95": [
        8.119,
        8.191
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-branched-offset/write": {
      "median": 0.999,
      "ci95": [
        0.996,
        1.011
      ],
      "winning_blocks": 10,
      "blocks": 21
    },
    "v3000_3d-branched-offset/roundtrip": {
      "median": 6.771,
      "ci95": [
        6.723,
        6.791
      ],
      "winning_blocks": 21,
      "blocks": 21
    }
  },
  "negative_intervals": [],
  "source_floor_failures": [
    {
      "lane": "sdf-5000/write",
      "speedup": {
        "median": 0.985,
        "ci95": [
          0.943,
          1.022
        ],
        "winning_blocks": 8,
        "blocks": 21
      },
      "unrounded_median": 0.9845157731632681
    }
  ],
  "all_source_snapshots_and_rewrites_exact": true,
  "source_valid_epochs": 10,
  "universal_superiority": false
}
```

## Evidence and reproduction

[Complete packet](2026-10-17-openbabel-degree-coordinate-artifacts/README.md) contains every raw matrix, exact input, complete common source output streams, retained Open Babel graph/rewrite outcomes, source/binary seals, compiler logs, disassembly, frozen-function proofs, full test logs and recomputed statistical audit. Original absolute paths in recipes require adaptation. Rebuild the named baseline and candidate in independent roots/targets using the original command in `baseline-reuse.json` and `build.py`; run `native.py`, `profile.py`, `seal.py`, then `formal.py`. `assemble.py` independently checks all source/ranked comparator ratios, the raw median floors, complete frozen functions and outcome hashes. Full local evidence: `/private/tmp/chematic-degree-filter-coordinate-evidence-20261017`. A later integration revision requires correctness validation and is not the timed candidate.
