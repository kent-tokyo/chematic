# Open Babel comparison: decimal-pair coordinate formatting

Decision: **FAIL: retain as an unpromoted draft**. The two diverse SDF write lanes are the predeclared primary. All source outputs must match and all 30 source median ratios must be at least 0.99. This is one macOS source experiment, not universal superiority or published-package evidence.

## Change and correctness

Measured baseline `34f141266a81264f42990f98acbe25c91d9074af` (the independent ring-partner repair #986 over original main `f20304f4`); measured production candidate `774ddc180d7389b95ed30ab7fbae6bcab65f2d64`. Both arms contain the identical correctness repair, so their production-source difference is the coordinate formatter. The shared fixed-four-decimal coordinate formatter converts integer and fractional digits two at a time with a 200-byte ASCII pair table. Common field padding is included in the same stack buffer and appended once. The rounding arithmetic, negative sign including negative zero, and standard-format fallback for ties, large and nonfinite values are unchanged. The helper moves into its own module to keep the existing writer within its refactor budget. No public API, dependency or unsafe additions. A separate baseline own-CPU profile observes prepared SDF write mode; its elapsed value is excluded from all formal results and preparation ancestry must not be attributed to the coordinate formatter.

The complete frozen old function is exact after renaming its declaration. Both it and Rust standard formatting serve as independent byte oracles for every fractional pair, 100,000 deterministic random f64 bit patterns and directed rounding/sign/padding boundaries. The pre-existing oracle for 400,000 generated coordinates plus directed cases remains. Native MOL library/integration validation: 1386 passed, 0 failed, 0 existing ignored across 67 suites. Release all-target Clippy, formatting, unchanged refactor ratchet and ten benchmark-runner tests pass. All 4,384 native source parse snapshots and rewritten records match byte for byte.

## Fixed protocol and provenance

Ten fixed profiles: two unfiltered 1,000-row diverse SDF slices; corresponding V3000 slices; 64 linear chains; continued versions; two synthetic 3D chain profiles and two fluorine-branched versions. Synthetic chains span 16..1024 carbon atoms. Their coordinates are arbitrary deterministic finite nonplanar coordinates, not physical conformers, embedding quality or representative holdout molecules. Both SDF write source CI lower bounds must exceed 1; every one of the 30 unrounded paired source median ratios must be >=0.99. Negative intervals, semantic rejections and floor failures remain visible. One valid epoch per profile, no adaptive timing reruns or pooling.

Twenty-one rotating three-arm fresh-process blocks per profile; three native batch repetitions, parse/prepared write/roundtrip. Input loading, process startup and plugin discovery are excluded. V3000 framing is outside the timer, SDF framing inside. Diagnostic Rust readers retain full batches. Prepared writing receives one untimed warm write; SDF reuses a Rust String while V3000 and Open Babel allocate. Whole-process high-water RSS includes preparation, retained inputs and warmup; it is not operation-only or streaming memory. External host load is uncontrolled and logged per block. Bootstrap intervals describe these paired measurements, not every host or molecule.

Separate source roots and initially empty targets, Rust opt3/full LTO/one codegen unit/one job/incremental=false/strip=none. Baseline files are verified against the named Git base. Both verbose build logs compile the production MOL library; actual `push_fixed4` machine-code bodies differ. Before any formal samples, a mandatory native test reproduced an existing main regression in both baseline and candidate: an explicit-H ring CH counted itself as its own partner. The separately committed #986 repair excludes self-partners and retains distinct partners. Baseline MOL2 tests and the complete candidate MOL suite pass after applying it identically to both arms. Original build/seal records and both failures remain in `unmeasured-main-failure/`. Inputs/hashes, primary modes, block counts and all adoption thresholds are unchanged. Two unused candidate files were removed from the initial baseline copy before source sealing; baseline production never referenced them, no source code changed and no formal samples existed. Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27 native C++ probes use -O3; RDKit 2025.09.3 is the independent output observer. Seals cover tracked sources, baseline files, exact corpora, executables and package-local Open Babel files before/after each profile. System and package-external runtime data are outside that seal.

Ratios are reference/candidate elapsed times, deterministic paired 10,000-resample bootstrap, seeds 0 (source) and 1 (Open Babel), reported to three decimals. Source retention is assessed using unrounded medians. Open Babel parse ranking requires equivalent accepted canonical graphs. Write/roundtrip ranking additionally requires declared atom/bond fields, coordinates to 1e-4, names/properties and explicit unknown stereo. Rejected or mismatched arms remain unranked; a strict external win requires a CI lower bound above 1 and 21/21 faster blocks. Earlier source optimization ratios are not multiplied into these measurements.

## Complete matrix

Milliseconds per full batch, native elapsed divided by three repetitions. Paired ratios need not equal ratios of the displayed marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 12.564 | 12.535 | 1.000 [0.999, 1.007]; 9/21 | unranked: semantic prerequisite fails |
| sdf-0 | write | 3.845 | 3.760 | 1.024 [1.018, 1.035]; 18/21 | 6.292 [6.266, 6.354]; 21/21 strict win |
| sdf-0 | roundtrip | 16.667 | 16.520 | 1.009 [1.004, 1.016]; 17/21 | 7.165 [7.149, 7.218]; 21/21 strict win |
| sdf-5000 | parse | 4.772 | 4.762 | 1.001 [0.983, 1.004]; 11/21 | 8.802 [8.737, 8.894]; 21/21 strict win |
| sdf-5000 | write | 1.621 | 1.573 | 1.022 [1.012, 1.034]; 17/21 | 8.869 [8.794, 8.908]; 21/21 strict win |
| sdf-5000 | roundtrip | 6.477 | 6.470 | 1.007 [0.997, 1.014]; 14/21 | 8.759 [8.717, 8.789]; 21/21 strict win |
| v3000-0 | parse | 13.310 | 13.407 | 0.996 [0.992, 1.000]; 8/21 | unranked: semantic prerequisite fails |
| v3000-0 | write | 2.819 | 2.794 | 1.009 [1.001, 1.014]; 16/21 | unranked: semantic prerequisite fails |
| v3000-0 | roundtrip | 16.351 | 16.374 | 1.001 [0.992, 1.007]; 11/21 | unranked: semantic prerequisite fails |
| v3000-5000 | parse | 5.091 | 5.047 | 1.012 [1.003, 1.014]; 16/21 | unranked: semantic prerequisite fails |
| v3000-5000 | write | 0.710 | 0.700 | 1.012 [1.012, 1.026]; 20/21 | unranked: semantic prerequisite fails |
| v3000-5000 | roundtrip | 5.794 | 5.781 | 1.002 [0.998, 1.006]; 13/21 | unranked: semantic prerequisite fails |
| v3000-chains | parse | 5.725 | 5.687 | 1.006 [1.002, 1.013]; 16/21 | 142.956 [142.050, 143.885]; 21/21 strict win |
| v3000-chains | write | 1.243 | 1.232 | 1.013 [0.989, 1.021]; 12/21 | 15.464 [15.300, 15.521]; 21/21 strict win |
| v3000-chains | roundtrip | 6.964 | 6.916 | 1.008 [1.003, 1.013]; 17/21 | 121.223 [120.814, 122.104]; 21/21 strict win |
| v3000-continued | parse | 7.953 | 7.912 | 1.005 [1.001, 1.009]; 16/21 | 104.191 [103.713, 104.808]; 21/21 strict win |
| v3000-continued | write | 1.240 | 1.234 | 1.005 [0.996, 1.008]; 13/21 | unranked: semantic prerequisite fails |
| v3000-continued | roundtrip | 9.254 | 9.140 | 1.006 [1.003, 1.017]; 16/21 | unranked: semantic prerequisite fails |
| v3000_3d-spiral | parse | 6.866 | 6.807 | 1.004 [0.999, 1.010]; 14/21 | 120.045 [119.498, 120.805]; 21/21 strict win |
| v3000_3d-spiral | write | 2.519 | 2.458 | 1.022 [1.012, 1.035]; 17/21 | 8.683 [8.538, 8.828]; 21/21 strict win |
| v3000_3d-spiral | roundtrip | 9.339 | 9.226 | 1.014 [1.009, 1.020]; 18/21 | 91.143 [90.211, 91.964]; 21/21 strict win |
| v3000_3d-branched-spiral | parse | 83.015 | 84.690 | 0.984 [0.981, 0.995]; 5/21 | unranked: semantic prerequisite fails |
| v3000_3d-branched-spiral | write | 2.558 | 2.482 | 1.009 [1.004, 1.036]; 15/21 | 8.631 [8.428, 8.664]; 21/21 strict win |
| v3000_3d-branched-spiral | roundtrip | 85.942 | 86.491 | 0.994 [0.983, 1.000]; 6/21 | 16.425 [16.282, 16.594]; 21/21 strict win |
| v3000_3d-offset | parse | 7.059 | 6.929 | 1.022 [1.018, 1.030]; 19/21 | 118.571 [117.745, 120.006]; 21/21 strict win |
| v3000_3d-offset | write | 2.507 | 2.437 | 1.030 [1.019, 1.040]; 21/21 | 9.075 [8.890, 9.111]; 21/21 strict win |
| v3000_3d-offset | roundtrip | 9.548 | 9.365 | 1.019 [1.012, 1.026]; 20/21 | 90.009 [89.337, 90.520]; 21/21 strict win |
| v3000_3d-branched-offset | parse | 82.780 | 83.502 | 0.991 [0.989, 0.997]; 5/21 | 16.790 [16.644, 16.826]; 21/21 strict win |
| v3000_3d-branched-offset | write | 2.489 | 2.442 | 1.025 [1.013, 1.039]; 18/21 | 8.922 [8.859, 8.997]; 21/21 strict win |
| v3000_3d-branched-offset | roundtrip | 85.313 | 85.445 | 0.998 [0.992, 0.999]; 5/21 | 16.613 [16.527, 16.699]; 21/21 strict win |

## Decision and all negative results

```json
{
  "local_adoption_gate_passed": false,
  "primary": {
    "sdf-0": {
      "median": 1.024,
      "ci95": [
        1.018,
        1.035
      ],
      "winning_blocks": 18,
      "blocks": 21
    },
    "sdf-5000": {
      "median": 1.022,
      "ci95": [
        1.012,
        1.034
      ],
      "winning_blocks": 17,
      "blocks": 21
    }
  },
  "source_ratios": {
    "sdf-0/parse": {
      "median": 1.0,
      "ci95": [
        0.999,
        1.007
      ],
      "winning_blocks": 9,
      "blocks": 21
    },
    "sdf-0/write": {
      "median": 1.024,
      "ci95": [
        1.018,
        1.035
      ],
      "winning_blocks": 18,
      "blocks": 21
    },
    "sdf-0/roundtrip": {
      "median": 1.009,
      "ci95": [
        1.004,
        1.016
      ],
      "winning_blocks": 17,
      "blocks": 21
    },
    "sdf-5000/parse": {
      "median": 1.001,
      "ci95": [
        0.983,
        1.004
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "sdf-5000/write": {
      "median": 1.022,
      "ci95": [
        1.012,
        1.034
      ],
      "winning_blocks": 17,
      "blocks": 21
    },
    "sdf-5000/roundtrip": {
      "median": 1.007,
      "ci95": [
        0.997,
        1.014
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "v3000-0/parse": {
      "median": 0.996,
      "ci95": [
        0.992,
        1.0
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "v3000-0/write": {
      "median": 1.009,
      "ci95": [
        1.001,
        1.014
      ],
      "winning_blocks": 16,
      "blocks": 21
    },
    "v3000-0/roundtrip": {
      "median": 1.001,
      "ci95": [
        0.992,
        1.007
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "v3000-5000/parse": {
      "median": 1.012,
      "ci95": [
        1.003,
        1.014
      ],
      "winning_blocks": 16,
      "blocks": 21
    },
    "v3000-5000/write": {
      "median": 1.012,
      "ci95": [
        1.012,
        1.026
      ],
      "winning_blocks": 20,
      "blocks": 21
    },
    "v3000-5000/roundtrip": {
      "median": 1.002,
      "ci95": [
        0.998,
        1.006
      ],
      "winning_blocks": 13,
      "blocks": 21
    },
    "v3000-chains/parse": {
      "median": 1.006,
      "ci95": [
        1.002,
        1.013
      ],
      "winning_blocks": 16,
      "blocks": 21
    },
    "v3000-chains/write": {
      "median": 1.013,
      "ci95": [
        0.989,
        1.021
      ],
      "winning_blocks": 12,
      "blocks": 21
    },
    "v3000-chains/roundtrip": {
      "median": 1.008,
      "ci95": [
        1.003,
        1.013
      ],
      "winning_blocks": 17,
      "blocks": 21
    },
    "v3000-continued/parse": {
      "median": 1.005,
      "ci95": [
        1.001,
        1.009
      ],
      "winning_blocks": 16,
      "blocks": 21
    },
    "v3000-continued/write": {
      "median": 1.005,
      "ci95": [
        0.996,
        1.008
      ],
      "winning_blocks": 13,
      "blocks": 21
    },
    "v3000-continued/roundtrip": {
      "median": 1.006,
      "ci95": [
        1.003,
        1.017
      ],
      "winning_blocks": 16,
      "blocks": 21
    },
    "v3000_3d-spiral/parse": {
      "median": 1.004,
      "ci95": [
        0.999,
        1.01
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "v3000_3d-spiral/write": {
      "median": 1.022,
      "ci95": [
        1.012,
        1.035
      ],
      "winning_blocks": 17,
      "blocks": 21
    },
    "v3000_3d-spiral/roundtrip": {
      "median": 1.014,
      "ci95": [
        1.009,
        1.02
      ],
      "winning_blocks": 18,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/parse": {
      "median": 0.984,
      "ci95": [
        0.981,
        0.995
      ],
      "winning_blocks": 5,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/write": {
      "median": 1.009,
      "ci95": [
        1.004,
        1.036
      ],
      "winning_blocks": 15,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/roundtrip": {
      "median": 0.994,
      "ci95": [
        0.983,
        1.0
      ],
      "winning_blocks": 6,
      "blocks": 21
    },
    "v3000_3d-offset/parse": {
      "median": 1.022,
      "ci95": [
        1.018,
        1.03
      ],
      "winning_blocks": 19,
      "blocks": 21
    },
    "v3000_3d-offset/write": {
      "median": 1.03,
      "ci95": [
        1.019,
        1.04
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-offset/roundtrip": {
      "median": 1.019,
      "ci95": [
        1.012,
        1.026
      ],
      "winning_blocks": 20,
      "blocks": 21
    },
    "v3000_3d-branched-offset/parse": {
      "median": 0.991,
      "ci95": [
        0.989,
        0.997
      ],
      "winning_blocks": 5,
      "blocks": 21
    },
    "v3000_3d-branched-offset/write": {
      "median": 1.025,
      "ci95": [
        1.013,
        1.039
      ],
      "winning_blocks": 18,
      "blocks": 21
    },
    "v3000_3d-branched-offset/roundtrip": {
      "median": 0.998,
      "ci95": [
        0.992,
        0.999
      ],
      "winning_blocks": 5,
      "blocks": 21
    }
  },
  "negative_intervals": [
    {
      "lane": "v3000_3d-branched-spiral/parse",
      "speedup": {
        "median": 0.984,
        "ci95": [
          0.981,
          0.995
        ],
        "winning_blocks": 5,
        "blocks": 21
      }
    },
    {
      "lane": "v3000_3d-branched-offset/parse",
      "speedup": {
        "median": 0.991,
        "ci95": [
          0.989,
          0.997
        ],
        "winning_blocks": 5,
        "blocks": 21
      }
    },
    {
      "lane": "v3000_3d-branched-offset/roundtrip",
      "speedup": {
        "median": 0.998,
        "ci95": [
          0.992,
          0.999
        ],
        "winning_blocks": 5,
        "blocks": 21
      }
    }
  ],
  "source_floor_failures": [
    {
      "lane": "v3000_3d-branched-spiral/parse",
      "speedup": {
        "median": 0.984,
        "ci95": [
          0.981,
          0.995
        ],
        "winning_blocks": 5,
        "blocks": 21
      },
      "unrounded_median": 0.9838152766910273
    }
  ],
  "all_source_snapshots_and_rewrites_exact": true,
  "source_valid_epochs": 10,
  "universal_superiority": false
}
```

## Evidence and reproduction

[Complete frozen packet](2026-10-17-openbabel-coordinate-pairs-artifacts/README.md), including all ten raw matrices, exact inputs, complete shared source outcome streams and retained Open Babel outcomes, source/binary/native seals, compiler logs, disassembly, test logs and independently recomputed statistics. Full local evidence: `/private/tmp/chematic-coordinate-pairs-evidence-20261017`. Rebuild the named revisions in independent roots/targets with `build.py`; run `native.py`, `seal.py` and `formal.py` with retained original absolute paths adapted. `assemble.py` independently verifies source ratios, eligible competitor ratios and the adoption decision. A later merge revision requires correctness validation and is not the timed candidate.
