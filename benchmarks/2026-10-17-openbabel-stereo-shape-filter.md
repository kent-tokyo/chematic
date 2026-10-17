# Open Babel comparison: reject impossible local stereocenter shapes before branch traversal

Decision: **FAIL: retain as an unpromoted draft**. Exact source outputs are retained; primary and retention results below determine adoption. This is a source experiment on one macOS host, not universal superiority or a published-package result.

Both branched synthetic 3D parse primaries improve 7.170–7.304× against the base and 119.243–120.589× against Open Babel under the graph-equivalence prerequisite, with 21/21 faster blocks. Three source write lanes fail the 0.99 median floor: SDF slice 0 (0.988), SDF slice 5000 (0.963), and continued V3000 (0.980). These correspond to about 1.2%, 3.9%, and 2.0% longer elapsed time. The candidate remains a draft and is not adopted into main.

## Change and correctness

Measured base `de41dadd8f123a4b1368937c7dd86df823434588`; measured candidate `812bf0a953789299cdb1c7cf4a81133a0afe14ba`. In `stereo_centers`, CH₂ atoms can satisfy `groups == 4` but fail the later local output-shape condition. The unchanged condition now rejects those atoms before `substituent_signatures_distinct` builds rooted branch signatures and clones visited HashSets. The branch algorithm and valid-candidate order remain unchanged. No public API/dependency/unsafe additions. A baseline own-CPU sample confirms branch traversal/allocation is a major load on the branched 3D primary; that diagnostic worker's timing is excluded from every formal sample.

The complete frozen classifier is byte-exact to the base after only renaming its declaration. All 14,072 pinned graphs, all 6,144 five-atom/hydrogen-override graph combinations, and directed ring-nitrogen/explicit-H/charged/isotope/annotated boundaries match. Native validation: 3295 passed, 0 failed, 3 ignored across 110 suites, covering perception, MOL, chemistry and fingerprints. Release all-target Clippy, formatting, refactor ratchet and ten benchmark-runner tests pass. All 4,384 profile source snapshots and rewrites are byte-exact. The frozen classifier supplies direct public-result retention; this experiment does not add a Python/RDKit/COS performance claim.

## Fixed protocol and build provenance

Ten profiles: two 1,000-row diverse SDF slices, two corresponding V3000 slices, 64 linear chains, continued versions, two synthetic nonplanar linear-chain profiles, and two branched versions with F attached to carbon 2. Each synthetic set spans 16..1024 carbons. Coordinates are deterministic finite arbitrary coordinates, not physical conformers, embedding quality, or a fresh representative holdout. Both branched 3D parse lanes are the predeclared primary. Adoption requires exact source results, both primary paired 95% CI lower bounds above 1, and all 30 source medians at least 0.99. Every negative interval is disclosed. No sample filtering, epoch pooling or adaptive reruns.

Twenty-one rotating three-arm fresh-process blocks per profile, three native batch repetitions, parse/prepared write/roundtrip. Startup, file loading and plugin discovery are excluded; V3000 delimiter splitting is outside the timer, SDF framing inside. Parse retains the full diagnostic batch. Prepared writing has one untimed warm write. Rust SDF writing reuses a String; V3000 and Open Babel allocate strings. RSS is whole-process high-water usage including setup, batch retention and warmup, not isolated-operation or streaming memory. External host load is uncontrolled and recorded per block.

Rust release opt3/full LTO/one codegen unit/one job, incremental=false, strip=none. Separate initially empty targets and separate source roots for both arms. Every baseline source file is checked against the base Git revision. Both verbose logs explicitly compile the production perception library. Full LTO inlines `stereo_centers` into `perceive_tetrahedral_from_3d_with_mode`; that actual production owner has different machine-code bodies in the two arms, independently of example executable hashes. The initial preflight rejected a generated continuation-line prefix with one space instead of two. That failed setup supplied zero formal timing samples; the defective corpus, plan, seal and log remain in `setup-failure/`. Only that framing prefix was corrected, with unchanged graphs, rows, source, hypothesis and gates. No reused timing samples. Open Babel 3.2.1 native C++ probes are built with -O3; RDKit 2025.09.3 observes graph and rewrite contracts. Source, baseline, corpus, binary and package-local Open Babel files are sealed before/after profiles; system and package-external runtime data are outside that seal.

Ratios are reference/candidate times with the unchanged deterministic 10,000-resample paired bootstrap. Open Babel parse ranking requires equivalent canonical graphs; write/roundtrip ranking requires accepted equivalent rewrites with declared atom/bond fields, coordinates to 1e-4, names/properties and unknown stereo. Mismatches/rejections remain present and unranked. A strict comparator win also requires 21/21 faster blocks. The earlier PR #966 whole-molecule guard and PR #960 writer optimization are already in this base; these ratios must not be multiplied by their historical ratios.

## Complete matrix

Milliseconds per full batch (native elapsed divided by three repetitions); ratios are paired medians, not ratios of the displayed marginal medians.

| Profile | Mode | Base ms | Candidate ms | Base/candidate [95% CI]; wins | Open Babel/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 24.787 | 24.771 | 1.001 [1.000, 1.006]; 14/21 | unranked: semantic prerequisite fails |
| sdf-0 | write | 7.716 | 7.901 | 0.988 [0.976, 0.995]; 6/21 | 5.867 [5.770, 5.998]; 21/21 strict win |
| sdf-0 | roundtrip | 32.815 | 32.919 | 0.997 [0.992, 1.001]; 8/21 | 7.100 [6.981, 7.149]; 21/21 strict win |
| sdf-5000 | parse | 9.461 | 9.424 | 1.004 [0.989, 1.013]; 11/21 | 8.862 [8.750, 8.969]; 21/21 strict win |
| sdf-5000 | write | 3.169 | 3.290 | 0.963 [0.957, 0.973]; 3/21 | 8.228 [8.071, 8.282]; 21/21 strict win |
| sdf-5000 | roundtrip | 12.781 | 12.816 | 1.004 [0.988, 1.016]; 11/21 | 8.561 [8.431, 8.671]; 21/21 strict win |
| v3000-0 | parse | 26.501 | 26.266 | 1.004 [0.987, 1.012]; 12/21 | unranked: semantic prerequisite fails |
| v3000-0 | write | 5.749 | 5.772 | 1.005 [0.968, 1.015]; 11/21 | unranked: semantic prerequisite fails |
| v3000-0 | roundtrip | 32.303 | 32.278 | 1.001 [0.979, 1.017]; 13/21 | unranked: semantic prerequisite fails |
| v3000-5000 | parse | 10.332 | 10.349 | 1.003 [0.974, 1.017]; 12/21 | unranked: semantic prerequisite fails |
| v3000-5000 | write | 1.378 | 1.401 | 0.991 [0.973, 1.003]; 9/21 | unranked: semantic prerequisite fails |
| v3000-5000 | roundtrip | 11.730 | 11.753 | 1.000 [0.984, 1.012]; 10/21 | unranked: semantic prerequisite fails |
| v3000-chains | parse | 11.134 | 11.110 | 0.999 [0.994, 1.007]; 9/21 | 141.988 [141.453, 143.359]; 21/21 strict win |
| v3000-chains | write | 2.435 | 2.462 | 0.993 [0.975, 1.012]; 7/21 | 15.102 [14.827, 15.317]; 21/21 strict win |
| v3000-chains | roundtrip | 13.681 | 13.613 | 1.001 [0.997, 1.012]; 12/21 | 119.276 [118.481, 119.756]; 21/21 strict win |
| v3000-continued | parse | 15.574 | 15.613 | 0.997 [0.988, 1.011]; 8/21 | 102.533 [100.905, 108.409]; 21/21 strict win |
| v3000-continued | write | 2.458 | 2.562 | 0.980 [0.952, 1.005]; 7/21 | unranked: semantic prerequisite fails |
| v3000-continued | roundtrip | 18.087 | 18.265 | 0.993 [0.988, 1.000]; 7/21 | unranked: semantic prerequisite fails |
| v3000_3d-spiral | parse | 6.896 | 6.948 | 1.000 [0.991, 1.004]; 11/21 | 117.504 [115.761, 118.262]; 21/21 strict win |
| v3000_3d-spiral | write | 2.484 | 2.502 | 0.999 [0.973, 1.008]; 9/21 | 8.503 [8.377, 8.627]; 21/21 strict win |
| v3000_3d-spiral | roundtrip | 9.359 | 9.440 | 0.990 [0.985, 0.997]; 6/21 | 88.192 [87.948, 89.759]; 21/21 strict win |
| v3000_3d-branched-spiral | parse | 84.721 | 11.638 | 7.304 [7.189, 7.339]; 21/21 | 120.589 [119.084, 121.818]; 21/21 strict win |
| v3000_3d-branched-spiral | write | 2.509 | 2.490 | 1.007 [1.000, 1.016]; 15/21 | 8.578 [8.466, 8.653]; 21/21 strict win |
| v3000_3d-branched-spiral | roundtrip | 87.037 | 14.038 | 6.202 [6.144, 6.239]; 21/21 | 101.001 [100.635, 101.749]; 21/21 strict win |
| v3000_3d-offset | parse | 6.883 | 6.939 | 0.992 [0.989, 0.999]; 6/21 | 118.317 [116.473, 119.378]; 21/21 strict win |
| v3000_3d-offset | write | 2.454 | 2.460 | 0.992 [0.984, 1.007]; 8/21 | 8.862 [8.841, 8.947]; 21/21 strict win |
| v3000_3d-offset | roundtrip | 9.407 | 9.416 | 0.995 [0.992, 1.003]; 7/21 | 89.442 [88.630, 90.349]; 21/21 strict win |
| v3000_3d-branched-offset | parse | 84.906 | 11.777 | 7.170 [7.104, 7.227]; 21/21 | 119.243 [117.582, 120.452]; 21/21 strict win |
| v3000_3d-branched-offset | write | 2.504 | 2.498 | 1.003 [0.982, 1.012]; 11/21 | 8.881 [8.780, 8.963]; 21/21 strict win |
| v3000_3d-branched-offset | roundtrip | 88.140 | 14.271 | 6.136 [6.105, 6.222]; 21/21 | 100.763 [99.598, 101.566]; 21/21 strict win |

## Decision, all negative intervals and floor failures

```json
{
  "local_adoption_gate_passed": false,
  "primary": {
    "v3000_3d-branched-spiral": {
      "median": 7.304,
      "ci95": [
        7.189,
        7.339
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-branched-offset": {
      "median": 7.17,
      "ci95": [
        7.104,
        7.227
      ],
      "winning_blocks": 21,
      "blocks": 21
    }
  },
  "source_ratios": {
    "sdf-0/parse": {
      "median": 1.001,
      "ci95": [
        1.0,
        1.006
      ],
      "winning_blocks": 14,
      "blocks": 21
    },
    "sdf-0/write": {
      "median": 0.988,
      "ci95": [
        0.976,
        0.995
      ],
      "winning_blocks": 6,
      "blocks": 21
    },
    "sdf-0/roundtrip": {
      "median": 0.997,
      "ci95": [
        0.992,
        1.001
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "sdf-5000/parse": {
      "median": 1.004,
      "ci95": [
        0.989,
        1.013
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "sdf-5000/write": {
      "median": 0.963,
      "ci95": [
        0.957,
        0.973
      ],
      "winning_blocks": 3,
      "blocks": 21
    },
    "sdf-5000/roundtrip": {
      "median": 1.004,
      "ci95": [
        0.988,
        1.016
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "v3000-0/parse": {
      "median": 1.004,
      "ci95": [
        0.987,
        1.012
      ],
      "winning_blocks": 12,
      "blocks": 21
    },
    "v3000-0/write": {
      "median": 1.005,
      "ci95": [
        0.968,
        1.015
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "v3000-0/roundtrip": {
      "median": 1.001,
      "ci95": [
        0.979,
        1.017
      ],
      "winning_blocks": 13,
      "blocks": 21
    },
    "v3000-5000/parse": {
      "median": 1.003,
      "ci95": [
        0.974,
        1.017
      ],
      "winning_blocks": 12,
      "blocks": 21
    },
    "v3000-5000/write": {
      "median": 0.991,
      "ci95": [
        0.973,
        1.003
      ],
      "winning_blocks": 9,
      "blocks": 21
    },
    "v3000-5000/roundtrip": {
      "median": 1.0,
      "ci95": [
        0.984,
        1.012
      ],
      "winning_blocks": 10,
      "blocks": 21
    },
    "v3000-chains/parse": {
      "median": 0.999,
      "ci95": [
        0.994,
        1.007
      ],
      "winning_blocks": 9,
      "blocks": 21
    },
    "v3000-chains/write": {
      "median": 0.993,
      "ci95": [
        0.975,
        1.012
      ],
      "winning_blocks": 7,
      "blocks": 21
    },
    "v3000-chains/roundtrip": {
      "median": 1.001,
      "ci95": [
        0.997,
        1.012
      ],
      "winning_blocks": 12,
      "blocks": 21
    },
    "v3000-continued/parse": {
      "median": 0.997,
      "ci95": [
        0.988,
        1.011
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "v3000-continued/write": {
      "median": 0.98,
      "ci95": [
        0.952,
        1.005
      ],
      "winning_blocks": 7,
      "blocks": 21
    },
    "v3000-continued/roundtrip": {
      "median": 0.993,
      "ci95": [
        0.988,
        1.0
      ],
      "winning_blocks": 7,
      "blocks": 21
    },
    "v3000_3d-spiral/parse": {
      "median": 1.0,
      "ci95": [
        0.991,
        1.004
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "v3000_3d-spiral/write": {
      "median": 0.999,
      "ci95": [
        0.973,
        1.008
      ],
      "winning_blocks": 9,
      "blocks": 21
    },
    "v3000_3d-spiral/roundtrip": {
      "median": 0.99,
      "ci95": [
        0.985,
        0.997
      ],
      "winning_blocks": 6,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/parse": {
      "median": 7.304,
      "ci95": [
        7.189,
        7.339
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/write": {
      "median": 1.007,
      "ci95": [
        1.0,
        1.016
      ],
      "winning_blocks": 15,
      "blocks": 21
    },
    "v3000_3d-branched-spiral/roundtrip": {
      "median": 6.202,
      "ci95": [
        6.144,
        6.239
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-offset/parse": {
      "median": 0.992,
      "ci95": [
        0.989,
        0.999
      ],
      "winning_blocks": 6,
      "blocks": 21
    },
    "v3000_3d-offset/write": {
      "median": 0.992,
      "ci95": [
        0.984,
        1.007
      ],
      "winning_blocks": 8,
      "blocks": 21
    },
    "v3000_3d-offset/roundtrip": {
      "median": 0.995,
      "ci95": [
        0.992,
        1.003
      ],
      "winning_blocks": 7,
      "blocks": 21
    },
    "v3000_3d-branched-offset/parse": {
      "median": 7.17,
      "ci95": [
        7.104,
        7.227
      ],
      "winning_blocks": 21,
      "blocks": 21
    },
    "v3000_3d-branched-offset/write": {
      "median": 1.003,
      "ci95": [
        0.982,
        1.012
      ],
      "winning_blocks": 11,
      "blocks": 21
    },
    "v3000_3d-branched-offset/roundtrip": {
      "median": 6.136,
      "ci95": [
        6.105,
        6.222
      ],
      "winning_blocks": 21,
      "blocks": 21
    }
  },
  "negative_intervals": [
    {
      "lane": "sdf-0/write",
      "speedup": {
        "median": 0.988,
        "ci95": [
          0.976,
          0.995
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    },
    {
      "lane": "sdf-5000/write",
      "speedup": {
        "median": 0.963,
        "ci95": [
          0.957,
          0.973
        ],
        "winning_blocks": 3,
        "blocks": 21
      }
    },
    {
      "lane": "v3000_3d-spiral/roundtrip",
      "speedup": {
        "median": 0.99,
        "ci95": [
          0.985,
          0.997
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    },
    {
      "lane": "v3000_3d-offset/parse",
      "speedup": {
        "median": 0.992,
        "ci95": [
          0.989,
          0.999
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    }
  ],
  "source_floor_failures": [
    {
      "lane": "sdf-0/write",
      "speedup": {
        "median": 0.988,
        "ci95": [
          0.976,
          0.995
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    },
    {
      "lane": "sdf-5000/write",
      "speedup": {
        "median": 0.963,
        "ci95": [
          0.957,
          0.973
        ],
        "winning_blocks": 3,
        "blocks": 21
      }
    },
    {
      "lane": "v3000-continued/write",
      "speedup": {
        "median": 0.98,
        "ci95": [
          0.952,
          1.005
        ],
        "winning_blocks": 7,
        "blocks": 21
      }
    }
  ],
  "all_source_snapshots_and_rewrites_exact": true,
  "source_valid_epochs": 10,
  "universal_superiority": false
}
```

## Evidence and reproduction

[Complete frozen packet](2026-10-17-openbabel-stereo-shape-filter-artifacts/README.md), including all ten raw timing matrices, inputs, source and executable seals, compiler logs, disassembly, tests and audit. Full local evidence: `/private/tmp/chematic-stereo-shape-evidence-20261017`. Build the two named source revisions using `build.py` in independent roots/targets, run `native.py`, `seal.py`, then `formal.py`; adapt retained original absolute paths. `assemble.py` independently recomputes every source/ranked comparator ratio and records every lane before assembling this report. A later merge revision, if any, requires fresh correctness validation and is not the timed candidate.
