# Scalar descriptor model reuse — 2026-10-16

Local adoption gate: **FAIL**. RDKit/COSMolKit-wide superiority remains unproved. All measured lanes are retained.

This candidate must remain unpromoted: the predeclared primary or retention gate did not pass. Correctness CI does not override this performance decision.

## Change and correctness

Base `2251a9f5d2e706087efaf55cb0f8ab1ff20b8016`, source version 1.1.5. Only scalar RDKit-model callbacks memoize `Option<Molecule>` after sanitization. Different scalar descriptors on the same original immutable graph can reuse that model. Generic callbacks keep their prior uncached model and aromatic-view behavior. The conversion algorithm, chemical gate, failure fallback and nested TLS bypass are unchanged. Clone and graph edits clear the derived model slot. No unsafe code, new dependency or altered descriptor formula.

This adds one `OnceLock<Arc<…>>` slot to every molecule and retains the corrected model (including any derived views) until that original molecule is dropped or edited. It trades memory for potential reuse. No universal or per-operation memory reduction is claimed.

The regression oracle freezes the complete original model callback function; `frozen-reference.json` verifies tokens after only renaming and removing whitespace. All 14,072 pinned input graphs preserve atom/bond data, adjacency, stereo neighbor order, bond directions, stereo groups/group bonds, tags and R-group labels, with both callback modes and repeated reads. Directed tests cover hydrogen/isotope handling, specified stereo, radicals/charges, dative cleanup, failed sanitization fallback, clone/edit invalidation, shared scalar model identity, ordinary graphs that need no model, generic callback isolation, nested calls, panic unwinding and concurrent cold readers.

The two source wheels also preserve all public descriptor float bits/integers, Morgan bytes and MACCS bytes across those 14,072 inputs: 112,576 descriptor outcomes, 56,288 Morgan outputs and 28,144 MACCS outcomes. Six pre-existing nonfinite descriptor outcomes are retained and excluded from accuracy claims. Development inputs are exposed, not a sealed holdout.

## Frozen protocol

Two unfiltered 1,000-row slices at offsets 0 and 5,000 of the pinned 10,000-input corpus. Four-arm rotating fresh processes: baseline, candidate, RDKit and COSMolKit. One 20-block standard matrix and one 20-block TPSA-then-ASA bundle matrix per slice. Standard operations: TPSA, Labute ASA, ring count, Morgan2 2048, chiral Morgan2, Morgan3 2048 and MACCS. Each has first-use, repeated-use and actual parse-inclusive modes: 42 standard plus six bundle source lanes. No adaptive timing rerun, sample pooling or substitution.

The bundle calls RDKit-default N/O-only TPSA followed by Labute ASA on the same molecule and retains both native return values inside timing. First use parses fresh objects outside timing; repeated use warms the entire operation/bundle once per object; parse inclusive parses and computes inside timing. All normalization and validation happen afterwards. The ordinary seven-operation runner is unchanged. Eager perception during comparator parsing differs from lazy first-use work. Imports/startup are excluded.

Rust 1.97.0; fresh CPython 3.13.6 wheels; opt3/full LTO/one codegen unit/one job, native-inchi, strip=none. RDKit 2026.09.1 and COSMolKit 0.5.0 native package files are rehashed before/after profiles. System/transitive libraries and package-external data are excluded, so this is not a full runtime closure. One shared macOS arm64 host, uncontrolled external activity; per-block load averages and whole-process RSS are retained. These are source results, not released-package or whole-toolkit wins.

Ratios are reference time / candidate time; brackets are deterministic paired-bootstrap 95% intervals (10,000 resamples, three-decimal rounding). Adoption requires exact public source outputs, CI lower >1 for both bundle first-use lanes, and all 48 source medians >=0.99. Every negative source interval is disclosed. A strict competitor win also requires matching outputs and 20/20 faster blocks. Different outputs remain unranked.

## Primary result

| Profile | Source ratio |
|---|---:|
| bundle-0 / tpsa_then_asa/first_use | 1.013 [1.010, 1.022]; 17/20 |
| bundle-5000 / tpsa_then_asa/first_use | 1.045 [1.040, 1.050]; 20/20 |

## Complete matrix: frontier-0

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.008 [0.996, 1.018]; 12/20 | 1.522 [1.489, 1.535]; 20/20 **strict win** | 0.550 [0.539, 0.580]; 0/20 |
| tpsa/repeated_use | 1.009 [1.001, 1.068]; 14/20 | 9.914 [9.619, 10.046]; 20/20 **strict win** | 1.736 [1.666, 2.015]; 20/20 **strict win** |
| tpsa/parse_inclusive | 0.978 [0.972, 0.985]; 3/20 | 12.755 [12.686, 12.852]; 20/20 **strict win** | 11.888 [11.838, 12.067]; 20/20 **strict win** |
| labute_asa/first_use | 1.006 [0.991, 1.034]; 12/20 | 0.750 [0.730, 0.762]; 0/20 | 0.367 [0.357, 0.383]; 0/20 |
| labute_asa/repeated_use | 0.996 [0.990, 1.042]; 8/20 | 2.752 [2.696, 2.790]; 20/20 **strict win** | 1.477 [1.372, 1.581]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.984 [0.970, 0.994]; 5/20 | 12.126 [12.018, 12.182]; 20/20 **strict win** | 11.426 [11.362, 11.492]; 20/20 **strict win** |
| num_rings/first_use | 1.003 [0.990, 1.015]; 11/20 | 0.038 [0.037, 0.039]; 0/20 | unranked: unequal output |
| num_rings/repeated_use | 1.002 [0.985, 1.011]; 10/20 | 1.731 [1.700, 1.772]; 20/20 **strict win** | unranked: unequal output |
| num_rings/parse_inclusive | 0.988 [0.973, 0.994]; 3/20 | 11.357 [11.284, 11.473]; 20/20 **strict win** | unranked: unequal output |
| morgan2_2048/first_use | 0.998 [0.992, 1.015]; 9/20 | 3.105 [3.077, 3.148]; 20/20 **strict win** | 2.032 [2.020, 2.056]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.997 [0.990, 1.012]; 10/20 | 5.200 [5.158, 5.279]; 20/20 **strict win** | 3.194 [3.163, 3.228]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.989 [0.985, 0.998]; 6/20 | 8.548 [8.472, 8.588]; 20/20 **strict win** | 7.471 [7.427, 7.517]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.997 [0.990, 1.004]; 8/20 | 1.267 [1.251, 1.286]; 20/20 **strict win** | 0.830 [0.828, 0.836]; 0/20 |
| morgan2_chiral/repeated_use | 0.997 [0.988, 1.012]; 9/20 | 4.998 [4.961, 5.095]; 20/20 **strict win** | 3.121 [3.081, 3.157]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.993 [0.986, 1.003]; 7/20 | 4.449 [4.404, 4.484]; 20/20 **strict win** | 3.904 [3.845, 3.938]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.997 [0.988, 1.013]; 9/20 | 3.360 [3.326, 3.403]; 20/20 **strict win** | 2.240 [2.228, 2.260]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.999 [0.990, 1.003]; 7/20 | 4.880 [4.850, 4.906]; 20/20 **strict win** | 3.119 [3.114, 3.153]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.987 [0.983, 1.003]; 6/20 | 7.874 [7.804, 8.016]; 20/20 **strict win** | 6.882 [6.823, 6.947]; 20/20 **strict win** |
| maccs/first_use | 1.004 [0.999, 1.010]; 13/20 | 1.807 [1.796, 1.812]; 20/20 **strict win** | 2.088 [2.062, 2.111]; 20/20 **strict win** |
| maccs/repeated_use | 1.004 [0.999, 1.009]; 13/20 | 1.864 [1.850, 1.874]; 20/20 **strict win** | 2.144 [2.135, 2.165]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.002 [0.995, 1.006]; 12/20 | 2.015 [2.007, 2.026]; 20/20 **strict win** | 2.274 [2.262, 2.291]; 20/20 **strict win** |

## Complete matrix: frontier-5000

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.999 [0.991, 1.007]; 9/20 | 0.454 [0.450, 0.457]; 0/20 | 0.125 [0.121, 0.131]; 0/20 |
| tpsa/repeated_use | 1.027 [0.989, 1.066]; 14/20 | 9.807 [9.510, 10.185]; 20/20 **strict win** | 1.347 [1.316, 1.392]; 20/20 **strict win** |
| tpsa/parse_inclusive | 0.990 [0.985, 0.999]; 5/20 | 5.218 [5.180, 5.265]; 20/20 **strict win** | 5.335 [5.291, 5.412]; 20/20 **strict win** |
| labute_asa/first_use | 1.005 [0.993, 1.009]; 11/20 | 0.211 [0.207, 0.218]; 0/20 | 0.093 [0.089, 0.096]; 0/20 |
| labute_asa/repeated_use | 1.017 [0.984, 1.039]; 13/20 | 2.747 [2.631, 2.802]; 20/20 **strict win** | 1.168 [1.138, 1.290]; 19/20 |
| labute_asa/parse_inclusive | 0.997 [0.987, 1.008]; 8/20 | 4.836 [4.807, 4.913]; 20/20 **strict win** | 5.229 [5.121, 5.254]; 20/20 **strict win** |
| num_rings/first_use | 1.021 [1.003, 1.029]; 15/20 | 0.107 [0.104, 0.111]; 0/20 | 0.090 [0.081, 0.103]; 0/20 |
| num_rings/repeated_use | 0.999 [0.990, 1.017]; 9/20 | 1.778 [1.721, 1.858]; 20/20 **strict win** | 1.113 [1.103, 1.132]; 19/20 |
| num_rings/parse_inclusive | 1.002 [0.995, 1.014]; 12/20 | 11.666 [11.517, 11.857]; 20/20 **strict win** | 12.704 [12.502, 12.793]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.000 [0.994, 1.009]; 10/20 | 1.566 [1.546, 1.589]; 20/20 **strict win** | 0.997 [0.994, 1.000]; 6/20 |
| morgan2_2048/repeated_use | 1.000 [0.996, 1.009]; 11/20 | 4.676 [4.599, 4.749]; 20/20 **strict win** | 2.839 [2.801, 2.895]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.998 [0.994, 1.008]; 9/20 | 4.600 [4.578, 4.638]; 20/20 **strict win** | 4.433 [4.409, 4.466]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.998 [0.995, 1.010]; 9/20 | 1.609 [1.593, 1.622]; 20/20 **strict win** | 1.025 [1.021, 1.059]; 19/20 |
| morgan2_chiral/repeated_use | 0.995 [0.993, 1.004]; 7/20 | 4.699 [4.674, 4.758]; 20/20 **strict win** | 2.919 [2.846, 2.934]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.003 [0.997, 1.012]; 12/20 | 4.614 [4.584, 4.643]; 20/20 **strict win** | 4.444 [4.390, 4.499]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.002 [0.990, 1.012]; 11/20 | 1.794 [1.761, 1.813]; 20/20 **strict win** | 1.176 [1.165, 1.211]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.999 [0.987, 1.002]; 8/20 | 4.542 [4.492, 4.586]; 20/20 **strict win** | 2.924 [2.889, 2.939]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.000 [0.992, 1.008]; 12/20 | 4.534 [4.497, 4.564]; 20/20 **strict win** | 4.285 [4.262, 4.322]; 20/20 **strict win** |
| maccs/first_use | 1.003 [0.999, 1.010]; 14/20 | 1.573 [1.566, 1.577]; 20/20 **strict win** | 1.680 [1.671, 1.695]; 20/20 **strict win** |
| maccs/repeated_use | 1.004 [0.998, 1.009]; 12/20 | 1.661 [1.653, 1.667]; 20/20 **strict win** | 1.780 [1.764, 1.789]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.003 [1.003, 1.009]; 15/20 | 1.734 [1.729, 1.743]; 20/20 **strict win** | 1.864 [1.849, 1.871]; 20/20 **strict win** |

## Complete matrix: bundle-0

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa_then_asa/first_use | 1.013 [1.010, 1.022]; 17/20 | 1.529 [1.520, 1.547]; 20/20 **strict win** | 0.536 [0.524, 0.543]; 0/20 |
| tpsa_then_asa/repeated_use | 1.006 [0.970, 1.020]; 11/20 | 4.915 [4.855, 4.990]; 20/20 **strict win** | 1.204 [1.193, 1.253]; 20/20 **strict win** |
| tpsa_then_asa/parse_inclusive | 1.016 [1.005, 1.019]; 16/20 | 11.706 [11.685, 11.784]; 20/20 **strict win** | 10.842 [10.802, 10.980]; 20/20 **strict win** |

## Complete matrix: bundle-5000

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa_then_asa/first_use | 1.045 [1.040, 1.050]; 20/20 | 0.583 [0.576, 0.590]; 0/20 | 0.162 [0.161, 0.164]; 0/20 |
| tpsa_then_asa/repeated_use | 0.994 [0.970, 1.014]; 9/20 | 4.817 [4.721, 4.903]; 20/20 **strict win** | 1.161 [1.117, 1.204]; 20/20 **strict win** |
| tpsa_then_asa/parse_inclusive | 1.029 [1.026, 1.040]; 20/20 | 4.882 [4.854, 4.892]; 20/20 **strict win** | 4.930 [4.906, 4.944]; 20/20 **strict win** |

## Every negative source interval

```json
{
  "frontier-0": {
    "tpsa/parse_inclusive": {
      "median": 0.978,
      "ci95": [
        0.972,
        0.985
      ],
      "winning_blocks": 3,
      "blocks": 20
    },
    "labute_asa/parse_inclusive": {
      "median": 0.984,
      "ci95": [
        0.97,
        0.994
      ],
      "winning_blocks": 5,
      "blocks": 20
    },
    "num_rings/parse_inclusive": {
      "median": 0.988,
      "ci95": [
        0.973,
        0.994
      ],
      "winning_blocks": 3,
      "blocks": 20
    },
    "morgan2_2048/parse_inclusive": {
      "median": 0.989,
      "ci95": [
        0.985,
        0.998
      ],
      "winning_blocks": 6,
      "blocks": 20
    }
  },
  "frontier-5000": {
    "tpsa/parse_inclusive": {
      "median": 0.99,
      "ci95": [
        0.985,
        0.999
      ],
      "winning_blocks": 5,
      "blocks": 20
    }
  },
  "bundle-0": {},
  "bundle-5000": {}
}
```

## Every source median below 0.99

```json
{
  "frontier-0": {
    "tpsa/parse_inclusive": {
      "median": 0.978,
      "ci95": [
        0.972,
        0.985
      ],
      "winning_blocks": 3,
      "blocks": 20
    },
    "labute_asa/parse_inclusive": {
      "median": 0.984,
      "ci95": [
        0.97,
        0.994
      ],
      "winning_blocks": 5,
      "blocks": 20
    },
    "num_rings/parse_inclusive": {
      "median": 0.988,
      "ci95": [
        0.973,
        0.994
      ],
      "winning_blocks": 3,
      "blocks": 20
    },
    "morgan2_2048/parse_inclusive": {
      "median": 0.989,
      "ci95": [
        0.985,
        0.998
      ],
      "winning_blocks": 6,
      "blocks": 20
    },
    "morgan3_2048/parse_inclusive": {
      "median": 0.987,
      "ci95": [
        0.983,
        1.003
      ],
      "winning_blocks": 6,
      "blocks": 20
    }
  },
  "frontier-5000": {},
  "bundle-0": {},
  "bundle-5000": {}
}
```

## Absolute times

Medians in milliseconds per complete 1,000-row batch.

| Profile / operation / mode | Baseline | Candidate | RDKit | COSMolKit |
|---|---:|---:|---:|---:|
| frontier-0 / tpsa/first_use | 1.337 | 1.324 | 2.019 | 0.725 |
| frontier-0 / tpsa/repeated_use | 0.037 | 0.036 | 0.355 | 0.067 |
| frontier-0 / tpsa/parse_inclusive | 4.326 | 4.408 | 56.274 | 52.868 |
| frontier-0 / labute_asa/first_use | 1.409 | 1.398 | 1.039 | 0.521 |
| frontier-0 / labute_asa/repeated_use | 0.035 | 0.036 | 0.099 | 0.053 |
| frontier-0 / labute_asa/parse_inclusive | 4.499 | 4.570 | 55.250 | 52.338 |
| frontier-0 / num_rings/first_use | 1.707 | 1.702 | 0.065 | 0.074 |
| frontier-0 / num_rings/repeated_use | 0.032 | 0.032 | 0.056 | 0.037 |
| frontier-0 / num_rings/parse_inclusive | 4.702 | 4.806 | 54.162 | 52.091 |
| frontier-0 / morgan2_2048/first_use | 5.002 | 5.008 | 15.588 | 10.202 |
| frontier-0 / morgan2_2048/repeated_use | 2.927 | 2.922 | 15.228 | 9.352 |
| frontier-0 / morgan2_2048/parse_inclusive | 8.214 | 8.264 | 70.889 | 62.112 |
| frontier-0 / morgan2_chiral/first_use | 12.397 | 12.369 | 15.814 | 10.378 |
| frontier-0 / morgan2_chiral/repeated_use | 3.061 | 3.062 | 15.504 | 9.567 |
| frontier-0 / morgan2_chiral/parse_inclusive | 15.826 | 15.941 | 71.156 | 62.454 |
| frontier-0 / morgan3_2048/first_use | 5.864 | 5.886 | 19.915 | 13.252 |
| frontier-0 / morgan3_2048/repeated_use | 3.996 | 4.010 | 19.540 | 12.479 |
| frontier-0 / morgan3_2048/parse_inclusive | 9.363 | 9.512 | 74.946 | 65.495 |
| frontier-0 / maccs/first_use | 227.539 | 226.986 | 408.564 | 470.926 |
| frontier-0 / maccs/repeated_use | 220.284 | 219.736 | 407.705 | 469.979 |
| frontier-0 / maccs/parse_inclusive | 231.838 | 231.948 | 467.010 | 526.488 |
| frontier-5000 / tpsa/first_use | 3.600 | 3.617 | 1.651 | 0.457 |
| frontier-5000 / tpsa/repeated_use | 0.037 | 0.037 | 0.358 | 0.050 |
| frontier-5000 / tpsa/parse_inclusive | 4.885 | 4.923 | 25.806 | 26.540 |
| frontier-5000 / labute_asa/first_use | 3.404 | 3.429 | 0.737 | 0.318 |
| frontier-5000 / labute_asa/repeated_use | 0.037 | 0.037 | 0.100 | 0.043 |
| frontier-5000 / labute_asa/parse_inclusive | 5.067 | 5.105 | 24.642 | 26.343 |
| frontier-5000 / num_rings/first_use | 0.576 | 0.568 | 0.061 | 0.051 |
| frontier-5000 / num_rings/repeated_use | 0.032 | 0.032 | 0.056 | 0.035 |
| frontier-5000 / num_rings/parse_inclusive | 2.045 | 2.068 | 23.968 | 25.996 |
| frontier-5000 / morgan2_2048/first_use | 5.487 | 5.507 | 8.581 | 5.480 |
| frontier-5000 / morgan2_2048/repeated_use | 1.754 | 1.762 | 8.231 | 5.022 |
| frontier-5000 / morgan2_2048/parse_inclusive | 7.102 | 7.087 | 32.736 | 31.498 |
| frontier-5000 / morgan2_chiral/first_use | 5.315 | 5.336 | 8.620 | 5.542 |
| frontier-5000 / morgan2_chiral/repeated_use | 1.761 | 1.762 | 8.300 | 5.119 |
| frontier-5000 / morgan2_chiral/parse_inclusive | 7.191 | 7.065 | 32.796 | 31.783 |
| frontier-5000 / morgan3_2048/first_use | 5.921 | 5.938 | 10.614 | 7.094 |
| frontier-5000 / morgan3_2048/repeated_use | 2.258 | 2.279 | 10.332 | 6.658 |
| frontier-5000 / morgan3_2048/parse_inclusive | 7.745 | 7.712 | 35.135 | 33.405 |
| frontier-5000 / maccs/first_use | 137.926 | 137.714 | 218.020 | 230.933 |
| frontier-5000 / maccs/repeated_use | 131.354 | 130.333 | 216.749 | 230.575 |
| frontier-5000 / maccs/parse_inclusive | 139.533 | 139.435 | 241.893 | 258.533 |
| bundle-0 / tpsa_then_asa/first_use | 1.908 | 1.883 | 2.888 | 1.008 |
| bundle-0 / tpsa_then_asa/repeated_use | 0.099 | 0.098 | 0.484 | 0.120 |
| bundle-0 / tpsa_then_asa/parse_inclusive | 4.937 | 4.868 | 56.895 | 52.911 |
| bundle-5000 / tpsa_then_asa/first_use | 4.255 | 4.071 | 2.388 | 0.662 |
| bundle-5000 / tpsa_then_asa/repeated_use | 0.098 | 0.100 | 0.485 | 0.118 |
| bundle-5000 / tpsa_then_asa/parse_inclusive | 5.568 | 5.407 | 26.340 | 26.546 |

## Whole-process memory

RSS includes imports, preparation, warmup and all worker operations; it does not isolate model retention or operation memory.

| Profile | Baseline MiB | Candidate MiB | RDKit MiB | COSMolKit MiB |
|---|---:|---:|---:|---:|
| frontier-0 | 46.6 | 46.8 | 106.6 | 106.2 |
| frontier-5000 | 40.0 | 40.4 | 88.4 | 83.8 |
| bundle-0 | 34.3 | 34.4 | 97.7 | 96.9 |
| bundle-5000 | 34.3 | 34.7 | 83.9 | 78.1 |

## Verification and retained evidence

Native core/smiles/chemistry tests: 1,939 passed, 0 failed, 3 ignored across 58 suites (`--lib --tests`, debug profile with debug info disabled). Core/smiles all-target Clippy, format and refactor-budget checks passed. Python: 2,302 passed, two skipped. Three bundle semantic claim-gate tests passed. These checks are distinct from the optimized timed wheels.

Pre-timing failures were corrected before any formal worker: a new test initially referenced the wrong path for DerivedScalar; disk exhaustion interrupted candidate compilation and native test linking; two incomplete, nonexecutable linker outputs were removed and relinked. The stopped build/test logs and targeted cleanup records are retained. No failed validation supplied timing observations. The baseline snapshot includes only named Git-tracked build inputs, excluding ignored Python caches.

[2026-10-16-python-scalar-model-reuse-artifacts](2026-10-16-python-scalar-model-reuse-artifacts/README.md) contains all raw matrices, fresh post-timing observer outputs and all COSMolKit mismatch rows, decisions, frozen plan, source/native/wheel manifests, candidate patch, drivers and logs. Full wheels, installed extensions, source snapshot and 14,072-input compressed output streams remain at `/private/tmp/chematic-scalar-model-reuse-evidence-20261016` with hashes. Source reconstruction is verified; bit-reproducible rebuilds are not claimed.

More hosts, sealed holdout data, other operations, 3D/reactions, CLI startup and published Python/WASM artifacts remain open. A negative interval on an untouched API does not alone establish a causal code regression, but still counts under the fixed retention rule.

