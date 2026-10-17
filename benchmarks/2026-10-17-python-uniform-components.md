# Python uniform-component preparation — 2026-10-17

Decision: **FAIL: keep an unpromoted draft**. This source experiment measures first-use Labute ASA, the complete 42-lane retention frontier, and output-qualified RDKit/COSMolKit comparisons. It does not establish toolkit-wide or published-package superiority.

NCI first-use ASA improves 1.052× [1.031,1.081] against the base (17/20 faster), while ChEMBL is 1.004× [0.994,1.013] (12/20), so the second primary is inconclusive. ChEMBL warm TPSA misses the median retention floor: 0.989 [0.966,1.010]; its interval includes 1, so this is a floor failure without a statistically established slowdown. No source interval is wholly below 1. First-use ASA remains slower than both comparators: RDKit/candidate ratios 0.819 and 0.229, COSMolKit/candidate 0.399 and 0.106. Keep this prototype unpromoted and leave main unchanged.

## Change

Measured base `ad0283a8101dfc6f7d5ea6b6bdea94a35f8404f0`; measured candidate `3222303a4a5c4a1956ee82cb0e35529273239a12`. Baseline own-CPU samples attribute significant initial ASA work through the actual Python getter to aromatic descriptor preparation. The profiling process also parses/frees batches; only getter ancestry supports this finding. The first attempt to attach the second diagnostic sampler failed during process startup; its log remains retained, and a separate retry after startup succeeded. No diagnostic worker supplies formal timings.

For a molecule already accepted by the original RDKit-model safety gate, the new structural predicate proves that each cyclic component is fully explicit aromatic or saturated without an incident multiple bond. Aromatic edges must be cyclic with flagged endpoints; non-aromatic cyclic edges must have unflagged endpoints; cyclic multiple bonds and exocyclic multiple bonds touching a non-aromatic ring atom keep the original path. Neutral hypercoordinate group-16 states and the existing macrocycle donor rule retain their checks. The shared cyclic flags/components are prepared as before, but temporary component fact tables can be skipped when the original complete-preperceived result is provably Identity. The full fallback is unchanged. No ring algorithm, public API, dependency, unsafe code, or eager descriptor parsing is added.

## Correctness and provenance

The complete frozen base shortcut is verified byte-exact after renaming its declaration. It matches cold/repeated decisions for all 14,072 pinned graphs and their attempted Kekulé forms, all 4,096 five-atom/flag/H combinations, and directed fused/mixed/exocyclic/macrocycle/unusual cases. The first new-test compilation used a nonexistent convenience Kekulé API; the retained failed log supplied no formal timing, and the test was corrected to the existing core API before committing/building the candidate. Production source did not change during that correction.

Independent wheel processes retain exact complete descriptor, Morgan and MACCS streams: 112,576 descriptor outcomes (float bits), 56,288 Morgan outputs and 28,144 MACCS outcomes over 14,072 inputs, zero parse/evaluation errors. Six existing nonfinite descriptor outcomes remain represented by their exact bits; they are not counted as valid external chemistry agreement. The ordered four-property dump and the separate fresh-per-operation timing workers check different preparation contexts.

Rust release validation: 3,835 passed, 0 failed, 4 existing ignored across 137 suites (perception, SMILES, chemistry, fingerprints and MOL). Release all-target Clippy, formatting, refactor ratchet and benchmark-runner tests pass. Full Python suites: `====================== 3473 passed, 32 skipped in 18.12s =======================` / `======================= 3473 passed, 32 skipped in 9.01s =======================`. Thirty skips reflect APFS EILSEQ for non-UTF-8 filenames; the other two reflect optional elixcee/OGB dependencies unavailable in this environment.

Both wheels use separate initially empty targets and separate source roots, Rust opt3/full LTO/one CGU/one job, incremental=false, strip=none, native-inchi. Every baseline source file matches its Git revision. Both verbose logs compile the production perception library, and the actual shortcut/callback instruction hashes differ. Wheel/extension and source hashes are sealed, as are RDKit, COSMolKit and NumPy package-local native files. System libraries and package-external runtime data are outside that seal. All local builds, checks and diagnostics finish before timing; no old wheel or timing samples are reused.

## Fixed protocol

CPython 3.13.6, RDKit 2026.09.1, COSMolKit 0.5.0, NumPy 2.5.3 on one shared macOS arm64 host. Two unfiltered 1,000-record slices (offsets 0 and 5,000) of the pinned browser 10k corpus. Twenty rotating four-arm fresh-process blocks per slice, sequential, one complete source-valid epoch each. Primary: both first-use Labute ASA paired CI lower bounds >1. Adoption also requires exact source outcomes and all 42 source median ratios >=0.99. Every negative interval is disclosed; no adaptive timing reruns, sample pooling or substitution.

Scalar operations (TPSA, ASA, ring count) use **20 copies of each full slice, parsed into 20,000 separate objects**. Fingerprints use the original single 1,000-record batch. This fixed choice was made before formal timing to reduce noise in historically roughly 40-microsecond warm scalar batches. A state-recording factory spy independently confirms 61,001 object creations, 40,001 objects evaluated once and 20,000 evaluated twice for one scalar worker: the initial parse probe/disjoint warm object, fresh first-use objects, once-warmed repeated-use objects, and genuine parse-inclusive objects. All preflight elapsed fields are discarded. These are repeated input records, not additional independent chemical cases or statistical blocks; comparisons with historical 1,000-object absolute timings are not claimed.

First-use parsing occurs outside the operation timer and may compare different eager/lazy preparation. Repeated-use warms the operation once per object. Parse-inclusive performs actual parsing and evaluation inside the same timer. Native result objects are retained; normalization/output comparisons run outside timing. Whole-process peak RSS includes imports, preparation and all operations; scalar batch sizes are larger than historical matrices. No isolated-operation memory claim. External host load is uncontrolled and recorded.

Ratios are reference/candidate elapsed times. Intervals use the unchanged deterministic paired bootstrap (10,000 resamples, three-decimal reporting). A strict external win requires equivalent outputs, CI lower>1 and all 20 blocks faster; mismatches remain unranked. Scalar mismatch counts include the 20 input copies. Source equivalence is exact; external scalar tolerance is 1e-9, fingerprints/rings require exact values. No Linux, Windows, WASM, independent holdout or published-artifact inference.

## Complete matrix

Displayed milliseconds normalize marginal batch medians to 1,000 calls by dividing scalar timings by 20. They are averages within a 20,000-object batch, not separate 1,000-object experiments. Ratios/intervals use the original raw paired batches.

| Offset | Operation / mode | Base ms per 1k | Candidate ms per 1k | Base/candidate [95% CI]; wins | RDKit/candidate [95% CI]; wins | COSMolKit/candidate [95% CI]; wins |
|---|---|---:|---:|---|---|---|
| 0 | tpsa/first_use | 1.331 | 1.268 | 1.078 [1.049, 1.110]; 16/20 | 1.905 [1.813, 2.042]; 20/20 strict win | 0.827 [0.746, 0.885]; 3/20 |
| 0 | tpsa/repeated_use | 0.099 | 0.094 | 1.114 [0.973, 1.197]; 13/20 | 4.091 [3.725, 4.395]; 20/20 strict win | 1.145 [1.050, 1.240]; 15/20 |
| 0 | tpsa/parse_inclusive | 4.334 | 4.217 | 1.019 [1.011, 1.041]; 16/20 | 13.629 [13.563, 13.890]; 20/20 strict win | 12.665 [12.561, 12.691]; 20/20 strict win |
| 0 | labute_asa/first_use | 1.411 | 1.339 | 1.052 [1.031, 1.081]; 17/20 | 0.819 [0.792, 0.848]; 1/20 | 0.399 [0.389, 0.425]; 0/20 |
| 0 | labute_asa/repeated_use | 0.107 | 0.102 | 1.061 [0.972, 1.112]; 13/20 | 1.008 [0.954, 1.049]; 11/20 | 0.948 [0.883, 1.017]; 7/20 |
| 0 | labute_asa/parse_inclusive | 4.588 | 4.487 | 1.019 [1.002, 1.034]; 15/20 | 12.670 [12.421, 12.813]; 20/20 strict win | 11.969 [11.789, 12.157]; 20/20 strict win |
| 0 | num_rings/first_use | 1.633 | 1.611 | 1.005 [0.991, 1.024]; 11/20 | 0.038 [0.036, 0.040]; 0/20 | unranked: output mismatch |
| 0 | num_rings/repeated_use | 0.098 | 0.093 | 1.078 [1.003, 1.163]; 15/20 | 0.643 [0.570, 0.682]; 0/20 | unranked: output mismatch |
| 0 | num_rings/parse_inclusive | 4.723 | 4.702 | 1.008 [0.993, 1.014]; 12/20 | 11.821 [11.623, 11.908]; 20/20 strict win | unranked: output mismatch |
| 0 | morgan2_2048/first_use | 5.305 | 5.136 | 1.023 [1.011, 1.038]; 18/20 | 3.142 [3.094, 3.172]; 20/20 strict win | 2.043 [2.027, 2.076]; 20/20 strict win |
| 0 | morgan2_2048/repeated_use | 3.048 | 3.003 | 1.022 [0.998, 1.039]; 13/20 | 5.287 [5.141, 5.415]; 20/20 strict win | 3.314 [3.203, 3.630]; 20/20 strict win |
| 0 | morgan2_2048/parse_inclusive | 8.481 | 8.347 | 1.011 [1.002, 1.034]; 15/20 | 8.844 [8.702, 9.495]; 20/20 strict win | 7.682 [7.575, 8.260]; 20/20 strict win |
| 0 | morgan2_chiral/first_use | 14.153 | 13.995 | 1.007 [0.993, 1.026]; 11/20 | 1.175 [1.155, 1.278]; 20/20 strict win | 0.776 [0.758, 0.841]; 3/20 |
| 0 | morgan2_chiral/repeated_use | 3.188 | 3.169 | 0.997 [0.986, 1.006]; 8/20 | 5.039 [4.984, 5.110]; 20/20 strict win | 3.252 [3.118, 3.313]; 20/20 strict win |
| 0 | morgan2_chiral/parse_inclusive | 16.186 | 16.076 | 1.005 [0.988, 1.018]; 12/20 | 4.520 [4.468, 4.556]; 20/20 strict win | 3.965 [3.932, 4.168]; 20/20 strict win |
| 0 | morgan3_2048/first_use | 6.116 | 6.001 | 1.018 [1.008, 1.033]; 16/20 | 3.410 [3.375, 3.500]; 20/20 strict win | 2.276 [2.243, 2.391]; 20/20 strict win |
| 0 | morgan3_2048/repeated_use | 4.214 | 4.230 | 1.008 [0.978, 1.040]; 11/20 | 4.856 [4.760, 4.955]; 20/20 strict win | 3.134 [3.046, 3.267]; 20/20 strict win |
| 0 | morgan3_2048/parse_inclusive | 9.681 | 9.544 | 1.009 [0.988, 1.026]; 13/20 | 8.105 [7.983, 8.153]; 20/20 strict win | 7.054 [6.898, 7.978]; 20/20 strict win |
| 0 | maccs/first_use | 231.627 | 231.715 | 0.999 [0.996, 1.003]; 9/20 | 1.807 [1.785, 1.830]; 20/20 strict win | 2.077 [2.059, 2.151]; 20/20 strict win |
| 0 | maccs/repeated_use | 222.551 | 222.618 | 1.000 [0.998, 1.015]; 11/20 | 1.876 [1.846, 1.895]; 20/20 strict win | 2.160 [2.138, 2.230]; 20/20 strict win |
| 0 | maccs/parse_inclusive | 234.819 | 234.234 | 1.001 [0.997, 1.011]; 12/20 | 2.027 [2.013, 2.052]; 20/20 strict win | 2.282 [2.261, 2.300]; 20/20 strict win |
| 5000 | tpsa/first_use | 3.429 | 3.420 | 1.004 [0.991, 1.008]; 12/20 | 0.515 [0.508, 0.520]; 0/20 | 0.175 [0.170, 0.189]; 0/20 |
| 5000 | tpsa/repeated_use | 0.100 | 0.102 | 0.989 [0.966, 1.010]; 8/20 | 3.562 [3.498, 3.844]; 20/20 strict win | 0.982 [0.949, 1.083]; 8/20 |
| 5000 | tpsa/parse_inclusive | 4.897 | 4.913 | 1.000 [0.991, 1.029]; 10/20 | 5.305 [5.235, 5.402]; 20/20 strict win | 5.449 [5.360, 5.504]; 20/20 strict win |
| 5000 | labute_asa/first_use | 3.575 | 3.583 | 1.004 [0.994, 1.013]; 12/20 | 0.229 [0.220, 0.235]; 0/20 | 0.106 [0.103, 0.108]; 0/20 |
| 5000 | labute_asa/repeated_use | 0.099 | 0.100 | 0.997 [0.937, 1.036]; 10/20 | 1.005 [0.945, 1.012]; 11/20 | 0.854 [0.825, 0.910]; 3/20 |
| 5000 | labute_asa/parse_inclusive | 5.096 | 5.099 | 1.007 [1.002, 1.017]; 15/20 | 4.940 [4.891, 5.004]; 20/20 strict win | 5.245 [5.157, 5.581]; 20/20 strict win |
| 5000 | num_rings/first_use | 0.534 | 0.545 | 1.000 [0.975, 1.017]; 10/20 | 0.109 [0.108, 0.112]; 0/20 | 0.174 [0.168, 0.178]; 0/20 |
| 5000 | num_rings/repeated_use | 0.070 | 0.068 | 1.023 [0.991, 1.049]; 13/20 | 0.838 [0.758, 0.860]; 2/20 | 0.837 [0.792, 0.849]; 3/20 |
| 5000 | num_rings/parse_inclusive | 2.039 | 2.039 | 1.005 [0.990, 1.033]; 11/20 | 11.963 [11.864, 12.389]; 20/20 strict win | 12.870 [12.787, 14.707]; 20/20 strict win |
| 5000 | morgan2_2048/first_use | 5.672 | 5.679 | 1.001 [0.985, 1.019]; 10/20 | 1.558 [1.526, 1.582]; 19/20 | 0.998 [0.975, 1.035]; 8/20 |
| 5000 | morgan2_2048/repeated_use | 1.794 | 1.728 | 1.026 [1.005, 1.044]; 15/20 | 4.869 [4.769, 4.926]; 20/20 strict win | 3.003 [2.930, 3.117]; 20/20 strict win |
| 5000 | morgan2_2048/parse_inclusive | 7.297 | 7.218 | 1.002 [0.997, 1.017]; 11/20 | 4.606 [4.572, 4.660]; 20/20 strict win | 4.443 [4.398, 4.563]; 20/20 strict win |
| 5000 | morgan2_chiral/first_use | 5.567 | 5.484 | 1.019 [0.993, 1.039]; 14/20 | 1.604 [1.580, 1.661]; 20/20 strict win | 1.048 [1.029, 1.077]; 17/20 |
| 5000 | morgan2_chiral/repeated_use | 1.768 | 1.774 | 1.006 [0.985, 1.017]; 12/20 | 4.841 [4.788, 4.923]; 20/20 strict win | 2.935 [2.918, 2.985]; 20/20 strict win |
| 5000 | morgan2_chiral/parse_inclusive | 7.220 | 7.288 | 0.997 [0.987, 1.012]; 9/20 | 4.618 [4.578, 4.870]; 20/20 strict win | 4.434 [4.374, 4.501]; 20/20 strict win |
| 5000 | morgan3_2048/first_use | 6.074 | 6.026 | 1.005 [0.994, 1.031]; 11/20 | 1.780 [1.766, 1.848]; 20/20 strict win | 1.207 [1.179, 1.213]; 19/20 |
| 5000 | morgan3_2048/repeated_use | 2.331 | 2.275 | 1.006 [0.993, 1.044]; 11/20 | 4.645 [4.607, 4.682]; 20/20 strict win | 2.990 [2.945, 3.156]; 20/20 strict win |
| 5000 | morgan3_2048/parse_inclusive | 7.891 | 7.762 | 1.021 [0.993, 1.048]; 14/20 | 4.575 [4.542, 4.671]; 20/20 strict win | 4.324 [4.275, 4.431]; 20/20 strict win |
| 5000 | maccs/first_use | 138.785 | 138.334 | 1.008 [0.996, 1.016]; 13/20 | 1.578 [1.564, 1.596]; 20/20 strict win | 1.688 [1.678, 1.716]; 20/20 strict win |
| 5000 | maccs/repeated_use | 131.154 | 130.760 | 1.002 [0.996, 1.013]; 11/20 | 1.667 [1.650, 1.710]; 20/20 strict win | 1.775 [1.770, 1.787]; 20/20 strict win |
| 5000 | maccs/parse_inclusive | 140.662 | 140.165 | 1.004 [0.999, 1.011]; 13/20 | 1.739 [1.728, 1.783]; 20/20 strict win | 1.866 [1.858, 1.911]; 20/20 strict win |

## Decision and every negative interval/floor failure

```json
{
  "local_adoption_gate_passed": false,
  "primary": {
    "0": {
      "median": 1.052,
      "ci95": [
        1.031,
        1.081
      ],
      "winning_blocks": 17,
      "blocks": 20
    },
    "5000": {
      "median": 1.004,
      "ci95": [
        0.994,
        1.013
      ],
      "winning_blocks": 12,
      "blocks": 20
    }
  },
  "source_floor_failures": [
    {
      "lane": "5000/tpsa/repeated_use",
      "speedup": {
        "median": 0.989,
        "ci95": [
          0.966,
          1.01
        ],
        "winning_blocks": 8,
        "blocks": 20
      }
    }
  ],
  "negative_intervals": [],
  "all_14072_source_outcomes_exact": true,
  "source_valid_epochs_each_slice": 1,
  "universal_superiority": false
}
```

## Evidence and reproduction

[Complete frozen packet](2026-10-17-python-uniform-components-artifacts/README.md): plan, both raw matrices, every retained public output, original validation/profile errors, compiler logs, production disassembly, common instrumentation patch, manifests and independent recomputation. Full local evidence: `/private/tmp/chematic-python-asa-evidence-20261017`. Build the named revisions with `build.py` in separate targets; run `native.py`, `verify.py`, `check-instrumentation.py`, `seal.py`, then `formal.py`. Adapt the retained original absolute paths. `assemble.py` independently recomputes every source and ranked external ratio before writing this report. A later integration head, if any, requires fresh correctness checks and is not the timed candidate.
