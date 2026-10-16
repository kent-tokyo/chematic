# Relevant-cycle expansion — 2026-10-17

Local performance adoption gate: **FAIL**. RDKit/COSMolKit-wide superiority is not established.

This candidate remains unpromoted. Correctness checks do not override the failed performance gate.

## Change and correctness boundary

Measured baseline `7d6b092f448913a4ef87f1a3d6353789dad69cb5`; measured candidate `69c6d39a1353ef450c1a8ca93d27858d0962be09`; source version 1.1.6. The candidate compares two shortest paths at equal distances from their shared root. A shared vertex has one root distance, so this changes the intersection check from a nested membership scan to a single aligned scan without changing the predicate.

`relevant_cycles_with_limits` separates the output-ring limit from the per-endpoint shortest-path limit. Existing relevant-cycle APIs continue to use the same limit for both. The RDKit ring descriptor bounds output expansion at its existing three-times-legacy rejection threshold, retaining the original one-million path limit and legacy fallback. Above that output threshold the final descriptor was already the legacy count. No ring rule, atom order, sorting, graph eligibility, sanitization or cache rule changes. The broader relevant-cycle algorithm and its legacy compatibility limitations remain.

The test oracle freezes the complete original expansion and intersection functions from the baseline, verified after restoring renamed identifiers. All 1,024 simple five-vertex graphs compare exact ring sequence, every output cap through count+1, equal-limit errors, and every shortest-path pair. Cages, disconnected graphs, long multiword macrocycles, and a path/output-cap separation regression are also checked. Existing fixed-10k ring-count, edit-invalidation and arbitrary-graph differential tests remain active.

All 14,072 exposed inputs preserve exact source descriptor float bits/integers, Morgan bytes, and MACCS bytes: 112,576 descriptor outcomes, 56,288 Morgan outputs, and 28,144 MACCS outcomes. Six pre-existing nonfinite descriptor outcomes are retained and excluded from accuracy claims. This is source preservation, not a claim that every outcome matches RDKit or COSMolKit.

## Frozen protocol

Two unfiltered 1,000-row slices of the pinned 10k corpus, at offsets 0 and 5,000. Each has one fixed epoch of 20 rotating four-arm fresh-process blocks: baseline wheel, candidate wheel, RDKit and COSMolKit. Seven operations and three modes give 42 source lanes. Native return values stay alive during timing; normalization is outside. First-use molecules are freshly parsed outside timing; repeated use warms each molecule once; parse-inclusive executes actual parsing and the operation together. RDKit parsing eagerly prepares ring information, while chematic can defer work; first-use timings alone do not represent equal end-to-end work. TPSA is RDKit-default N/O.

Fresh source wheels use Rust 1.97.0, opt3/full LTO/one codegen unit/one job, native-inchi, strip=none. Comparators are RDKit 2026.09.1 and COSMolKit 0.5.0 on CPython 3.13.6, one shared macOS arm64 host. Complete source/input hashes, wheel/extension hashes and package-local comparator native hashes are checked before/after each profile. These hashes exclude system libraries and package-external runtime data. Load averages and whole-process peak RSS are retained; external host activity is uncontrolled.

Ratios are reference time / candidate time. Intervals use the unchanged deterministic paired-bootstrap procedure with 10,000 resamples. The primary is offset-0 first-use ring count, with 95% interval lower bound >1. All 42 source-relative medians must also be >=0.99 and all 14,072 source outcomes exact. Every negative interval and failed floor is disclosed. Strict competitor wins additionally require equivalent outputs and 20/20 faster blocks; unequal outputs are unranked. No adaptive reruns, pooled epochs or sample substitutions.

The initial pre-timing plan covered output limiting alone. Code inspection then identified the quadratic path intersection; the candidate description was amended before the candidate wheel and before any formal timing. The primary, corpora, block count and acceptance gates did not change. Both plans are retained. An older source wheel was used only for preliminary per-row diagnosis.

## Primary decision

Offset 0, `num_rings/first_use`, source ratio: **0.997 [0.992, 1.006]; 9/20**.

## Complete matrix: offset 0

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.990 [0.981, 1.000]; 6/20 | 1.481 [1.470, 1.522]; 20/20 **strict win** | 0.530 [0.527, 0.538]; 0/20 |
| tpsa/repeated_use | 1.002 [0.972, 1.019]; 10/20 | 9.764 [9.588, 10.115]; 20/20 **strict win** | 1.718 [1.675, 1.789]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.005 [0.998, 1.009]; 13/20 | 13.116 [13.024, 13.237]; 20/20 **strict win** | 12.230 [12.118, 12.267]; 20/20 **strict win** |
| labute_asa/first_use | 0.985 [0.975, 0.992]; 3/20 | 0.726 [0.703, 0.746]; 0/20 | 0.347 [0.342, 0.350]; 0/20 |
| labute_asa/repeated_use | 1.002 [0.986, 1.007]; 11/20 | 2.730 [2.691, 2.758]; 20/20 **strict win** | 1.479 [1.436, 1.571]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.003 [0.998, 1.008]; 12/20 | 12.323 [12.225, 12.410]; 20/20 **strict win** | 11.709 [11.659, 11.767]; 20/20 **strict win** |
| num_rings/first_use | 0.997 [0.992, 1.006]; 9/20 | 0.038 [0.037, 0.039]; 0/20 | unranked: unequal outputs |
| num_rings/repeated_use | 1.010 [0.994, 1.028]; 12/20 | 1.693 [1.680, 1.755]; 20/20 **strict win** | unranked: unequal outputs |
| num_rings/parse_inclusive | 1.003 [0.998, 1.015]; 12/20 | 11.554 [11.467, 11.653]; 20/20 **strict win** | unranked: unequal outputs |
| morgan2_2048/first_use | 1.000 [0.991, 1.009]; 10/20 | 3.117 [3.078, 3.174]; 20/20 **strict win** | 2.020 [2.013, 2.040]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.999 [0.993, 1.008]; 10/20 | 5.186 [5.142, 5.213]; 20/20 **strict win** | 3.189 [3.152, 3.233]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.005 [0.998, 1.009]; 13/20 | 8.602 [8.526, 8.677]; 20/20 **strict win** | 7.562 [7.491, 7.605]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.002 [1.000, 1.008]; 14/20 | 1.273 [1.256, 1.292]; 20/20 **strict win** | 0.827 [0.823, 0.832]; 0/20 |
| morgan2_chiral/repeated_use | 0.998 [0.992, 1.008]; 10/20 | 4.992 [4.969, 5.052]; 20/20 **strict win** | 3.098 [3.075, 3.110]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.004 [1.000, 1.009]; 14/20 | 4.467 [4.428, 4.491]; 20/20 **strict win** | 3.921 [3.905, 3.943]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.998 [0.990, 1.002]; 7/20 | 3.366 [3.329, 3.391]; 20/20 **strict win** | 2.242 [2.217, 2.253]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.998 [0.992, 1.002]; 9/20 | 4.810 [4.795, 4.875]; 20/20 **strict win** | 3.081 [3.072, 3.108]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.000 [0.998, 1.006]; 10/20 | 7.936 [7.880, 7.978]; 20/20 **strict win** | 6.949 [6.897, 6.973]; 20/20 **strict win** |
| maccs/first_use | 0.999 [0.995, 1.004]; 9/20 | 1.796 [1.791, 1.800]; 20/20 **strict win** | 2.061 [2.057, 2.070]; 20/20 **strict win** |
| maccs/repeated_use | 1.000 [0.992, 1.004]; 10/20 | 1.850 [1.842, 1.856]; 20/20 **strict win** | 2.128 [2.123, 2.135]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.002 [0.999, 1.007]; 13/20 | 2.005 [1.997, 2.015]; 20/20 **strict win** | 2.261 [2.252, 2.269]; 20/20 **strict win** |
## Complete matrix: offset 5000

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.000 [0.997, 1.004]; 10/20 | 0.444 [0.439, 0.453]; 0/20 | 0.120 [0.118, 0.122]; 0/20 |
| tpsa/repeated_use | 0.995 [0.983, 1.005]; 9/20 | 9.812 [9.701, 10.013]; 20/20 **strict win** | 1.347 [1.307, 1.356]; 20/20 **strict win** |
| tpsa/parse_inclusive | 0.995 [0.989, 1.005]; 7/20 | 5.184 [5.170, 5.235]; 20/20 **strict win** | 5.347 [5.326, 5.381]; 20/20 **strict win** |
| labute_asa/first_use | 0.996 [0.990, 1.001]; 6/20 | 0.205 [0.199, 0.209]; 0/20 | 0.086 [0.084, 0.088]; 0/20 |
| labute_asa/repeated_use | 0.998 [0.994, 1.022]; 8/20 | 2.661 [2.616, 2.702]; 20/20 **strict win** | 1.210 [1.179, 1.221]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.993 [0.986, 1.001]; 6/20 | 4.863 [4.845, 4.894]; 20/20 **strict win** | 5.189 [5.150, 5.224]; 20/20 **strict win** |
| num_rings/first_use | 0.993 [0.985, 1.000]; 6/20 | 0.104 [0.099, 0.107]; 0/20 | 0.087 [0.082, 0.092]; 0/20 |
| num_rings/repeated_use | 1.011 [0.982, 1.027]; 11/20 | 1.685 [1.667, 1.786]; 20/20 **strict win** | 1.128 [1.091, 1.148]; 19/20 |
| num_rings/parse_inclusive | 1.004 [0.996, 1.013]; 13/20 | 11.720 [11.582, 11.788]; 20/20 **strict win** | 12.795 [12.722, 12.844]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.000 [0.995, 1.007]; 10/20 | 1.555 [1.544, 1.569]; 20/20 **strict win** | 0.998 [0.988, 1.005]; 7/20 |
| morgan2_2048/repeated_use | 1.003 [0.998, 1.006]; 14/20 | 4.671 [4.653, 4.741]; 20/20 **strict win** | 2.866 [2.841, 2.879]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.001 [0.993, 1.009]; 11/20 | 4.609 [4.577, 4.653]; 20/20 **strict win** | 4.444 [4.435, 4.453]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.000 [0.995, 1.007]; 10/20 | 1.609 [1.594, 1.620]; 20/20 **strict win** | 1.027 [1.022, 1.030]; 20/20 **strict win** |
| morgan2_chiral/repeated_use | 0.998 [0.995, 1.005]; 8/20 | 4.748 [4.692, 4.786]; 20/20 **strict win** | 2.901 [2.870, 2.914]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.000 [0.993, 1.008]; 9/20 | 4.617 [4.589, 4.659]; 20/20 **strict win** | 4.454 [4.434, 4.465]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.000 [0.996, 1.004]; 10/20 | 1.791 [1.782, 1.803]; 20/20 **strict win** | 1.182 [1.173, 1.196]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.004 [1.002, 1.009]; 16/20 | 4.513 [4.497, 4.561]; 20/20 **strict win** | 2.877 [2.866, 2.937]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.001 [0.996, 1.009]; 10/20 | 4.514 [4.478, 4.530]; 20/20 **strict win** | 4.281 [4.274, 4.304]; 20/20 **strict win** |
| maccs/first_use | 1.000 [0.996, 1.006]; 10/20 | 1.571 [1.569, 1.579]; 20/20 **strict win** | 1.680 [1.676, 1.686]; 20/20 **strict win** |
| maccs/repeated_use | 1.000 [0.994, 1.003]; 9/20 | 1.661 [1.652, 1.664]; 20/20 **strict win** | 1.772 [1.765, 1.783]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.996 [0.993, 1.003]; 7/20 | 1.726 [1.721, 1.734]; 20/20 **strict win** | 1.851 [1.848, 1.856]; 20/20 **strict win** |

## Every negative source interval

```json
{
  "0": {
    "labute_asa/first_use": {
      "median": 0.985,
      "ci95": [
        0.975,
        0.992
      ],
      "winning_blocks": 3,
      "blocks": 20
    }
  },
  "5000": {}
}
```

## Every retention-floor failure

```json
{
  "0": {
    "labute_asa/first_use": {
      "median": 0.985,
      "ci95": [
        0.975,
        0.992
      ],
      "winning_blocks": 3,
      "blocks": 20
    }
  },
  "5000": {}
}
```

## Absolute batch times and process memory

Median milliseconds per complete 1,000-input batch.

| Offset / lane | Baseline | Candidate | RDKit | COSMolKit |
|---|---:|---:|---:|---:|
| 0 / tpsa/first_use | 1.315 | 1.326 | 1.987 | 0.704 |
| 0 / tpsa/repeated_use | 0.036 | 0.036 | 0.351 | 0.062 |
| 0 / tpsa/parse_inclusive | 4.295 | 4.269 | 56.040 | 52.204 |
| 0 / labute_asa/first_use | 1.387 | 1.404 | 1.018 | 0.489 |
| 0 / labute_asa/repeated_use | 0.035 | 0.035 | 0.097 | 0.052 |
| 0 / labute_asa/parse_inclusive | 4.484 | 4.453 | 54.946 | 52.133 |
| 0 / num_rings/first_use | 1.700 | 1.701 | 0.065 | 0.069 |
| 0 / num_rings/repeated_use | 0.032 | 0.032 | 0.055 | 0.037 |
| 0 / num_rings/parse_inclusive | 4.688 | 4.660 | 53.939 | 51.743 |
| 0 / morgan2_2048/first_use | 4.969 | 4.974 | 15.590 | 10.122 |
| 0 / morgan2_2048/repeated_use | 2.930 | 2.922 | 15.178 | 9.334 |
| 0 / morgan2_2048/parse_inclusive | 8.223 | 8.180 | 70.341 | 61.748 |
| 0 / morgan2_chiral/first_use | 12.330 | 12.309 | 15.658 | 10.217 |
| 0 / morgan2_chiral/repeated_use | 3.054 | 3.053 | 15.339 | 9.469 |
| 0 / morgan2_chiral/parse_inclusive | 15.817 | 15.741 | 70.436 | 61.846 |
| 0 / morgan3_2048/first_use | 5.909 | 5.896 | 19.884 | 13.175 |
| 0 / morgan3_2048/repeated_use | 4.019 | 4.024 | 19.408 | 12.443 |
| 0 / morgan3_2048/parse_inclusive | 9.389 | 9.355 | 74.347 | 65.087 |
| 0 / maccs/first_use | 226.295 | 226.496 | 406.910 | 467.350 |
| 0 / maccs/repeated_use | 219.017 | 219.608 | 405.438 | 467.230 |
| 0 / maccs/parse_inclusive | 232.021 | 231.626 | 464.043 | 523.570 |
| 5000 / tpsa/first_use | 3.581 | 3.577 | 1.589 | 0.428 |
| 5000 / tpsa/repeated_use | 0.035 | 0.036 | 0.351 | 0.048 |
| 5000 / tpsa/parse_inclusive | 4.868 | 4.886 | 25.430 | 26.170 |
| 5000 / labute_asa/first_use | 3.386 | 3.395 | 0.699 | 0.292 |
| 5000 / labute_asa/repeated_use | 0.035 | 0.035 | 0.094 | 0.042 |
| 5000 / labute_asa/parse_inclusive | 5.019 | 5.035 | 24.518 | 26.092 |
| 5000 / num_rings/first_use | 0.559 | 0.564 | 0.059 | 0.049 |
| 5000 / num_rings/repeated_use | 0.032 | 0.032 | 0.053 | 0.036 |
| 5000 / num_rings/parse_inclusive | 2.042 | 2.027 | 23.717 | 25.934 |
| 5000 / morgan2_2048/first_use | 5.426 | 5.419 | 8.447 | 5.416 |
| 5000 / morgan2_2048/repeated_use | 1.735 | 1.731 | 8.112 | 4.954 |
| 5000 / morgan2_2048/parse_inclusive | 7.062 | 7.056 | 32.542 | 31.301 |
| 5000 / morgan2_chiral/first_use | 5.278 | 5.282 | 8.509 | 5.437 |
| 5000 / morgan2_chiral/repeated_use | 1.737 | 1.740 | 8.269 | 5.052 |
| 5000 / morgan2_chiral/parse_inclusive | 7.035 | 7.046 | 32.593 | 31.281 |
| 5000 / morgan3_2048/first_use | 5.838 | 5.836 | 10.463 | 6.892 |
| 5000 / morgan3_2048/repeated_use | 2.263 | 2.254 | 10.199 | 6.495 |
| 5000 / morgan3_2048/parse_inclusive | 7.689 | 7.692 | 34.700 | 32.870 |
| 5000 / maccs/first_use | 136.779 | 136.743 | 214.954 | 229.746 |
| 5000 / maccs/repeated_use | 129.033 | 129.024 | 214.313 | 229.124 |
| 5000 / maccs/parse_inclusive | 138.458 | 138.861 | 239.997 | 257.031 |

RSS includes imports, preparation, warmup and all operations; it is not isolated operation memory.

| Offset | Baseline MiB | Candidate MiB | RDKit MiB | COSMolKit MiB |
|---|---:|---:|---:|---:|
| 0 | 46.9 | 46.8 | 106.5 | 106.3 |
| 5000 | 40.3 | 40.3 | 88.4 | 83.8 |

## Validation and retained evidence

Final native perception/chemistry tests: 1,543 passed, 0 failed, 7 ignored in 36 suites. Both affected crates passed release all-target Clippy. Format and refactor budget passed. Python: 3,314 passed, 32 skipped on each source wheel; 30 non-UTF8 filename cases are unsupported on this filesystem, plus two existing skips. Native tests use opt3/LTO=false/16 codegen units; these binaries supply no timing samples.

Pre-timing setup failures are retained: the initial sparse build-input snapshot lacked data files and was completed before a successful wheel build; native/Python validation needed restoration of tracked nested fixtures. Both wheels first hit 30 sandbox denials while creating non-UTF8 filename fixtures; rerunning with appropriate permissions passed all supported cases and explicitly skipped those 30 filesystem-unsupported cases. Initial denied logs are retained. The initial focused test pass covered output limiting; final native validation was rerun after the linear intersection and expanded differential tests. No failed or superseded check supplied observations.

[2026-10-17-python-ring-expansion-artifacts](2026-10-17-python-ring-expansion-artifacts/README.md) retains the frozen plans, complete raw timing matrices, every post-timing observer output and mismatch row, exact source dumps, source/wheel/native manifests, reconstructable candidate patch, build/validation/launch drivers and logs. Local source snapshots and the exact wheel binaries remain at `/private/tmp/chematic-ring-expansion-evidence-20261017`. Source reconstruction does not guarantee bit-reproducible binaries.

Additional hosts, sealed holdouts, other operations, 3D/reaction workloads and published Python/WASM packages remain unmeasured. Source intervals in untouched APIs do not by themselves establish a causal regression. A finite source comparison cannot prove universal or package-wide superiority.

