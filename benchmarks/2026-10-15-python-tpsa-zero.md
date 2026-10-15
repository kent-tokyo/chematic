# Python default TPSA zero-contribution guard — 2026-10-15

The production change is adopted in the attached source candidate. This is unmerged source-version 1.1.4 evidence, without current-main verification or hosted CI. Complete superiority over RDKit or COSMolKit is not established.

Decision: Adopt the default TPSA zero-contribution guard: NCI first-use median 1.028 [1.022,1.039], 18/20 faster, and parse-inclusive 1.026 [1.019,1.035], 17/20 faster, establish meaningful intended median improvements with exact output parity. Neither is a strict all-block win. ChEMBL first-use TPSA 0.996 [0.991,1.010] remains an uncertain small difference. No initial source-relative CI upper is below 1 across 42 lanes, so the predeclared adaptive rule selects no follow-ups. The rejected indexed-normalization prototype and its separate observations remain excluded. This does not establish complete superiority over either comparator.

## Change and correctness

The RDKit-default TPSA implementation checks the selected model view for nitrogen or oxygen before requesting its aromatic descriptor view. If neither element occurs, every term of the original positive-zero fold is zero, including an empty molecule, so it returns the same positive zero. The full original N/O computation, its floating order and model scalar cache are unchanged. Parsing and the S/P-inclusive TPSA function are unchanged. No public API, dependency or unsafe code is added. Rustfmt also formats the adjacent existing descriptor-view test; an independent rustfmt comparison proves no logic change. An initial manifest preflight caught that previously unlisted test formatting before any timing block; its failure is retained and the exact already-built/tested source is listed in the final manifest.

The immediate baseline includes all fourteen accepted source candidates; rejected prototypes are excluded. Its exact freshly built wheel is reused from the preceding separately archived rejected indexed-normalization experiment, with all 10,061 current source/input hashes verified unchanged. The baseline was built before either prototype edit; the new candidate wheel and every timing here are new. baseline-reuse.json records the wheel identity.

An independent frozen copy of the entire original default TPSA implementation checks exact floating bits on 14,072 pinned inputs, parsed and attempted Kekule forms (original clone on Kekule error), with cold and repeated calls: 56,288 comparisons. Empty inputs, wildcard, carbon-only aromatic/Kekule forms, S/P-only, metal, explicit-H and late-N/O cases also match. Independent processes produce byte-identical complete descriptor and Morgan dumps: 112,576 descriptor outcomes and 56,288 Morgan outputs, with zero parse/evaluation errors. Existing pinned corpora are reused.

Scoped release Rust tests: 957 passed, 0 failed, 1 existing ignored. Python: 1,389 passed, one skipped. Release all-target Clippy for chemistry, formatting and the unchanged refactor budget (487 production files, 8,063 functions, 18 file allowances, 88 function allowances) passed. Full logs are retained. Builds and CPU sampling completed before formal timing; only owned completed build targets were removed.

Pre-timing structural coverage finds 43/1,000 NCI inputs with no N/O and 0/1,000 ChEMBL inputs. All 43 original baseline default TPSA outputs are exact positive zero. The full original unfiltered slices remain the timed populations; no record is selected or removed. The preceding experiment and its initial/adaptive negatives remain in [the prior rejected experiment](2026-10-15-python-aromatic-indexed.md); observations are not pooled and this guard does not retroactively establish a benefit for that rejected prototype.

## Protocol

Pinned environment: CPython 3.13.6, RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, Rust 1.97.0, maturin 1.13.3, macOS arm64. Two fixed unfiltered 1,000-input slices, offsets 5,000 then 0, use 20 rotating fresh-process four-arm blocks, seven operations and all three preparation modes. This is one exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.

Ratios are reference elapsed / candidate elapsed; greater than 1 is faster. Brackets are paired median-ratio 95% bootstrap intervals. A strict win requires semantic agreement, CI lower > 1 and all 20 blocks faster. Output mismatches suppress ranking; all observations and failures are retained. Whole-process RSS includes imports and preparation.

First-use prepares fresh parsed objects outside timing and can compare different eager/lazy preparation. Repeated-use warms the operation once on each object; parse-inclusive times fresh parsing and evaluation. Return objects are retained during timing; normalization and semantic checks occur outside timing. Every block validates output determinism and wheel/extension identity.

An external checkpoint wrapper leaves frozen benchmark workers, timers, output checks and selection policy unchanged. It reserves 16 MiB of owned storage and on ENOSPC releases that reserve and retries writing the same computed checkpoint for at most 60 seconds. It never removes or remeasures a performance observation. Wrapper logs retain any such event. Matrix completion and any checkpoint events are recorded in the packet; all initial and selected adaptive matrices are retained.

Initial comparator strict wins: first_use: 12/27, repeated_use: 23/27, parse_inclusive: 27/27. COSMolKit has six unequal ring counts at offset 0; the corresponding modes remain unranked. Every remaining deficit is shown below.

## Initial matrix: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.996 [0.991, 1.010]; 7/20 | 1.494 [1.485, 1.531]; 19/20 | 0.551 [0.532, 0.559]; 0/20 |
| tpsa/repeated_use | 1.000 [0.985, 1.034]; 10/20 | 9.768 [9.620, 10.087]; 20/20 **strict win** | 1.726 [1.665, 1.877]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.001 [0.995, 1.027]; 11/20 | 13.121 [13.070, 13.259]; 20/20 **strict win** | 12.428 [12.302, 12.582]; 20/20 **strict win** |
| labute_asa/first_use | 0.998 [0.990, 1.005]; 8/20 | 0.737 [0.721, 0.752]; 0/20 | 0.359 [0.346, 0.372]; 0/20 |
| labute_asa/repeated_use | 1.001 [0.987, 1.061]; 11/20 | 2.708 [2.628, 2.808]; 20/20 **strict win** | 1.521 [1.424, 1.623]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.004 [0.995, 1.017]; 11/20 | 12.377 [12.265, 12.477]; 20/20 **strict win** | 11.860 [11.739, 12.059]; 20/20 **strict win** |
| num_rings/first_use | 1.004 [0.997, 1.013]; 13/20 | 0.037 [0.036, 0.040]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 0.991 [0.983, 1.043]; 9/20 | 1.727 [1.687, 1.795]; 20/20 **strict win** | unranked (output mismatch) |
| num_rings/parse_inclusive | 1.006 [0.997, 1.022]; 12/20 | 11.638 [11.539, 11.708]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 1.008 [0.999, 1.013]; 14/20 | 3.130 [3.113, 3.149]; 20/20 **strict win** | 2.074 [2.047, 2.145]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 1.010 [1.005, 1.033]; 15/20 | 5.208 [5.156, 5.249]; 20/20 **strict win** | 3.252 [3.190, 3.341]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.003 [0.996, 1.024]; 12/20 | 8.689 [8.582, 8.765]; 20/20 **strict win** | 7.683 [7.615, 7.785]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.004 [0.997, 1.010]; 12/20 | 1.280 [1.256, 1.297]; 20/20 **strict win** | 0.834 [0.827, 0.856]; 3/20 |
| morgan2_chiral/repeated_use | 1.010 [1.006, 1.018]; 17/20 | 5.058 [5.007, 5.103]; 20/20 **strict win** | 3.161 [3.103, 3.317]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.002 [0.997, 1.008]; 11/20 | 4.490 [4.479, 4.528]; 20/20 **strict win** | 3.965 [3.925, 4.001]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.008 [1.003, 1.013]; 15/20 | 3.398 [3.356, 3.424]; 20/20 **strict win** | 2.271 [2.239, 2.302]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.009 [1.004, 1.017]; 17/20 | 4.868 [4.807, 4.970]; 20/20 **strict win** | 3.167 [3.123, 3.274]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.004 [0.997, 1.022]; 12/20 | 8.033 [7.926, 8.139]; 20/20 **strict win** | 7.028 [6.987, 7.134]; 20/20 **strict win** |
| maccs/first_use | 1.000 [0.996, 1.003]; 9/20 | 1.886 [1.877, 1.896]; 20/20 **strict win** | 2.174 [2.165, 2.208]; 20/20 **strict win** |
| maccs/repeated_use | 1.002 [0.994, 1.003]; 13/20 | 1.950 [1.941, 1.959]; 20/20 **strict win** | 2.250 [2.234, 2.270]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.000 [0.996, 1.003]; 10/20 | 2.112 [2.104, 2.124]; 20/20 **strict win** | 2.388 [2.374, 2.400]; 20/20 **strict win** |

## Initial matrix: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.028 [1.022, 1.039]; 18/20 | 0.451 [0.445, 0.457]; 0/20 | 0.120 [0.117, 0.122]; 0/20 |
| tpsa/repeated_use | 1.006 [0.983, 1.017]; 11/20 | 9.857 [9.520, 10.009]; 20/20 **strict win** | 1.288 [1.256, 1.336]; 19/20 |
| tpsa/parse_inclusive | 1.026 [1.019, 1.035]; 17/20 | 5.241 [5.150, 5.308]; 20/20 **strict win** | 5.421 [5.349, 5.446]; 20/20 **strict win** |
| labute_asa/first_use | 0.999 [0.991, 1.010]; 8/20 | 0.210 [0.206, 0.214]; 0/20 | 0.086 [0.082, 0.088]; 0/20 |
| labute_asa/repeated_use | 1.011 [0.992, 1.020]; 12/20 | 2.597 [2.540, 2.648]; 19/20 | 1.161 [1.128, 1.234]; 18/20 |
| labute_asa/parse_inclusive | 1.007 [0.996, 1.012]; 12/20 | 4.873 [4.782, 4.930]; 20/20 **strict win** | 5.217 [5.128, 5.275]; 20/20 **strict win** |
| num_rings/first_use | 0.998 [0.985, 1.014]; 9/20 | 0.098 [0.097, 0.103]; 0/20 | 0.082 [0.075, 0.087]; 0/20 |
| num_rings/repeated_use | 1.023 [0.997, 1.039]; 12/20 | 1.670 [1.658, 1.743]; 20/20 **strict win** | 1.066 [1.038, 1.097]; 17/20 |
| num_rings/parse_inclusive | 1.003 [0.998, 1.022]; 12/20 | 11.797 [11.693, 11.970]; 20/20 **strict win** | 12.934 [12.808, 13.025]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.005 [0.999, 1.024]; 14/20 | 1.558 [1.548, 1.595]; 20/20 **strict win** | 0.997 [0.990, 1.015]; 9/20 |
| morgan2_2048/repeated_use | 1.006 [1.000, 1.014]; 14/20 | 4.742 [4.714, 4.776]; 20/20 **strict win** | 2.882 [2.864, 2.889]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.004 [0.999, 1.014]; 14/20 | 4.664 [4.586, 4.698]; 20/20 **strict win** | 4.482 [4.443, 4.516]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.001 [0.997, 1.017]; 10/20 | 1.605 [1.593, 1.620]; 20/20 **strict win** | 1.022 [1.016, 1.039]; 17/20 |
| morgan2_chiral/repeated_use | 1.010 [1.000, 1.023]; 14/20 | 4.749 [4.698, 4.875]; 20/20 **strict win** | 2.902 [2.888, 2.932]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.007 [1.001, 1.010]; 15/20 | 4.666 [4.618, 4.736]; 20/20 **strict win** | 4.515 [4.478, 4.537]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.003 [0.996, 1.008]; 12/20 | 1.800 [1.785, 1.816]; 20/20 **strict win** | 1.181 [1.171, 1.186]; 19/20 |
| morgan3_2048/repeated_use | 1.007 [1.004, 1.014]; 15/20 | 4.544 [4.498, 4.593]; 20/20 **strict win** | 2.880 [2.859, 2.917]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.009 [0.998, 1.010]; 13/20 | 4.580 [4.542, 4.603]; 20/20 **strict win** | 4.338 [4.312, 4.362]; 20/20 **strict win** |
| maccs/first_use | 0.998 [0.993, 1.001]; 8/20 | 1.639 [1.629, 1.649]; 20/20 **strict win** | 1.749 [1.741, 1.759]; 20/20 **strict win** |
| maccs/repeated_use | 1.003 [1.001, 1.005]; 17/20 | 1.736 [1.725, 1.746]; 20/20 **strict win** | 1.856 [1.853, 1.865]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.006 [1.003, 1.009]; 18/20 | 1.805 [1.793, 1.815]; 20/20 **strict win** | 1.934 [1.928, 1.942]; 20/20 **strict win** |

## Adaptive observations

The plan was saved before timing: an initial source-relative CI upper < 1 selects its entire operation for one separate 20-block follow-up in all modes. Initial and follow-up observations remain separate, without pooling or recursive retries. All new follow-up negatives remain recorded.

```json
{
  "5000": {
    "negative_initial_lanes": {},
    "operations": [],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  },
  "0": {
    "negative_initial_lanes": {},
    "operations": [],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  }
}
```

## Reproduction and status

[The evidence packet](2026-10-15-python-tpsa-zero-artifacts/README.md) contains both measured wheels, every initial/adaptive matrix, compressed exact output dumps, build/test logs and structural coverage, original implementation, exact patch, and all frozen source/input files. Archive readback checks all 10062 candidate hashes; baseline overrides and removal of the added test reconstruct all 10061 baseline hashes. The patch independently reconstructs changed files in a tiny Git repository. Hashes demonstrate source reconstruction, not bit-reproducibility of rebuilt wheels.

Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared Git metadata was not refreshed; this is not merged, published or hosted-CI-validated. Prior unresolved regressions against older source candidates, including V3000 reading and Morgan/MACCS observations, are not resolved by an immediate-baseline comparison. OpenBabel format lanes are outside this descriptor experiment. Separate-run win counts do not establish attributable algorithm gains.
