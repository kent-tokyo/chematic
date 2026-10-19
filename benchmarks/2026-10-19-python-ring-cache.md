# Morgan, TPSA and Python ring-count caches — 2026-10-19

**Decision: adopt this combined source candidate.** The preregistered three primary lanes must each have an unrounded paired 95% lower bound above 1; all 21 source operation/boundary medians must be at least 0.99. This is one-host source-wheel evidence, limited to seven selected operations and an exposed development corpus.

| Source primary | Paired speedup | Paired 95% CI | Faster blocks |
|---|---:|---|---:|
| `morgan2_chiral/first_use` | 1.007530619x | [1.004598281, 1.010795465] | 22/24 |
| `tpsa/first_use` | 1.166398766x | [1.157876831, 1.206684777] | 24/24 |
| `num_rings/repeated_use` | 1.048943389x | [1.021057357, 1.064191510] | 23/24 |

Minimum source median: 0.996752074 in `maccs/parse_inclusive`. The [CSV](2026-10-19-python-ring-cache-lanes.csv) retains all 63 comparisons. The [numeric packet](2026-10-19-python-ring-cache-numeric.json.gz) retains all 24 blocks, every raw nanosecond duration, output digests, exact intervals, source floors, pilot observations and the decision.

## External comparisons

A strict external win requires equivalent output, an unrounded paired lower bound above 1 and 24/24 faster blocks. Source adoption uses its separate preregistered rule.

- rdkit: 20/21 comparable lanes meet the strict external rule.

- cosmolkit: 17/18 comparable lanes meet the strict external rule.

External non-wins and ineligible comparisons:

| Lane | Comparator | Speedup | Paired 95% CI | Faster blocks |
|---|---|---:|---|---:|
| `labute_asa/first_use` | cosmolkit | 1.239826877x | [1.106474506, 1.282254477] | 20/24 |
| `num_rings/first_use` | rdkit | 1.913711267x | [1.899340170, 1.987492296] | 23/24 |
| `num_rings/first_use` | cosmolkit | unranked: unequal outputs | — | — |
| `num_rings/repeated_use` | cosmolkit | unranked: unequal outputs | — | — |
| `num_rings/parse_inclusive` | cosmolkit | unranked: unequal outputs | — | — |

COSMolKit ring counts differ on six normal-corpus rows, so its ring lanes are unranked. Selected-operation results cannot establish toolkit-wide superiority. The [published v1.1.8 ASA result](2026-10-19-v118-published-labute-asa-cosmolkit.md) separately retains a COSMolKit first-use tie; a source wheel cannot supersede that package result.

## Source controls

Every source lane below 1 is listed here, including uncertain decreases.

| Lane | Source speedup | Paired 95% CI |
|---|---:|---|
| `num_rings/parse_inclusive` | 0.999383217x | [0.996614725, 1.005448561] |
| `maccs/first_use` | 0.997525775x | [0.994913644, 1.000436319] |
| `maccs/repeated_use` | 0.997335422x | [0.995532286, 0.998230902] |
| `maccs/parse_inclusive` | 0.996752074x | [0.994977377, 1.000220031] |

## Changes and correctness

- Morgan hashing folds the same ordered components directly into the existing u32 hash instead of building and traversing a stack-backed component buffer. Casts, hash order and ring-flag behavior remain identical.
- TPSA computes three-membered-ring membership only for the two neutral nitrogen environments and one neutral oxygen environment that consume it. Other atom types skip the unused bond-pair checks. A frozen old environment builder verifies the optimization independently.
- A private numeric cache keeps four scalar cells and the ring-count atomic inline. The 19 type-erased cells move into a lazily allocated boxed array, initialized without running perception. The inline cache shrinks from 528 to 88 bytes on this arm64 build; the existing scalar offsets and numeric cache semantics stay unchanged, while ring count moves from offset 520 to 64. Used type-erased arrays are separately allocated; whole-process RSS is reported below. Payload computations still run outside initialization locks. Empty peeks and numeric-only use do not allocate payload storage; mutation and clone invalidation remain unchanged.
- The immutable Python molecule wrapper also memoizes its numeric ring count. Existing bounded SMILES preparation seeds the already computed count without extra eager work or a larger graph budget. Other graph-producing constructors start cold; clones share the immutable graph and may copy the numeric value. Mutable properties do not affect graph state; pickle omits the internal cache and restores it cold. The hot read avoids following the inner Arc to the core cache.
- Native library tests: chemistry 989 passed / 1 existing ignore, core 165 passed, fingerprints 352 passed, Python Rust unit tests 10 passed. 173 focused Python regression cases pass, including seven new cases covering property edits, pickle, cage and large cold graphs, editable snapshots and graph-producing transforms. Modified-library Clippy passes with warnings denied. The core suite includes the private layout contract and existing cache/invalidation checks.
- Fresh candidate output: all 56,288 R2/R3 chiral/nonchiral fingerprint cells equal both the frozen source baseline and RDKit; all 42,216 scalar value/refusal states equal the source baseline. All 14,072 TPSA values agree with RDKit within absolute tolerance 1e-9. No row is dropped.

Retained broad scalar audit states other than numeric equality:

| Arm | Corpus / operation | Counts |
|---|---|---|
| candidate | `stress-v1.smi/labute_asa` | both_unavailable=3 |
| candidate | `stress-v1.smi/num_rings` | mismatch=4 |
| cosmolkit | `rdkit-js-browser-10k-v1.smi/num_rings` | mismatch=6 |
| cosmolkit | `stress-v1.smi/labute_asa` | both_unavailable=3, source_unavailable=6 |
| cosmolkit | `stress-v1.smi/num_rings` | mismatch=4, source_unavailable=6 |
| cosmolkit | `stress-v1.smi/tpsa` | source_unavailable=6 |

## Preserved rejected epochs

1. [Standalone direct Morgan hash](2026-10-19-morgan-direct-connectivity.md): chiral Morgan first-use improved, but TPSA first-use 0.977689664 missed the 0.99 control floor. Rejected.
2. [Combined Morgan/TPSA typing](2026-10-19-morgan-tpsa-environments.md): both first-use primaries improved, but ring-count repeated-use 0.989148273 missed the same control floor. Rejected.

3. [Layout-only pilot](2026-10-19-morgan-tpsa-hot-scalars-pilot.md): eight diagnostic blocks showed TPSA repeated-use 0.928339526 and ring repeated-use 0.975118018. Rejected before formal timing.

4. [Compact core cache pilot](2026-10-19-morgan-tpsa-compact-cache-pilot.md): first-use scalar gains remain diagnostic; ring repeated-use 0.997178657 and ring parse-inclusive 0.989811712 prevented formal adoption.

This fifth candidate adds immutable Python-wrapper ring-count memoization to the compact-cache variant. It retains the same three required improvement primaries and all 21 control floors as the two rejected diagnostic candidates. No elapsed observations are pooled or reused. The immutable original baseline wheel is reused; its elapsed measurements are fresh. Eight balanced AB/BA pilot blocks cover all 21 cells before formal measurement. Pilot observations are separately tagged diagnostic-only and never decide formal adoption.

## Absolute batch and preparation observations

Median milliseconds per complete 10,000-row batch; paired ratios above are not ratios of these aggregate medians.

| Operation / boundary | Baseline | Candidate | RDKit | COSMolKit |
|---|---:|---:|---:|---:|
| `morgan2_chiral/first_use` | 79.418291 | 78.600250 | 126.243813 | 82.888937 |
| `tpsa/first_use` | 3.738000 | 3.195854 | 18.856458 | 6.252021 |
| `num_rings/repeated_use` | 0.312104 | 0.296792 | 0.546063 | 0.389104 |

Parse-only and whole-process RSS, including imports and every operation, are separate descriptive observations. Parse-only is outside the three source improvement primaries and 21 descriptor floors.

| Arm | Median parse-only ms / 10,000 | Median process peak RSS MiB |
|---|---:|---:|
| baseline | 64.208750 | 168.875 |
| candidate | 63.990874 | 168.906 |
| rdkit | 405.585958 | 372.156 |
| cosmolkit | 422.317500 | 465.547 |

Source parse-only paired ratio: 0.998572591, pointwise 95% CI [0.980755845, 1.005061309], 10/24 faster blocks. This observation is retained even if unfavorable.


## Method and provenance

- Frozen baseline source `00596fb1ea5397524f98e3789c41f17ee8616fdc`; frozen candidate source `4560db00b3056aec6197c051239da89497c90bf8`. The candidate native module is built from an empty independent target; the baseline module was originally built from a different empty target. Release opt3, full LTO, one codegen unit, one build job, native-inchi, no incremental build, debug0, no stripping. Loaded native hashes match the supplied wheel artifacts.
- Host `macOS-27.0.1-arm64-arm-64bit-Mach-O`, CPython `3.13.6`, 10 logical CPUs. Installed comparators: RDKit `2026.09.1` and COSMolKit `0.5.0`.
- Corpus: all 10,000 rows, offset 0, `rdkit-js-browser-10k-v1.smi`, SHA-256 `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`. It is exposed development data, not a blind holdout.
- All 24 permutations of baseline/candidate/RDKit/COSMolKit occur once, shuffled with seed 20261019. Every arm occupies every position six times; pair precedence is balanced twelve/twelve. Sequential fresh processes use PYTHONHASHSEED=20261019.
- Each timed cell is one full 10,000-row batch. First-use parses new objects outside timing; repeated-use primes fresh objects outside timing; parse-inclusive includes parsing and the operation. Native results are retained and normalized outside timing. Parse-only and process peak RSS remain separate observations in the packet.
- Paired pointwise 95% bootstrap intervals use 10,000 fixed-seed resamples. The three primaries and all 21 floors use unrounded numbers. Host activity is not isolated; block-start load averages are retained. No failure, slow block or unfavorable observation is filtered.
- Source and benchmark-script hashes are checked before and after timing. Local wheels, full Morgan/scalar correctness streams, source manifests and complete logs remain private. Public export is limited to numeric timings/accounting, digests and selected provenance. CI recomputes exact decisions and exercises corruption controls (39 local method/audit tests pass); the final packet uses `--require-adoption`, while archived rejected epochs do not. CI does not repeat the local timing experiment.

## Integration boundary

After this source-wheel epoch, main `5bd2d99c90a39ebac3c5131497061f8c249d92bf` was merged as `67878c937d2d8158b16de3de2702db17f7a02a10`. It brings NumPy fragment-counter support, published-package/channel records and dependency updates to smallvec 1.16.3, cc 1.8.0 and serde_json 1.0.154. The integrated source passes 1,516 local library tests (one existing chemistry ignore), modified-library Clippy with warnings denied and formatting. The index is checked against all 533 tracked benchmark records using path placeholders for sparse files; hosted CI checks the full checkout separately.

The integrated combination has not been timed again. All speed ratios, confidence intervals and adoption floors above belong to the frozen `00596fb1` / `4560db00` wheel pair, not the later main combination or a published package. The 173 focused Python tests also refer to the frozen candidate wheel; hosted Python CI validates the integrated source separately.

## Reproduction

Run `scripts/bench_python_speed_frontier.py` with the frozen source wheels, all seven operations, --limit 10000 --blocks 24 --balanced-permutations and PYTHONHASHSEED=20261019, as in the preceding combined experiment. Validate the public packet with `python scripts/check_python_speed_frontier_packet.py benchmarks/2026-10-19-python-ring-cache-numeric.json.gz --require-adoption`.
