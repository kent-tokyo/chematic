# Morgan hashing and TPSA environment typing — 2026-10-19

**Decision: reject the combined source candidate.** Two fixed full-corpus primaries must both have unrounded paired 95% lower bounds above 1, and all 21 source operation/boundary medians must be at least 0.99. This is one-host source-wheel evidence, not a published-package, second-host or toolkit-wide claim.

| First-use source primary | Paired speedup | Paired 95% CI | Faster blocks |
|---|---:|---|---:|
| `morgan2_chiral` | 1.013350x | [1.009352, 1.015810] | 24/24 |
| `tpsa` | 1.031682x | [1.021296, 1.035423] | 24/24 |

The minimum source median is 0.989148273 in `num_rings/repeated_use`. There are 6 medians below 1 and 1 intervals entirely below 1. The largest observed speed decrease is 1.085%. Every lane, including non-wins and unequal-output cases, remains in the [63-comparison CSV](2026-10-19-morgan-tpsa-environments-lanes.csv) and [numeric packet](2026-10-19-morgan-tpsa-environments-numeric.json.gz).

## Primary comparison

Times are median milliseconds per full 10,000-row batch. Speedups use paired block ratios rather than ratios of aggregate medians.

| First-use operation | Baseline | Candidate | RDKit | COSMolKit |
|---|---:|---:|---:|---:|
| `morgan2_chiral` | 79.269146 | 78.281709 | 125.689500 | 82.877875 |
| `tpsa` | 3.704187 | 3.590729 | 18.869000 | 6.132459 |

| First-use operation | Comparator | Speedup | Paired 95% CI | Faster blocks | Strict external win |
|---|---|---:|---|---:|---|
| `morgan2_chiral` | rdkit | 1.603902x | [1.600272, 1.616004] | 24/24 | True |
| `morgan2_chiral` | cosmolkit | 1.059537x | [1.053530, 1.062300] | 24/24 | True |
| `tpsa` | rdkit | 5.248701x | [5.181135, 5.303068] | 24/24 | True |
| `tpsa` | cosmolkit | 1.704841x | [1.696011, 1.724586] | 24/24 | True |

A strict external win requires equivalent output, an unrounded lower bound above 1 and 24/24 faster blocks. A favorable median alone is insufficient.

- `rdkit`: 20/21 comparable operation/boundary lanes pass the strict external rule.
- `cosmolkit`: 17/18 comparable operation/boundary lanes pass the strict external rule.

COSMolKit ring counts differ on six normal-corpus rows, so all three ring-count boundaries remain unranked against it. Other external non-wins and every raw duration remain in the packet. Seven selected operations on an exposed corpus cannot establish universal superiority. The [published v1.1.8 ASA record](2026-10-19-v118-published-labute-asa-cosmolkit.md) separately retains its COSMolKit first-use tie; this source experiment does not supersede that package result.

## Negative source observations

| Lane | Ratio | Paired 95% CI |
|---|---:|---|
| `labute_asa/parse_inclusive` | 0.990993353 | [0.988464337, 0.998748771] |

## Changes and correctness

- Morgan connectivity hashing folds the same five components directly into the existing 32-bit hash, then appends the same ring flag when needed. The stack-backed component builder and traversal are removed; all values, integer casts and hash order remain identical.
- TPSA only computes three-membered-ring membership for the two neutral nitrogen environments and one neutral oxygen environment that consume it. Other table cases and fallbacks ignore that flag. Their unnecessary neighbor-pair bond checks are skipped; the original checks remain for every ring-sensitive type. Parser preparation and cache rules follow the source baseline.
- 352 fingerprint tests pass on the unchanged fingerprint source, including the frozen component-vector oracle over every atom of all 14,072 normal/stress rows with actual and empty ring flags. The new chemistry variant passes 989 library tests (1 existing ignore), three targeted frozen-TPSA rechecks and modified-crate Clippy. Frozen old environment typing covers normal/stress raw and Kekule forms; targeted amine/ether/charged/isotopic-H/P/S cases also match.
- All 56,288 complete R2/R3 chiral/nonchiral fingerprints equal both source baseline and RDKit. All 42,216 scalar values/refusal states remain source-exact. TPSA agrees with RDKit on all 14,072 rows. Refusals are states, not numeric matches.

Retained non-equal scalar states from the broad audit:

| Arm | Corpus / operation | Counts other than numeric equality |
|---|---|---|
| candidate | `stress-v1.smi/labute_asa` | both_unavailable=3 |
| candidate | `stress-v1.smi/num_rings` | mismatch=4 |
| cosmolkit | `rdkit-js-browser-10k-v1.smi/num_rings` | mismatch=6 |
| cosmolkit | `stress-v1.smi/labute_asa` | both_unavailable=3, source_unavailable=6 |
| cosmolkit | `stress-v1.smi/num_rings` | mismatch=4, source_unavailable=6 |
| cosmolkit | `stress-v1.smi/tpsa` | source_unavailable=6 |

## Earlier rejected experiment and diagnostic

The [standalone direct-hash epoch](2026-10-19-morgan-direct-connectivity.md) is rejected: chiral Morgan improves 1.011160x but first-use TPSA falls to 0.977690x, below the same 0.99 floor. All 24 blocks and 63 comparisons remain available. The combined candidate adds a genuinely different TPSA change and a second improvement primary. No earlier elapsed samples are reused or pooled. The baseline wheel is reused as an immutable artifact.

Before formal timing, eight balanced two-source diagnostic blocks checked only TPSA and chiral Morgan. Their samples remain separately tagged as diagnostic-only in the packet and never decide adoption.

| Diagnostic lane | Source median ratio |
|---|---:|
| `tpsa/first_use` | 1.021462182 |
| `tpsa/repeated_use` | 1.034163452 |
| `tpsa/parse_inclusive` | 1.007265946 |
| `morgan2_chiral/first_use` | 1.012323930 |
| `morgan2_chiral/repeated_use` | 1.056826815 |
| `morgan2_chiral/parse_inclusive` | 1.011425813 |

## Method and provenance

- Baseline source: `00596fb1ea5397524f98e3789c41f17ee8616fdc`. Candidate source: `373568afa5fb2ae5ce4b86994d365c11db3f4ea2`. Each native module was originally built from an independent empty target with native-inchi, release opt3, full LTO, one codegen unit, one build job, no incremental build and no stripping. Loaded extension hashes match the supplied wheels. The baseline artifact is reused; candidate build and all formal timings are new.
- Host: `macOS-27.0.1-arm64-arm-64bit-Mach-O`, CPython `3.13.6`, 10 logical CPUs. Installed pinned comparators are RDKit `2026.09.1` and COSMolKit `0.5.0`.
- Corpus: all 10,000 rows at offset 0 of `rdkit-js-browser-10k-v1.smi`, SHA-256 `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`. This is exposed development data, not a blind holdout.
- All 24 four-arm permutations occur once, shuffled with seed 20261019. Each arm occupies each position six times; pair precedence is balanced twelve/twelve. Arms run sequentially in fresh processes. All workers use `PYTHONHASHSEED=20261019`; the preceding rejected epoch used its inherited Python hash seed. Epochs remain separate.
- Each operation/boundary cell times one full 10,000-row batch per block. First-use parses fresh objects outside timing; repeated-use primes fresh objects outside timing; parse-inclusive puts parsing and the actual operation in one loop. Native return objects are retained and normalized outside timing. No row, error, negative observation or slow block is filtered.
- Pointwise paired 95% bootstrap intervals use 10,000 fixed-seed resamples. Both source primaries and all 21 median floors use unrounded values. Other host activity was not isolated; block-start load averages are retained. This limits transfer to other workloads, hosts and artifacts.
- Source/runtime/script seals are checked before and after the complete epoch. Local source snapshots, modules, full molecular outputs and logs remain local. The public packet exports numeric timing/accounting, output digests and selected provenance. Required benchmark CI recomputes its decisions and exercises corruption controls; it does not repeat this local speed experiment.

Parse-only and whole-process peak RSS (including imports and all operations) are separate descriptive observations:

| Arm | Median parse-only ms / 10,000 | Median process peak RSS MiB |
|---|---:|---:|
| baseline | 64.196604 | 168.922 |
| candidate | 64.199167 | 168.828 |
| rdkit | 403.078583 | 372.109 |
| cosmolkit | 419.464625 | 465.125 |

## Reproduction

```text
env PYTHONHASHSEED=20261019 python scripts/bench_python_speed_frontier.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --limit 10000 --blocks 24 --balanced-permutations \
  --ops tpsa,labute_asa,num_rings,morgan2_2048,morgan2_chiral,morgan3_2048,maccs \
  --base-revision 00596fb1ea5397524f98e3789c41f17ee8616fdc \
  --baseline-wheel <baseline-wheel> --baseline-path <baseline-site> \
  --candidate-wheel <candidate-wheel> --candidate-path <candidate-site> \
  --output <raw-report.json>
python scripts/check_python_speed_frontier_packet.py \
  benchmarks/2026-10-19-morgan-tpsa-environments-numeric.json.gz
```
