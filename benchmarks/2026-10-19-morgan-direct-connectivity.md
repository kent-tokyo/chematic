# Direct Morgan connectivity hashing — 2026-10-19

**Decision: reject the source candidate.** First-use chiral radius-2 Morgan on all 10,000 fixed rows has source speedup 1.011160x, paired 95% interval [1.005350, 1.014816], and 19/24 faster blocks.

This is one-host source-wheel evidence, not a published-package, second-host or toolkit-wide claim. The preregistered primary is the full-corpus first-use chiral Morgan lane; adoption also requires all 21 source operation/boundary median ratios to be at least 0.99. The primary must have a 95% lower bound above 1. Unrounded values decide both gates.

## Source results

The minimum source ratio is 0.977689664 in `tpsa/first_use`. There are 7 ratios below 1 and 4 intervals entirely below 1; the largest observed speed decrease is 2.231%. All observations remain in the [63-comparison CSV](2026-10-19-morgan-direct-connectivity-lanes.csv) and [numeric packet](2026-10-19-morgan-direct-connectivity-numeric.json.gz).

| First-use chiral Morgan | Median ms / 10,000 rows |
|---|---:|
| baseline | 79.205687 |
| candidate | 78.476208 |
| rdkit | 126.531312 |
| cosmolkit | 83.002813 |

| Comparator | Paired speedup | Paired 95% CI | Faster blocks | Strict win |
|---|---:|---|---:|---|
| baseline | 1.011160x | [1.005350, 1.014816] | 19/24 | False |
| rdkit | 1.605929x | [1.595746, 1.619714] | 24/24 | True |
| cosmolkit | 1.054173x | [1.045795, 1.063985] | 22/24 | False |

Source lanes with an entirely negative interval:

| Lane | Ratio | 95% CI |
|---|---:|---|
| `tpsa/first_use` | 0.977689664 | [0.967888435, 0.990089120] |
| `maccs/first_use` | 0.995100186 | [0.993129825, 0.996875476] |
| `maccs/repeated_use` | 0.992857777 | [0.988197368, 0.994953178] |
| `maccs/parse_inclusive` | 0.990410450 | [0.984853738, 0.996580956] |

## External scope

- `rdkit`: 19/21 comparable operation/boundary lanes meet the strict lower-bound-above-1 and 24/24-faster rule.
- `cosmolkit`: 14/18 comparable operation/boundary lanes meet the strict lower-bound-above-1 and 24/24-faster rule.

Every unequal-output lane is unranked. The numeric packet retains all normal-corpus mismatches and denominators, including COSMolKit ring-count differences. First-use, repeated-use, and parse-inclusive results are distinct. Seven selected operations on a development corpus cannot establish universal superiority.

## Implementation and correctness

The runtime change folds the five connectivity-invariant components directly into the existing 32-bit hash, then conditionally appends the same ring flag. It removes a stack-backed SmallVec builder and its traversal; component values, integer casts, hash order and chemical algorithms are preserved. Parser preparation and cache behavior follow the source baseline.

- 352 fingerprint library tests pass, including an independent copy of the old component-vector implementation on all atoms of 14,072 normal/stress inputs, with both actual and empty ring flags.
- All 56,288 complete Morgan fingerprints (radii 2/3, chirality off/on) are byte-exact against both source baseline and RDKit. All 14,072 rows are retained; there are no source/RDKit refusals in this lane.
- Modified-crate Clippy passes. Six benchmark-method tests cover output/tolerance boundaries and the full permutation schedule. A boundary test proves that unrounded primary/control values, rather than three-decimal display values, decide adoption. Packet corruption checks and hosted CI are separate integration gates.

## Method and provenance

- Baseline source: `00596fb1ea5397524f98e3789c41f17ee8616fdc`. Candidate source: `d7d3e7de6ba7e2f3e47bd3f4e567a29126e2c4d7`. Both native modules are built from independent empty targets using native-inchi, release opt3, full LTO, one codegen unit, one build job, no incremental build and no stripping. Installed native hashes are verified against the supplied wheels. Wheel/native/source hashes and the registration hash are retained in the numeric packet.
- Host: `macOS-27.0.1-arm64-arm-64bit-Mach-O`, CPython `3.13.6`, 10 logical CPUs. This experiment uses installed pinned RDKit `2026.09.1` and COSMolKit `0.5.0`.
- Corpus: every row of `rdkit-js-browser-10k-v1.smi`, SHA-256 `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`, offset 0, limit 10,000. The corpus is exposed development data, not a blind holdout.
- All 24 four-arm permutations occur once, shuffled with seed 20261019. Each arm occupies each position six times; each pair occurs in both precedence orders twelve times. Arms run sequentially in fresh processes. Every operation/boundary cell times one full 10,000-row batch per block.
- First-use parses fresh objects outside the timer; repeated-use primes fresh objects outside the timer; parse-inclusive measures parsing and the actual operation in the same loop. Native API returns are retained and normalized outside timing. No output, error, negative observation or slow block is filtered.
- Pointwise paired bootstrap intervals use 10,000 resamples and fixed comparator-specific seeds. All 21 source floors and the primary are recomputed from raw durations without rounding. Other host activity was not isolated; block-start load averages are retained. This limits transfer to other workloads, hosts and exact packages.
- Local source snapshots, installed modules, full molecular output values and build logs are retained locally. The committed packet contains numeric timing/accounting, output digests and explicitly selected provenance only. Hosted CI recomputes its accounting and decisions; it does not reproduce these local speed measurements.

Whole-process peak RSS includes imports and all operations; it is descriptive, not a per-operation ranking. Parse-only time is reported independently of the 21 operation/boundary gates.

| Arm | Median parse-only ms / 10,000 | Median process peak RSS MiB |
|---|---:|---:|
| baseline | 64.139062 | 168.898 |
| candidate | 64.981750 | 169.086 |
| rdkit | 405.075000 | 372.352 |
| cosmolkit | 420.266604 | 465.789 |

## Reproduction

```text
python scripts/bench_python_speed_frontier.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --limit 10000 --blocks 24 --balanced-permutations \
  --ops tpsa,labute_asa,num_rings,morgan2_2048,morgan2_chiral,morgan3_2048,maccs \
  --base-revision 00596fb1ea5397524f98e3789c41f17ee8616fdc \
  --baseline-wheel <baseline-wheel> --baseline-path <baseline-site> \
  --candidate-wheel <candidate-wheel> --candidate-path <candidate-site> \
  --output <raw-report.json>
python scripts/check_python_speed_frontier_packet.py \
  benchmarks/2026-10-19-morgan-direct-connectivity-numeric.json.gz
```
