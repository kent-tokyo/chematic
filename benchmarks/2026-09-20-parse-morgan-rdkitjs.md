# Parse + compatible Morgan vs official RDKit.js — 2026-09-20

## Scope

This is a **source-candidate** comparison for merge commit
[`7d98dcd3`](https://github.com/kent-tokyo/chematic/commit/7d98dcd33e7f587f93fed6b015da807bacaf9e8e),
not a measurement of a published chematic package. The comparator is pinned to
`@rdkit/rdkit@2026.3.6` (runtime `2026.03.6`).

Each timed row parses a SMILES, computes a radius-2 2048-bit compatible Morgan
fingerprint, packs and consumes the same 256-byte LSB-first output, then
releases it. RDKit's bit-string packing is inside its timed arm. This excludes
download-to-ready, initialization, prepared-object reuse, similarity search,
and browser memory.

## Correctness boundary

The fixed 10k corpus has 9,999 supported molecules and one unchanged typed
refusal: a degree-10 Fe(II) coordination structure. Every supported row had
exact configured bits in Chromium, Firefox, and WebKit. An independent ChEMBL
5k Chromium gate is 5,000/5,000 exact. The candidate also passes all three
Rust/Python/Node-WASM pairs on 5,000 molecules.

These are configured compatible-Morgan bits, not a claim that chematic's
native ECFP4 or every Morgan option is universally bit-identical to RDKit.

## Local macOS browser records

Fresh browser processes were used per arm and repetition. Values are median
milliseconds per molecule; the speedup is RDKit divided by chematic.

| Browser / corpus | Repetitions | chematic | RDKit.js | Median speedup |
|---|---:|---:|---:|---:|
| Chromium / fixed 10k | 20 | 0.039839 | 0.105116 | 2.63x |
| Firefox / fixed 10k | 10 | 0.314531 | 0.647615 | 2.06x |
| WebKit / fixed 10k | 10 | 0.037204 | 0.110011 | 2.96x |
| Chromium / independent ChEMBL 5k | 20 | 0.022650 | 0.139910 | 6.17x |

Raw records: [Chromium 10k](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-chromium-10k.json),
[Firefox 10k](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-firefox-10k.json),
[WebKit 10k](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-webkit-10k.json),
and [ChEMBL 5k](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-chromium-chembl-5k.json).

## GitHub-hosted confirmation

The same fixed-10k operation ran in the three-browser
[GitHub Actions gate](https://github.com/kent-tokyo/chematic/actions/runs/35484834092).
The gate uses ten paired repetitions and requires the paired log-speedup 95%
lower bound to exceed 1.0. Its independently recomputed results were:

| Browser | Geometric mean speedup | 95% lower bound |
|---|---:|---:|
| Chromium | 3.14x | 3.11x |
| Firefox | 2.36x | 2.36x |
| WebKit | 3.73x | 3.67x |

The workflow stores raw JSON as CI artifacts and calls
[`scripts/check_parse_morgan_rdkit_speed_gate.py`](../scripts/check_parse_morgan_rdkit_speed_gate.py).
It is a source-build confirmation, not a release-package benchmark. A future
published package must be rerun through the same contract before it receives a
release-specific speed claim.

## Reproduce

Build a fresh Web-WASM candidate and run the isolated runner with the pinned
package. The CI workflow is the canonical three-browser command set:
[`parse-morgan-rdkitjs-gate.yml`](../.github/workflows/parse-morgan-rdkitjs-gate.yml).

For a local JSON record, use
[`scripts/bench_browser_wasm_vs_rdkit_isolated.py`](../scripts/bench_browser_wasm_vs_rdkit_isolated.py)
with the fixed `rdkit-js-browser-10k-v1.smi` corpus, then validate it with
`python3 scripts/check_parse_morgan_rdkit_speed_gate.py --input RESULT.json`.
