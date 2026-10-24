# RDKit search and adaptive 3D follow-up

Date: 2026-10-25. Source candidate based on `origin/main`
`b2d6d3a4f9e3b100992d2fab86aeaee32df3d26a`; comparator RDKit 2026.09.1
on macOS arm64. These are source-built results on existing development
corpora. The records retain exact comparator, corpus, script, and binary
hashes. They do not establish toolkit-wide or published-package superiority.

## Aromatic triple-bond ring in RDKit-parity search

The remaining two boolean SMARTS differences in the 4,072-row stress corpus
came from one molecule containing a six-carbon aromatic cycle with an
explicit triple bond. Global aromaticity re-perception removed the other
five aromatic bonds when another component of the molecule needed full
re-perception. RDKit kept the aromatic atom and bond flags. The new narrow
rule retains those flags for an isolated six-carbon, six-bond aromatic cycle
with one triple bond; it leaves the triple bond order intact. Direct regression
tests cover the simple ring and the multi-component stress molecule.

The [search packet](../validation/results/rdkit-search-aromatic-triple-2026-10-25.json)
compares RDKit `HasSubstructMatch` with 25 SMARTS queries on all 4,072
RDKit-accepted stress rows. The RDKit-parity mode now agrees on
**101,800/101,800** booleans, with zero refusals. The CLI binary SHA-256 is
`68beac9036172dbd4a96e0281ea5a4d504021e468231aa26bedaad22be96ad5e`.
The default native search mode has 96,225/101,800 agreements; its different
aromaticity semantics remain visible in the packet. Neither count covers
all valid SMARTS or molecules.

After this change, the existing fingerprint/stereo gate was rerun with RDKit
2026.09.1. It found 95,360/95,360 exact 2048-bit fingerprint cells across
five algorithms and three fixed corpora, plus zero mismatches on 324 stereo
search pairs and 1,653 stereo self-searches. Its output is byte-identical to
the [prior machine packet](2026-10-25-rdkit-core-usr-accuracy.json) (SHA-256
`c09a14d33071ef398553ef7c8c84e25546c0e5a550a870523879b111dd576fc9`).
This is a regression check on fixed data, not proof of superiority in search
or fingerprint accuracy.

## Adaptive MMFF94 retry and stereo

Adaptive and gradient-rescue retries now accept replacement coordinates only
when all declared stereochemical constraints can be evaluated and remain
satisfied, in addition to their existing finite-energy and non-increase
rules. A focused test rejects mirrored and planar geometry at a declared
tetrahedral center. The same gate is used by the existing quality profile.

The [12-row comparison](../validation/results/etkdg-adaptive-stereo-safe-chembl-hash24-35-2026-10-25.json)
uses the next 12 eligible ChEMBL SMILES in SHA-256 order, positions 24–35.
Both engines initially generate 20 ETKDGv3 seed-42 candidates and run up to
200 MMFF94 iterations per candidate. CheMatic's adaptive arm alone continues
its lowest-energy unfinished candidate for up to 400 more iterations. The
metric is the lowest **converged** MMFF94 energy; a difference within 1e-6
kcal/mol is a tie.

| Arm | Comparable | CheMatic lower | RDKit lower | Ties |
| --- | ---: | ---: | ---: | ---: |
| CheMatic compatibility, 20 initial candidates | 12 | 1 | 1 | 10 |
| CheMatic adaptive, plus one bounded retry | 12 | 3 | 0 | 9 |

The sole RDKit-lower compatibility row, `chembl_corpus_4885`, moved from
50.253775 to 35.897340 kcal/mol after the bounded retry; RDKit's best
converged result was 44.431771 kcal/mol. RDKit also found a 35.898842
kcal/mol **unconverged** conformer for this row, which is excluded from the
stated converged-energy metric. The adaptive arm improved its own converged
minimum on 3/12 rows. The scorer checks atom-indexed graph and bond orders,
independently recomputes MMFF94 energy on returned CheMatic coordinates, and
checks declared tetrahedral and double-bond stereo. The largest absolute
same-coordinate energy difference is below 6e-14 kcal/mol.

The [opt-in torsion profile on the same 12 rows](../validation/results/etkdg-torsion-search-chembl-hash24-35-2026-10-25.json)
found five CheMatic-lower converged minima and seven ties against 60 RDKit
candidates, with no RDKit-lower row. That profile can spend 60–80 initial
embeddings, additional retry iterations, and torsion minimizations. It is a
separate, higher-compute search-quality observation. Neither experiment
measures experimental geometry, physical force-field accuracy, conformer
diversity, throughput, or equal-work superiority.

## Reproduce

Run with RDKit 2026.09.1 in the selected Python environment:

```sh
cargo build --locked -p chematic-cli
cargo build --locked -p chematic-3d --example etkdg_quality_seed_budget_probe
python scripts/validation/rdkit_cli_search_differential.py \
  --binary target/debug/chematic \
  --output /tmp/rdkit-search-aromatic-triple.json
python scripts/validation/rdkit_cli_stereo_fingerprint_gate.py \
  --binary target/debug/chematic \
  --stereo-self-corpus scripts/fixtures/descriptor_census_corpus.smi \
  --fingerprint-corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --fingerprint-corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --fingerprint-corpus scripts/fixtures/descriptor_census_corpus.smi
python scripts/benchmarks/etkdg_adaptive_stereo_safe_probe.py \
  --binary target/debug/examples/etkdg_quality_seed_budget_probe \
  --chembl-start 24 --output /tmp/etkdg-adaptive-stereo-safe.json
cargo build --locked -p chematic-3d --example etkdg_torsion_search_probe
python scripts/benchmarks/etkdg_torsion_search_probe.py \
  --binary target/debug/examples/etkdg_torsion_search_probe \
  --tier chembl --chembl-start 24 \
  --output /tmp/etkdg-torsion-search-chembl-hash24-35.json
```

The reported source binary SHA-256 values distinguish results from any later
rebuild. The benchmark scripts use fixed selection rules and fail on graph,
energy, or declared-stereo discrepancies.
