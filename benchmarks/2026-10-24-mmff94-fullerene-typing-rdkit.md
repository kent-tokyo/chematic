# MMFF94 atom typing on fused fullerene cages

Date: 2026-10-24. Source candidate based on `origin/main` at `aeae3f53`.
Comparator: RDKit 2026.09.1, Python 3.13, macOS arm64. Both implementations
parse the same SMILES, add explicit hydrogens and compare numeric MMFF94 types
at the same atom index. The Rust probe also emits atomic numbers, so an
element/order mismatch is counted separately. RDKit molecules for which
`MMFFGetMoleculeProperties` returns no properties are excluded from type
comparison and counted.

| Corpus | Rows | Comparable | RDKit MMFF unsupported | Different atoms | Different rows |
| --- | ---: | ---: | ---: | ---: | ---: |
| Browser 10k | 10,000 | 9,796 | 204 | 0 | 0 |
| ChEMBL 5k | 5,000 | 4,990 | 10 | 0 | 0 |

The [browser row record](../validation/results/mmff94-rust-rdkit-browser10k-2026-10-24.json)
and [ChEMBL row record](../validation/results/mmff94-rust-rdkit-chembl5k-2026-10-24.json)
give corpus and source binary SHA-256 hashes, host, versions, status counts,
and difference categories. The unmodified 1.2.1 wheel measured immediately
before this source fix had 72 differing heavy-atom types across exactly three
browser inputs (indices 3837, 3976 and 4004), and zero differences on ChEMBL.
Every browser difference was fullerene carbon type 2 in RDKit versus type 37
in CheMatic. The source diagnostic before the fix gave 48 MMFF-aromatic atoms
on each cage; after the fix it gives RDKit's 24.
The CI source-wheel census against its pinned RDKit 2026.03.6 comparator
also produced zero atom-type differences on the 9,796 comparable browser
inputs; its exact-count gate now requires zero rather than the former 72.

The cause was in the RDKit-compatible Kekulé walk: CheMatic stopped expanding
neighbors once the current atom had spent its double-bond candidacy. RDKit
continues that walk even when the current atom cannot receive a double bond.
On these fused cages, the different traversal selected a matching whose
MMFF aromatic re-perception covered too many rings. The fix retains the
candidate check for forming a double bond and always queues unvisited
neighbors. The three fullerene structures are now regression fixtures.

## Reproduce

With RDKit 2026.09.1 installed in the Python interpreter:

```sh
cargo build --locked -p chematic-ff --example mmff94_type_probe
python3 scripts/diagnostics/mmff94_rust_rdkit_type_census.py \
  --binary target/debug/examples/mmff94_type_probe \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --output /tmp/mmff94-browser10k.json --require-exact
python3 scripts/diagnostics/mmff94_rust_rdkit_type_census.py \
  --binary target/debug/examples/mmff94_type_probe \
  --corpus scripts/fixtures/chembl_accuracy_corpus_4999.smi \
  --output /tmp/mmff94-chembl5k.json --require-exact
cargo test -p chematic-perception -p chematic-ff --lib
```

This result establishes atom-type agreement for these two fixed corpora. It
does not measure MMFF parameter coverage, energy, gradients, minimization,
experimental geometry or conformer search, so it does not establish overall
3D or force-field superiority over RDKit.
