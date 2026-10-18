# RDKit reaction, Cleanup, and labelled R-group source check

Date: 2026-10-18

CheMatic base: `540086f978612df1e21baa731c0884dea7441f46` plus the changes in this branch

Oracle: RDKit `2026.09.1`

## Cleanup

The source-built `cleanup_rdkit_dump` example was run over all 4,072 rows of
`validation/cosmolkit_comparison/corpora/stress-v1.smi`. RDKit accepted 4,069
Cleanup operations on inputs with the same parsed canonical SMILES. CheMatic
matched **3,885/4,069 (95.48%)**, leaving 184 differences, including seven
CheMatic output errors. RDKit's default `MetalDisconnector` element scope was
used for metal bond selection. The pre-change local extension matched
3,836/4,069 on the same compatible row set: 49 differences were fixed and no
previously matching rows regressed. This is a parity measure, not an
independent measure of chemical correctness.

Reproduce from this branch:

```sh
cargo build -p chematic-chem --example cleanup_rdkit_dump --offline
target/debug/examples/cleanup_rdkit_dump \
  < validation/cosmolkit_comparison/corpora/stress-v1.smi \
  > /tmp/chematic-cleanup-source.jsonl
python scripts/rdkit_cleanup_source_gate.py \
  --dump /tmp/chematic-cleanup-source.jsonl \
  --baseline benchmarks/2026-10-18-rdkit-cleanup-before-metal-scope.json \
  --output /tmp/chematic-cleanup-source-gate.json
```

The [machine-readable source result](2026-10-18-rdkit-cleanup-source-gate.json)
contains row-level differences, input/corpus/source hashes, and the 49 fixed
indices. The [pre-change extension result](2026-10-18-rdkit-cleanup-before-metal-scope.json)
records its binary hash. The extension predates the present source changes;
the candidate measurement uses a compiled Rust example.

## Reaction application

The existing 83-case graph, atom-origin, and template-map gate was rerun with
the source-built CheMatic `1.1.7` extension and RDKit `2026.09.1`: 80 cases had exact
semantic agreement and three inputs were invalid in both engines. Separately,
the pinned opposite-enantiomer case produced one raw product set in RDKit while
CheMatic refused it, as required by that case's stereo contract. The
[83-case result](2026-10-18-rdkit-reaction-83-source-gate.json) and
[stereo safety result](2026-10-18-rdkit-reaction-stereo-safety.json) are bounded
checks; they do not measure reaction yield, selectivity, or broad SMIRKS
precision/recall. Reaction source was unchanged in this branch.

## Labelled R-group decomposition

An explicit mapped core site without a substituent used to omit an R column
even when other input molecules established that column. RDKit's
`RGroupDecompose` yields a hydrogen group in that situation. The source fix
adds the hydrogen group only to an established column and gives the rows a
matching labelled core. A sole unsubstituted core still has no R column. The
Rust regression covers phenyl, cyclohexyl, oxygen, and nitrogen core sites.
The source-built Python extension (SHA-256
`8344e5e29734502ecad68177e6002fddcadc09281d396409e095c9958d8db221`)
matched RDKit's core and R1 graphs with the same labels on all four
two-molecule cases in the
[machine-readable R-group gate](2026-10-18-rdkit-rgroup-hydrogen-source-gate.json).
The five labelled R-group Python tests also passed. This
change does not add multi-core alignment, nonterminal attachment support, or
an unbounded match search.

## Remaining work

Cleanup still differs on 184 scored rows. The labelled R-group surface is
single-core and bounded. The 83-case reaction gate is a selected compatibility
set. None of these measurements establishes general accuracy above RDKit.
