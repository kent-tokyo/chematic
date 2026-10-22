# RDKit automatic R-group columns and shared-site fragments

Date: 2026-10-22. This is a source-built Python extension based on `654993b8`, compared with RDKit 2026.09.1 on macOS arm64. The exact source hash, extension hash, inputs, and both output tables are in the [machine record](2026-10-22-rdkit-rgroup-auto-columns.json). The comparator is [`scripts/rdkit_rgroup_auto_label_gate.py`](../scripts/rdkit_rgroup_auto_label_gate.py).

## Result

Four selected batches now match RDKit on all **16/16 rows** after canonicalizing each core and R-group graph while retaining atom-map labels. The source before this change matched 2/5 benzene-position rows and 0/3 pyridine rows; the `CC` and `CO` batches failed with `MultipleAttachments` instead of producing four rows each. The existing mapped-hydrogen gate also passes **16/16** with the new extension.

The change numbers only automatic sites established by the batch, in reverse core-query atom order, while preserving explicit labels. It also collects multiple substituents at one unlabelled site into a single R column. Disconnected components share the wildcard label; connected components are emitted once with all their boundary wildcards. Ambiguous multiple cuts at an explicitly mapped or stereogenic core atom still fail closed.

## Reproduce

```sh
cd crates/chematic-py
CARGO_INCREMENTAL=0 CARGO_TARGET_DIR=/private/tmp/chematic-rgroup-target maturin build --profile dev --no-default-features --features extension-module --interpreter python3 --offline
cd ../..
PYTHONPATH=<source-wheel-install>:<RDKit-2026.09.1-install> python3 scripts/rdkit_rgroup_auto_label_gate.py --output /private/tmp/rgroup-auto-columns.json
PYTHONPATH=<source-wheel-install>:<RDKit-2026.09.1-install> python3 scripts/rdkit_rgroup_hydrogen_gate.py --output /private/tmp/rgroup-hydrogen.json
```

Rust validation: 1,030 `chematic-chem` unit tests passed, one ignored; package Clippy with `-D warnings` and workspace formatting passed.

## Boundary

These four batches were selected to expose automatic label gaps and shared-site fragments. They are development examples, not a representative holdout. An exploratory 20-molecule benzene batch still has **0/20 exact labelled rows** because the engines choose different alignments of symmetric ring sites; 18/20 can be made exact with one global R-label permutation, while two rows use a different equivalent embedding. The current change therefore does not establish general R-group accuracy superiority. Reaction and standardization comparisons, and chemistry, stereochemistry, search, and fingerprint quality remain separate gates.
