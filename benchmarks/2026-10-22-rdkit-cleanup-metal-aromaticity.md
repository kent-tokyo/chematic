# Cleanup metal-ligand aromaticity after sanitization (2026-10-22)

Pinned comparison: RDKit 2026.09.1, CheMatic source based on `83bbd279`. The source, corpus, and dump hashes are in the machine reports. The fixed stress row (index 1686) has a carbon marked aromatic in the input but kekulized to non-aromatic by RDKit sanitization. Cleanup now uses perceived aromaticity when deciding whether to retain a carbon-metal bond. When that formerly aromatic carbon loses its metal bond, it takes the implicit hydrogen written by RDKit. A benzene-bound Pt negative control retains its bond.

| Corpus | Comparable rows | Previous Cleanup exact | Candidate Cleanup exact | Fixed | Regressed |
| --- | ---: | ---: | ---: | ---: | ---: |
| `stress-v1.smi` | 4,069 | 4,056 | 4,057 | 1 | 0 |
| Independent 10k | 10,000 | 9,997 | 9,997 | 0 | 0 |

Normalize and Reionize remain exact on all comparable rows in both corpora. The stress Cleanup result still differs on 12 rows: 3 acid/base outputs that RDKit itself cannot reparse, 1 canonical representation difference, 1 metal charge assignment, and 7 metal disconnection cases. Agreement with RDKit is a compatibility measure, not independent chemical truth.

Evidence: [stress machine report](../validation/results/rdkit-cleanup-metal-aromaticity-stress-2026-10-22.json) and [10k machine report](../validation/results/rdkit-cleanup-metal-aromaticity-10k-2026-10-22.json). Generate each source dump with `cargo run --locked --offline -q -p chematic-chem --example cleanup_rdkit_dump < CORPUS > DUMP`, then run `scripts/rdkit_cleanup_source_gate.py --dump DUMP --output REPORT --previous-source PREVIOUS_REPORT` using the pinned RDKit 2026.09.1 environment. Pass `--corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi` for the 10k comparison.

`python3 scripts/check_standardization_source_evidence.py`, `cargo test --locked -p chematic-chem`, `cargo clippy --locked -p chematic-chem --all-targets -- -D warnings`, and `cargo fmt --all -- --check` passed locally.

The reaction and R-group comparison boundaries remain as described in [the workflow packet](2026-10-18-rdkit-workflows-followup.md). This result does not establish overall superiority to RDKit.
