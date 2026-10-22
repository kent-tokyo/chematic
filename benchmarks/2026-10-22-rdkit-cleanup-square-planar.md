# Cleanup square-planar metal stereo (2026-10-22)

Pinned comparison: RDKit 2026.09.1, CheMatic source based on `b7817c64`, with the source file hashes and corpus hashes recorded in the machine reports. The candidate remaps the square-planar tag when Cleanup disconnects all three bonds of a metal center with one vacant site. RDKit cycles `SP1 → SP2 → SP3 → SP1` for this case; a four-ligand center retains its tag. Direct fixtures cover all three permutations, Pt and Pd, and the four-ligand control.

| Corpus | Comparable rows | Previous Cleanup exact | Candidate Cleanup exact | Fixed | Regressed |
| --- | ---: | ---: | ---: | ---: | ---: |
| `stress-v1.smi` | 4,069 | 4,055 | 4,056 | 475 | 0 |
| Independent 10k | 10,000 | 9,997 | 9,997 | 0 | 0 |

Candidate Normalize and Reionize match RDKit on all comparable rows in both corpora. The stress comparison uses the later Reionize source report as its baseline, so the one fixed row is attributable to this square-planar patch. The current stress source has 13 Cleanup mismatches. Three of these have RDKit Cleanup SMILES that RDKit itself cannot parse again (indices 1497, 1630, 2677); CheMatic retains reparsable charged forms. Exact-string disagreement for these rows is not evidence of lower chemical validity.

Evidence: [stress machine report](../validation/results/rdkit-cleanup-square-planar-stress-2026-10-22.json) and [10k machine report](../validation/results/rdkit-cleanup-square-planar-10k-2026-10-22.json). Each report records the source dump hash, input hash, RDKit version, fixed and regressed indices, and non-reparsable outputs. The gate was run against the previous checked-in reports with `scripts/rdkit_cleanup_source_gate.py --previous-source`. `cargo test --locked -p chematic-chem`, `cargo clippy --locked -p chematic-chem --all-targets -- -D warnings`, and `cargo fmt --all -- --check` passed.

This change only addresses Cleanup. The reaction and R-group comparison boundaries in [the workflow packet](2026-10-18-rdkit-workflows-followup.md) remain open; this result does not establish toolkit-wide superiority to RDKit.
