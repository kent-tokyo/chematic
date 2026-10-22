# Symmetric R-group column alignment (2026-10-22)

This is a source-level labelled R-group correction based on `origin/main` at
`e5284b23c3d5bd7208b8e15816093e495086378b`. The comparator is RDKit
2026.09.1 on local macOS arm64, CPython 3.13.6. The imported RDKit package
came from a wheel with SHA-256
`29b6176ace255989abde23d40d2187b5ed5effdb2f00dadaced24f181ab27046`;
the loaded R-group module hash is in the machine report.

For unlabelled benzene and the ordered molecules `Cc1c(Cl)cccc1`,
`Cc1ccc(Cl)cc1`, `Cc1cc(Cl)ccc1`, the previous CheMatic source assigned
methyl to a constant column and moved chloro among columns. RDKit keeps
chloro in `R4` and moves methyl through `R3`, `R1`, and `R2`.

The candidate breaks the first-row symmetry tie by assigning the higher
atomic-number substituent root to the earlier automatic core site. Once a
column exists, the existing fragment agreement and fingerprint similarity
scores still choose later rows. Explicitly labelled cores keep their previous
tie-break policy.

| Fixed source panel | Previous | Candidate | RDKit 2026.09.1 |
| --- | ---: | ---: | ---: |
| 11 small three-molecule batches, exact core/R-group graphs and labels | 10/11 in a prior 1.2.0 wheel diagnostic; the newly added Rust test failed on starting source | **11/11** | reference |

The previous wheel diagnostic is separately identified by SHA-256
`aa635ac5c5ce1a18e5dcd4edd88dbd5fa1352168fb9d7078c8995f6e5f8a6a00`;
it is not used as proof of exact starting-source identity. The source-built
candidate result and every row are recorded in the [machine report](../validation/results/rdkit-rgroup-symmetry-source-2026-10-22.json), which hashes
the source files, example binary, and RDKit module. Reproduce it by building
`cargo build --locked -p chematic-chem --example rgroup_gate_probe` and running
`scripts/validation/rdkit_rgroup_source_gate.py --binary PATH --output PATH`
with RDKit 2026.09.1.

The new Rust reproducer failed before the correction and passes after it.
All 33 R-group unit tests and the full `chematic-chem` suite passed (1,032
unit tests passed, one existing ignored, plus integration and documentation
tests). Clippy with warnings denied and the refactor-budget gate passed.

This panel is small and exposed. The earlier 20-row symmetric-ring batch
remains a separate unresolved comparison. No broad R-group accuracy, speed,
Python-wheel, reaction, or Cleanup superiority over RDKit is established by
this result.
