# Maintenance scripts

Scripts are retained when they are used by CI, reproduce a checked-in evidence
record, generate a committed fixture, or provide a documented release operation.
One-off diagnostics should be removed after their result is captured unless a
current issue or benchmark links to them.

## Main entry points

| Purpose | Entry point |
|---|---|
| Release/version checks | `check_release_*.py`, `check_publish_graph.py`, `bump_version.py` |
| Compatibility dashboard | `check_compatibility_profiles.py`, `generate_compatibility_dashboard.py` |
| Benchmark index | `check_benchmark_index.py` |
| Parser security | `run_isolated_parser_security.py` |
| Browser comparison | `bench_browser_wasm_vs_rdkit_isolated.py` |
| Sealed cohort | `prepare_sealed_accuracy_cohort.py` |
| V3000 external readers | `v3000_*_gate.py` |
| Stereo gates | `stereo_torture_suite_gate.py`, `stereo_spelling_invariance_gate.py` |
| Cross-binding Node dump | `binding_dump.mjs <mode>` |
| Published v1.0.30 artifact packet | `check_v1030_artifact_packet.py`, `emit_rdkit_smarts_oracle.py`, `run_published_npm_chemistry_lane.mjs`, `run_published_npm_63op_outputs.mjs`, `run_reaction_compatibility_v2.py`, and the `tools/published_rust_gate` lockfile |
| Published-wheel output and speed comparison | `bench_python_op_matrix_vs_rdkit.py --outputs-only/--paired`, `check_published_python_version_outputs.py`, `bench_published_python_versions.py` |

## Organization rule

- Shared mechanics belong in one helper; callers select a named profile/mode.
- Evidence scripts write versioned JSON and keep raw failures visible.
- Historical measurements keep their original script contract in Git history.
- Historical generated-package reports are immutable evidence, not a check
  against the current source export surface.
- A script with no CI, documentation, issue, benchmark, or fixture reference is
  a deletion candidate.
- Generated packages, build output, and temporary corpora do not belong in this
  directory.
