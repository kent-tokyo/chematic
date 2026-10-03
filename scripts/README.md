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
| Published browser replication | `check_published_browser_paired.py`, `check_v1030_published_browser_morgan_rows.py`, and `check_published_linux_browser_replication.py` validate the pinned Chromium and Ubuntu three-engine records |
| Sealed cohort | `prepare_sealed_accuracy_cohort.py` |
| V3000 external readers | `v3000_*_gate.py` |
| Stereo gates | `stereo_torture_suite_gate.py`, `stereo_spelling_invariance_gate.py` |
| A6 source-wheel termination diagnostic | `public_package_3d_chematic.py` runs the fixed 265 rows; `summarize_a6_mmff94_termination.py` accounts for every outcome; `compare_a6_source_published_rows.py` checks same-host success retention without claiming the full quality gate |
| Cross-binding Node dump | `binding_dump.mjs <mode>` |
| Published v1.0.30 artifact packet | `check_v1030_artifact_packet.py`, `check_v1030_published_chemistry_residuals.py`, `emit_rdkit_smarts_oracle.py`, `run_published_npm_chemistry_lane.mjs`, `run_published_npm_63op_outputs.mjs`, `run_reaction_compatibility_v2.py`, `reaction_template_map_gate.py`, `check_published_rust_63op_outputs.py`, and the version-pinned `tools/published_rust_gate*` lockfiles |
| Published-wheel output and speed comparison | `bench_python_op_matrix_vs_rdkit.py --outputs-only/--paired`, `check_published_python_version_outputs.py`, `bench_published_python_versions.py` |
| Isolated Python speed and process RSS | `bench_published_python_isolated_paired.py` (20 fresh-process AB/BA blocks), `check_published_python_isolated_paired.py` (archived six-lane integrity gate) |
| Isolated published-Rust speed and process RSS | `bench_published_rust_isolated_paired.py` builds two `--locked` crates.io graphs and runs fresh-process AB/BA blocks; `check_published_rust_isolated_paired.py` verifies archived outputs, paired intervals and RSS |
| Isolated published npm/WASM speed and process RSS | `bench_published_wasm_paired.mjs` pins v1.0.29/v1.0.30/RDKit.js tarballs and runs 20 fresh-process ABBA/BAAB blocks; `check_published_wasm_paired.py` verifies the archived output gate, raw block accounting, intervals and whole-process RSS |

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
