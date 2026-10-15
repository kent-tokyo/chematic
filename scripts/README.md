# Maintenance scripts

Scripts are retained when they are used by CI, reproduce a checked-in evidence
record, generate a committed fixture, or provide a documented release operation.
One-off diagnostics should be removed after their result is captured unless a
current issue or benchmark links to them.

## Main entry points

| Purpose | Entry point |
|---|---|
| Open Babel SDF batch throughput | `build_openbabel_sdf_speed_corpus.py --offset` generates fixed unfiltered V2000 slices and refuses format changes; `bench_openbabel_sdf_throughput.py` validates parse and rewrite semantics separately and records 21 native timing/RSS blocks. See [the original packet](../benchmarks/2026-10-14-openbabel-sdf-speed.md) and [the reader follow-up](../benchmarks/2026-10-14-openbabel-sdf-read-followup.md). |
| Open Babel V3000 block throughput | `build_openbabel_v3000_speed_corpus.py` retains unfiltered force-V3000 slices; `bench_openbabel_v3000_throughput.py` compares pre-split native MOL API batches in 21 rotating blocks. Rejected outputs keep their row positions and fail the semantic gate. The [packet](../benchmarks/2026-10-14-openbabel-v3000-speed.md) separates diverse chemistry from controlled chain-size scaling. |
| Python descriptor/fingerprint speed | `bench_python_speed_frontier.py` compares exact source wheels, RDKit and COSMolKit across first-use, repeated-use and actual parse-inclusive loops. `verify_python_descriptor_frontier.py` and `verify_python_morgan_frontier.py` retain float bits, errors and fingerprints. See the [core-model preparation follow-up](../benchmarks/2026-10-14-python-model-view-speed.md) and [sanitized-view reuse](../benchmarks/2026-10-14-python-sanitized-view-reuse.md). |
| Source Python speed frontier | `bench_python_speed_frontier.py` compares baseline/candidate wheels, RDKit and COSMolKit with first-use, repeated-use and parse-inclusive boundaries; `verify_python_morgan_frontier.py` dumps the fixed 14,072-input Morgan parity corpus. See [the dated packet](../benchmarks/2026-10-14-python-speed-frontier.md). |
| Scalar/stereo speed follow-up | `bench_python_speed_frontier.py --offset` selects a fixed corpus slice; `verify_python_descriptor_frontier.py` preserves first/repeated float bits, errors and nonfinite values for four descriptors. See [the scalar/stereo packet](../benchmarks/2026-10-14-python-scalar-stereo.md) and [the first-use perception follow-up](../benchmarks/2026-10-14-python-first-use-followup.md). |
| Long scalar reuse batches | `bench_python_scalar_reuse.py` retains and validates every result of one million prepared TPSA/ASA calls per arm, independently of the short-batch frontier protocol. Both original and additional protocols are in [the follow-up packet](../benchmarks/2026-10-14-python-scalar-stereo.md). |
| Release/version checks | `check_release_*.py`, `check_publish_graph.py`, `bump_version.py` |
| Compatibility dashboard | `check_compatibility_profiles.py`, `generate_compatibility_dashboard.py` |
| Benchmark index | `check_benchmark_index.py` |
| Refactoring debt ratchet | `check_refactor_budget.py` checks `validation/refactor_budget_v1.json`; `--print-baseline` prints a reviewed replacement without writing files |
| Open Babel file-I/O comparison | `check_openbabel_file_io_contract.py`, `check_openbabel_file_io_gate.py`, `check_openbabel_file_io_semantics.py`, `bench_openbabel_file_io_paired.py`, `bench_openbabel_file_io_same_process.py`, and `check_openbabel_file_io_evidence.py` |
| Parser security | `run_isolated_parser_security.py` |
| Browser comparison | `bench_browser_wasm_vs_rdkit_isolated.py` |
| Published browser replication | `check_published_browser_paired.py`, `check_v1030_published_browser_morgan_rows.py`, and `check_published_linux_browser_replication.py` validate the pinned Chromium and Ubuntu three-engine records |
| Sealed cohort | `prepare_sealed_accuracy_cohort.py` |
| V3000 external readers | `v3000_*_gate.py` |
| Stereo gates | `stereo_torture_suite_gate.py`, `stereo_spelling_invariance_gate.py` |
| Open Babel 0D MOL stereo source check | `compare_openbabel_0d_stereo.py` compares before and after CLI binaries on identical Open Babel-generated records; accepts JSONL or `.smi` corpora and records generation failures and source-gate exclusions |
| A6 source-wheel termination diagnostic | `public_package_3d_chematic.py` runs the fixed 265 rows; `summarize_a6_mmff94_termination.py` accounts for every outcome; `compare_a6_source_published_rows.py` checks same-host success retention without claiming the full quality gate |
| A6 platform gate (#739) | `pipeline_v2_vs_rdkit_common_scorer` scores one platform's 265 rows; `a6_platform_gate.py` requires 265 successes, all independently sound, stereo-clean and clash-free, and reports whether the coordinates equal another platform's row by row |
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
