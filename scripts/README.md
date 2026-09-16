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
