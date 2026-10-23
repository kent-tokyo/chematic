# Maintenance scripts

This directory contains current repository automation. Keep a script here only
when CI calls it, it regenerates a maintained fixture, it reproduces a current
evidence record, or it is part of the documented release process. One-off
investigations should be removed after their conclusion is recorded.

## Where scripts belong

| Purpose | Directory |
|---|---|
| Benchmarks | [`benchmarks/`](benchmarks/) |
| Validation and compatibility gates | [`validation/`](validation/) |
| Reproducible diagnostics | [`diagnostics/`](diagnostics/) |
| Fixture and report generation | [`generators/`](generators/) |
| Release automation | [`release/`](release/) |
| Shared Python helpers | [`lib/`](lib/) |
| Script tests | [`tests/`](tests/) |
| General repository checks | [`checks/`](checks/) |
| Evidence-record checks | [`checks/evidence/`](checks/evidence/) |

The root remains a reviewed legacy list, not the default destination for new
files. `tools/script_root_baseline.list` records its allowed direct children,
and `tools/check_script_layout.py` enforces that list in CI.

## Finding the right command

- Search `.github/workflows/` for CI and release entry points.
- Use [`scripts/check.sh`](check.sh) for the local repository gate bundle.
- Start from [`docs/reference/validation.md`](../docs/reference/validation.md) for current chemistry
  evidence and [`docs/reference/benchmark.md`](../docs/reference/benchmark.md) for measurements.
- Follow the dated record in `benchmarks/` or `validation/` before rerunning an
  evidence script; it defines the corpus, versions, hashes, and claim boundary.
- Use `--help` for manually run validation and diagnostic commands.

Frequently used entry points include:

| Task | Entry point |
|---|---|
| Release and version checks | `checks/check_release_*.py`, `checks/check_publish_graph.py`, `release/bump_version.py` |
| Documentation | `mkdocs build --strict`; the CI workflow is the canonical command |
| Benchmark index | `checks/check_benchmark_index.py` |
| Compatibility dashboard | `checks/check_compatibility_profiles.py`, `generators/generate_compatibility_dashboard.py` |
| Evidence records | `checks/evidence/check_*_evidence.py` |
| Cross-binding fingerprints | `validation/check_fingerprint_cross_binding.py` |
| Same-process format contracts | `validation/check_same_process_contracts.py` |
| Parser security | `diagnostics/run_isolated_parser_security.py` |
| Refactoring budget | `checks/check_refactor_budget.py` |
| Frozen trust holdouts | `checks/check_trust_holdouts.py`, `validation/score_trust_holdouts.py` |

## Retention rules

- Update CI, docs, imports, and `__file__`-based root lookup when moving a
  retained script.
- Put generated results under `benchmarks/` or `validation/`, never beside the
  script.
- Keep raw failures and row accounting visible in evidence output.
- Consolidate repeated mechanics into a shared helper or a named profile.
- Do not rewrite historical measurements when reorganizing current code.
- Delete a script with no current CI, documentation, issue, benchmark, fixture,
  or import reference after checking that Git history preserves it.

Historical public-artifact runners and one-off audits are kept in their release
tags rather than the current tree. Their measured records remain versioned
evidence and do not qualify the current release.
