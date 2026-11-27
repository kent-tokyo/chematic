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
- Use [`bash scripts/check.sh`](check.sh) for the local gate bundle. Use the workflow for its runtime matrix and required CI checks.
- Start from [`docs/reference/validation.md`](../docs/reference/validation.md) for current chemistry
  evidence and [`docs/reference/benchmark.md`](../docs/reference/benchmark.md) for measurements.
- Follow the dated record in `benchmarks/` or `validation/` before rerunning an
  evidence script; it defines the corpus, versions, hashes, and claim boundary.
- Use `--help` for manually run validation and diagnostic commands.

Run these entry points from the repository root.

| Task | Entry point |
|---|---|
| Documentation and release wording | `python3 scripts/checks/check_release_docs_consistency.py` and `python3 scripts/checks/check_product_names.py` |
| Documentation site | `mkdocs build --strict` in the [pinned environment](../CONTRIBUTING.md#documentation-build) |
| Benchmark and roadmap indexes | `python3 scripts/checks/check_benchmark_index.py` and `python3 scripts/checks/check_roadmap_disposition.py` |
| Refactoring budget | `python3 scripts/checks/check_refactor_budget.py` |
| Binding contracts | `python3 scripts/checks/check_cross_binding_manifest.py` and `python3 scripts/checks/check_binding_surface.py` |
| Compatibility dashboard | [Profile checker](checks/check_compatibility_profiles.py) and [dashboard generator](generators/generate_compatibility_dashboard.py) |
| Release and version checks | [Release checklist](release/README.md) |
| Recorded evidence and holdouts | [Evidence checks](checks/evidence/) and [validation runners](validation/) |

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
