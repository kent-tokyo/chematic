# Parent API documentation and typing — 2026-10-17

RDKit [#9651](https://github.com/rdkit/rdkit/issues/9651) reports missing generated
Python documentation for the XyzParent convenience functions. The live issue
is open. chematic has a directly measurable analogous gap: seven native Parent
methods exist but are omitted from its Python type stub and generated reference.

This change publishes those existing methods' signatures, defaults and behavior;
it does not add or change a chemical transformation. `ParentStage` and
`SuperParentReport` provide runtime-exported TypedDict definitions for the existing
ordinary report dictionaries. The guide distinguishes native chemistry from
related RDKit names and explains fragment/charge/isotope/stereo scope, explicit
H retention, composition order, unsigned budgets, partial/timeout status and
property copying.

## Before and after

- `baseline-source-gate.json` / `.log`: the new gate fails against exact baseline
  files extracted from the base commit using `git show HEAD:path` before this
  change was committed. Native methods exist, but their typed declarations and
  report type exports/reference entries/guide are missing.
- `typing-baseline.log`: a valid seven-method consumer fails with 9 errors
  against that baseline stub. `typing-valid.log`: the same consumer and the new
  report-types module pass a scoped strict mypy 1.8.0 check after the fix.
- `typing-invalid.log`: 6 expected diagnostics cover wrong limit types, incorrect
  Mol/status assignments, invalid stage literal and a nonexistent dictionary key.
  These are scoped checks, not a full-package strict mypy claim.
- `baseline-html-gate.json` / `.log`: the previously built baseline documentation
  lacks all seven method anchors, both report-type anchors and the Parent guide.
  The old site came from the preceding crystal-OBJ worktree based on the same
  main commit; its Parent surface is unchanged from main.
- `native-default-negative-control.json`: changing the native PyO3 default from
  16 to 17 in an isolated copy is rejected as a signature mismatch.
- `api-gate.json`: the final generated site contains all seven method and both
  report-type anchors plus the guide. The source gate also verifies native
  existence, return annotations, budget defaults and required documentation.
- `pytest.log`: 10 passed, 145 deselected. The new executable-guide test runs all
  three Python blocks (four simple methods, composed stage output, limits,
  timeout and empty graph); nine existing Parent contract tests also pass.
- `mkdocs.log`: strict build passed with the repository's pinned documentation
  dependencies. Source gates for stubs, binding surface, cross-binding manifest,
  unsafe surface, refactor budget, naming and diff whitespace pass.
- Workflow YAML parses and the new job/gate steps are present. Hosted CI is
  separate from these local results.

## Recurrence guard

`scripts/check_parent_api_docs.py` runs without importing the native package.
The existing Python 3.9/3.13 stub jobs now run its source checks. A new lightweight
CI job builds MkDocs with pinned dependencies and checks the actual HTML. The
existing Pages workflow checks the rendered result before deployment and now
also responds to Python package/doc-gate changes. No deployment was invoked.
The guide's executable examples run with the normal Python test suite.

## Native evidence boundary

No Rust implementation changed. Runtime examples/contracts use the existing
native fixture from source commit `71d9c12b0c16a8e6b23e4078479d5754230c831b`,
version 1.1.6, SHA-256
`c6db880b6f5ca4ea0388f238c24e8a5b55751daab1084cf0d1d7b56e2427257b`.
The candidate Python package contains the candidate Python files and that
extension. This is not a fresh native build, wheel, or release validation.

`git diff --quiet 71d9c12b0c16a8e6b23e4078479d5754230c831b HEAD -- crates/chematic-chem/src/parent.rs crates/chematic-chem/src/standardize.rs crates/chematic-chem/src/tautomer.rs crates/chematic-chem/src/tautomer crates/chematic-py/src/mol_methods/transformations.rs`
succeeds at the candidate base. These Parent implementation/binding paths are
unchanged from the fixture source. `provenance.json` records source/base/tool
versions and hashes. No RDKit chemistry-parity claim is made or needed for this
documentation/typing correction.

## Reproduction

```sh
python3.13 scripts/check_parent_api_docs.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-parent-api-docs-site
python3.13 scripts/check_parent_api_docs.py --site-dir /private/tmp/chematic-parent-api-docs-site --output validation/results/parent-api-docs-2026-10-17/api-gate.json
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-parent-api-python python3.13 -m pytest crates/chematic-py/tests/test_parent_api_docs.py crates/chematic-py/tests/test_cross_binding_contract.py -k parent -q -p no:cacheprovider
```

Strict mypy configuration: `strict=True`, `follow_imports=silent`,
`follow_imports_for_stubs=True`, `ignore_missing_imports=True`; explicitly check
`_parent_types.py` and `validation/typing/parent_api_valid.py`, with
`MYPYPATH=crates/chematic-py/python` and `--no-incremental`. For the negative
consumer use `--follow-imports=silent --ignore-missing-imports --no-incremental`
and `validation/typing/parent_api_invalid.py`; expect six errors and exit 1.
