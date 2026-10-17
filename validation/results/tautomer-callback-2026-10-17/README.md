# Tautomer callback evidence — 2026-10-17

RDKit [#9601](https://github.com/rdkit/rdkit/issues/9601) is open and includes
`TautomerEnumerator.SetCallback` among its remaining typed boundaries. This
change adds a chematic-native typed callable API and actual in-loop cancellation,
not RDKit's class interface. Existing Rust transforms, stereo normalization,
iteration order and default enumeration entry points are retained.

Results:

- `pytest.log`: 224 passed (62 added callback cases, 162 existing compatibility).
- `rust-tests.log`: full chematic-smiles library, 403 passed, 1 existing ignored.
- `rdkit-oracle.json`: RDKit 2026.09.1, 24 molecules × 12 policies = 288/288
  matched; final SMILES, final map keys, final status and checkpoint-key traces.
  557 candidate checkpoints. Transform-rule counts are candidate diagnostics;
  RDKit does not expose that count in its Python callback result.
- `typing-valid.log`: strict check of the new helper and positive public consumer.
- `typing-invalid.log`: exactly five expected errors (callback return, molecule,
  limit, result assignment, frozen result). This is not a full-package strict
  mypy claim.
- `native-build.log`: fresh local dev extension, including native InChI, built
  with Python 3.13.6 and macOS dynamic symbol lookup. This is not a wheel build
  or hosted-CI result. Native checksum/tool versions/source hashes are recorded
  in `provenance.json`; no earlier native fixture is used.
- `mkdocs.log`: strict documentation build passed.
- Source gates: stubs, cross-binding manifest, binding surface, unsafe surface,
  refactor budget and product names passed. The resolved `enumerate` complexity
  allowance was removed; new functions meet the existing soft limits.

Commands, from the candidate worktree:

```sh
CARGO_TARGET_DIR=/private/tmp/chematic-sparse-tanimoto-target CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 CARGO_INCREMENTAL=0 cargo test -p chematic-smiles --lib
PYO3_PYTHON=/opt/homebrew/bin/python3.13 CARGO_TARGET_DIR=/private/tmp/chematic-sparse-tanimoto-target CARGO_PROFILE_DEV_DEBUG=0 CARGO_INCREMENTAL=0 cargo rustc -p chematic-py --features native-inchi -- -C link-arg=-undefined -C link-arg=dynamic_lookup
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-tautomer-callback-python python3.13 -m pytest crates/chematic-py/tests/test_tautomer_control.py crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/chematic-tautomer-callback-python:/private/tmp/chematic-sparse-tanimoto-rdkit python3.13 validation/validate_tautomer_callback.py --output validation/results/tautomer-callback-2026-10-17/rdkit-oracle.json
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-tautomer-callback-site
```

The candidate Python package directory contains the candidate `.py`/`.pyi` files
and the freshly built `libchematic.dylib` linked as `chematic.so`. The oracle
path contains the installed RDKit 2026.09.1 wheel. Both are local tooling, not
runtime dependencies added to chematic.

The strict mypy configuration uses `strict=True`, `follow_imports=silent`,
`follow_imports_for_stubs=True` and `ignore_missing_imports=True`; the positive
check explicitly includes `_tautomer_control.py` and
`validation/typing/tautomer_control_valid.py`. The negative check uses
`MYPYPATH=crates/chematic-py/python`, `--follow-imports=silent` and
`--ignore-missing-imports`. All checks use mypy 1.8.0 with `--no-incremental`.

Scope: concrete molecule profile only, positive Python limits ≤10000, graph
atoms ≤1024, synchronous bool-returning callbacks. Cancellation is cooperative
at transform checkpoints and still normalizes partial-result stereo. It is not
preemption during matching. Oracle agreement on these cases does not establish
universal chemistry parity or support every RDKit option.
