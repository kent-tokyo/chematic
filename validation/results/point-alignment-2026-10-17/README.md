# Typed point alignment — 2026-10-17

This is a Python/NumPy implementation of the point-input item in RDKit #9601.
It adds `rdkit_compat.rdAlignment.GetAlignmentTransform`; it does not change
native Rust alignment or claim complete wrapper compatibility. Only the
default `maxIterations=50` is accepted by the SVD backend. Native binding
inventory is unchanged.

- Clean baseline Python source at `f20304f4` has no `rdAlignment` namespace.
- 224 selected Python tests pass, including 62 new alignment tests.
- 576 live RDKit 2026.09.1 comparisons pass for lists, arrays and actual RDKit
  Point3D objects; positive weights, default weights and reflections are
  included. Maximum absolute SSD difference: 2.3877788635218167e-10.
- 1,152 coordinate/weight-scaling checks pass. Returned rotations are checked
  for orthogonality and determinant sign, and SSD against transformed residuals.
- Maximum nondegenerate matrix-element difference from RDKit's quaternion
  solver: 3.641913534035979e-5. This is a tolerance diagnostic; degenerate
  transforms may be nonunique. No exact-transform parity claim is made.
- The helper passes isolated strict mypy. Valid public calls type-check;
  wrong point/weight/reflect/iteration types produce four expected arg-type
  errors. Package-wide strict mypy remains unsuccessful due to existing
  generic annotations in `__init__.pyi`; its complete diagnostic is retained.
  Those unrelated stubs are not changed or claimed as strict-clean.
- Refactor, source-stub, binding inventory, cross-binding manifest, unsafe
  surface, product naming and strict MkDocs checks pass. The two production
  Python files parse under Python 3.9 grammar (runtime validation used 3.13).

The native extension is reused from source commit
`e2041226b67030108acb00309da5f0659f636c1d` as import scaffolding. The new solver
calls NumPy only; this is not a fresh native build, wheel packaging or registry
publication validation. Source/fixture hashes and versions are in the
manifest and corpus report.

Reproduce with candidate Python source and a compatible chematic extension:

```sh
export PYTHONPATH=/private/tmp/chematic-point-alignment-python:/private/tmp/chematic-sparse-tanimoto-rdkit
export PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
python3.13 -m pytest crates/chematic-py/tests/test_point_alignment.py crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
python3.13 validation/validate_point_alignment.py
MYPYPATH=crates/chematic-py/python python3.13 -m mypy --follow-imports=silent --ignore-missing-imports --no-incremental validation/typing/point_alignment_valid.py
# Expected exit 1 with exactly four arg-type errors:
MYPYPATH=crates/chematic-py/python python3.13 -m mypy --follow-imports=silent --ignore-missing-imports --no-incremental validation/typing/point_alignment_invalid.py
# Isolate helper strictness from legacy package stubs:
mkdir -p /private/tmp/chematic-point-alignment-typecheck
cp crates/chematic-py/python/chematic/_point_alignment.py /private/tmp/chematic-point-alignment-typecheck/point_alignment.py
python3.13 -m mypy --strict --follow-imports=silent --ignore-missing-imports --no-incremental /private/tmp/chematic-point-alignment-typecheck/point_alignment.py
```
