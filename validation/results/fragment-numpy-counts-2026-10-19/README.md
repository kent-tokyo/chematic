# Typed NumPy cut-count destinations

Base: `00596fb1ea5397524f98e3789c41f17ee8616fdc`. Checked 2026-10-19.
RDKit [#9601](https://github.com/rdkit/rdkit/issues/9601) remains open and
lists `FragmentOnBonds.cutsPerAtom` among mutable-output typing tasks.
The existing bounded chematic API accepted lists only. This addition accepts
plain, writable, non-overlapping 1D uint32 arrays with exactly one item per atom.
This is a chematic extension inspired by that boundary, not a nanobind fix,
new fragmentation algorithm, RDKit NumPy-output parity or whole-issue resolution.

## Results

- Existing and new Python fragmentation tests: **80 passed**, including 25 new
  cases. Four valid layout cases fail against the prior list-only helper.
- Live RDKit **2026.09.1**: **260/260** nonempty cut/layout combinations over
  ten supplied molecules agree on every atom's count and canonical fragment
  graph. The reference uses RDKit list outputs; NumPy outputs are candidate-only.
  The four layouts are plain, strided, reversed and big-endian uint32.
- Scoped strict mypy passes the changed helper and valid public consumer.
  The negative consumer reports exactly two expected argument errors.
  Imported legacy declarations are followed silently; package-wide strict
  typing is not claimed.
- Python source/stub grammar, source binding/unsafe/refactor/naming gates,
  Black, executed guide examples, strict MkDocs and whitespace checks pass.

Unit tests cover dtype/shape/read-only/overlap/subclass rejection, no native
call on invalid output, unchanged array on native refusal and injected staging
allocation failure, isolated result counts, empty cuts and empty molecules.
Existing list semantics and bridge-single-bond restrictions are unchanged.
No concurrent-caller synchronization or input-callback rollback is promised.

## Runtime and reproducibility

Local source Python is loaded from
`/private/tmp/chematic-fragment-numpy-20261019`, with the existing v1.1.7 native
extension from `/private/tmp/chematic-blog-v117-venv`. Its SHA256 is recorded
in `oracle-summary.json`. Fragmentation Rust core and binding sources have no
diff from tag `v1.1.7` to the base. This is not a fresh current-main native build
or a published wheel test. No Rust code changed. Unit tests use NumPy 2.5.3;
the RDKit environment supplies NumPy 2.5.4 for the reference comparison.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
PYTHONPATH=/private/tmp/chematic-fragment-numpy-20261019 \
python3.13 -m pytest crates/chematic-py/tests/test_bond_fragment.py \
  crates/chematic-py/tests/test_bond_fragment_numpy.py -q -p no:cacheprovider

PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH=/private/tmp/chematic-fragment-numpy-20261019:/private/tmp/chematic-rdkit-2026091/lib/python3.13/site-packages \
python3.13 validation/results/fragment-numpy-counts-2026-10-19/oracle.py

MYPYPATH=/private/tmp/chematic-fragment-numpy-20261019 \
python3.13 -m mypy \
  --config-file validation/results/fragment-numpy-counts-2026-10-19/mypy.ini \
  crates/chematic-py/python/chematic/_bond_fragment.py \
  validation/results/fragment-numpy-counts-2026-10-19/typing_valid.py
```
