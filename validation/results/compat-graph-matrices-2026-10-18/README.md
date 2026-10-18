# Bounded NumPy compatibility matrices — 2026-10-18

Related open issue: https://github.com/rdkit/rdkit/issues/9647

The native Rust matrix implementations already have checked dimensions and
allocation/resource limits. The missing capability was the RDKit-style facade
entry points. This adds typed NumPy `GetAdjacencyMatrix`/`GetDistanceMatrix`,
without claiming the same C++ memory-safety defect existed in chematic or
resolving every affected upstream API.

- Base: `f99c7b8ab33098a2e1fd9fd568aee370bf79ab38`.
- `pytest.log`: 245 tests pass, including 60 new cases. Actual graphs at
  1025/46341/65535/65536/65537 atoms are refused, with NumPy dense allocation
  functions instrumented to fail if reached. Adjacency succeeds at 1024 atoms.
  Existing graph-query, native matrix-limit and compatibility tests pass.
- `oracle.py`, `oracle-summary.json`, `oracle.log`: independent RDKit 2026.09.1
  calls, 32 original minimal structures, 12 explicit retained-H profiles and
  110 atom-order variants. All 1320 matrix calls / 51552 cells agree:
  adjacency values and dtype exactly; maximum observed distance error zero.
  Both force values, both bond-weight options and both atom-weight options
  are checked. No unsafe large upstream matrix call is executed.
- `typing-valid.log`, `typing-invalid.log`: strict helper/public consumer pass;
  negative consumer exits 1 with four expected public argument errors. Legacy
  imported modules/stubs use silent imports; package-wide strict typing is not
  claimed.
- `docs-and-typing.log`, `mkdocs.log`: Python 3.9 source grammar, executed guide,
  strict build, rendered guide and two public API anchors pass.
- `source-gates.log`: source stubs, cross-binding manifest, binding inventory,
  unsafe surface, refactor budget and naming gates pass. No wheel was provided.
- `native-source-comparison.log`: matrix source modules, graph-query bindings
  and the native atom_count/distance methods match the verified old extension
  source. This is not whole-extension equivalence or a current-main rebuild.
- `provenance.json`: source, runtime fixture, reference and evidence hashes.

Profile: at most 1024 raw input atoms, including H nodes; native adjacency
also limits bonds to 16384. The prepared row count must agree with the facade
Mol's stored atom count; H-removing inputs refuse rather than return wrong
indices. Native prepared-index APIs remain unchanged. Both force values
recompute; nonempty cache prefixes and nonzero adjacency fills explicitly
refuse. Flags require Python bool; zero fill supports integer index values
excluding Python/NumPy bool. Every NumPy output owns independent storage.

No 3D/conformer facade, cache storage, shortest-path wrapper, new native
algorithm, published-wheel qualification, performance comparison or all-input
RDKit parity claim is included.

## Reproduce

Add a copy of this revision's Python package with a compatible native extension
to PYTHONPATH. For the oracle also add an independent RDKit 2026.09.1 install.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_compat_graph_matrices.py \
  crates/chematic-py/tests/test_rdkit_graph_queries.py \
  crates/chematic-py/tests/test_rdkit_matrix_limits.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 python3.13 validation/results/compat-graph-matrices-2026-10-18/oracle.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-matrices-site
```

Strict mypy uses `follow_imports=silent`, `follow_imports_for_stubs=True`,
`ignore_missing_imports=True`, and MYPYPATH set to the candidate package.
Check `_graph_matrices.py` with `typing_valid.py`, then `typing_invalid.py`
separately and expect four errors. Unit tests require neither RDKit nor any
upstream PoC file. All oracle structures are authored minimal SMILES, not
copied molecular datasets or issue attachments.
