# Bounded NumPy graph matrices

`chematic.rdkit_compat` now exposes `GetAdjacencyMatrix` and
`GetDistanceMatrix`, using the existing native RDKit-compatible matrix
implementations. These functions were absent from the facade. They carry the
native size boundary into RDKit-style NumPy workflows associated with open
[RDKit #9647](https://github.com/rdkit/rdkit/issues/9647).
Chematic's existing Rust implementations already have checked dimensions and
fallible allocation; this adds the Python facade, rather than fixing the same
C++ overflow in chematic. RDKit is not a runtime dependency.

```python
from chematic import rdkit_compat as Chem

mol = Chem.MolFromSmiles("C=CO")
assert mol is not None
adjacency = Chem.GetAdjacencyMatrix(mol, useBO=True)
distance = Chem.GetDistanceMatrix(mol, useBO=True, useAtomWts=True)
assert adjacency.tolist() == [[0, 2, 0], [2, 0, 1], [0, 1, 0]]
assert distance.tolist() == [[1, .5, 1.5], [.5, 1, 1], [1.5, 1, .75]]
assert adjacency.dtype.name == distance.dtype.name == "float64"
```

The [RDKit API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html)
argument names and supported numeric definitions are retained:

| Function | Result |
| --- | --- |
| `GetAdjacencyMatrix(mol, useBO=False, emptyVal=0, force=False, prefix="")` | Owned C-contiguous NumPy int32 binary adjacency; `useBO=True` gives float64 endpoint bond contributions. Aromatic entries are 1.5; weighted dative entries are asymmetric. |
| `GetDistanceMatrix(mol, useBO=False, useAtomWts=False, force=False, prefix="")` | Owned C-contiguous float64 topological distances, with `1e8` between fragments. Weighted edges use reciprocal bond order, aromatic edges 2/3; atom-weighted diagonal is `6 / atomic_number`. |

Empty molecules produce shape `(0, 0)` with the correct dtype. Every call
returns fresh mutable storage; changing a result cannot change the molecule,
other results or later calls. The existing native methods still return lists.

All flags must be Python bools. `emptyVal` accepts integer index values,
including NumPy integers, but excludes Python/NumPy bools. Nonzero fill values
raise `ValueError`: RDKit's nonzero byte-fill behavior is not reproduced.
Both values of `force` compute a fresh result. There is no matrix cache,
so only empty `prefix` is supported; nonempty strings raise
`NotImplementedError`, while non-string prefixes raise `TypeError`.
Molecule properties are unchanged, including no cached matrix properties.

## Size and atom-index boundaries

Both functions refuse more than 1,024 **input** atoms before calling the native
matrix operation or allocating NumPy dense output. Explicit H atoms count
toward this limit even if preparation could remove them. The native adjacency
profile also caps bonds at 16,384 and checks prepared dense dimensions; the
distance profile retains its checked dimensions and fallible native buffers.
No limit override is added. This guards dense quadratic storage and the native
Floyd-Warshall operation, not just integer overflow.

The argument must be a `chematic.rdkit_compat.Mol`; native Mol objects and
RWMol editors must be wrapped or converted to a snapshot first. Inputs must
satisfy the existing [native graph-query](rdkit-graph-queries.md) and
[distance profile](rdkit-distance-matrix-limits.md). There is no new query,
unsanitized or general CTAB profile.

The native prepared model can remove explicit H atoms, whereas the facade
still exposes stored graph atom indices. If preparation changes the row count,
these functions raise `ValueError` rather than returning a matrix whose indices
disagree with `mol.GetNumAtoms()`. The named native methods remain available
for callers deliberately working with prepared-model indices. Explicit H
atoms added by native `rdkit_add_hydrogens()` are retained and supported.
General hydrogen-normalization changes to facade parsers are outside scope.

The facade does not add 3D/conformer distance storage, a shortest-path wrapper,
cache prefixes or matrix caching. It preserves the existing model's numeric
edge cases; it is not a new finite-coordinate validation API.

## Validation

The new tests cover manual matrices, aromatic/dative weights, empty graphs,
dtypes and ownership, independent repeated/forced results, source properties,
permutations and editor snapshots, atom-index mismatch refusal and strict
option types. Actual graphs at 1,025, 46,341, 65,535, 65,536 and 65,537 input
atoms are refused, with NumPy dense allocation functions instrumented to fail
if reached. A 1,024-atom chain succeeds for adjacency. Existing native graph,
matrix-limit and compatibility regressions also pass.

The independent RDKit 2026.09.1 oracle covers 32 original minimal structures,
12 added-H profiles and 110 atom-order variants. Across 1,320 matrix calls
and 51,552 compared cells, adjacency values/dtypes match exactly and the maximum
observed distance difference is zero. Both force settings, all supported weight
combinations, empty graphs, charged fragments and isotope/stereo examples are
included. This bounded corpus is not an all-input compatibility claim; large
unsafe upstream matrix calls were not executed.

Reproduction, source/import-fixture hashes and logs are under
`validation/results/compat-graph-matrices-2026-10-18/`. Tests use a verified
existing native extension, not a fresh current-main native build or published
wheel. No performance comparison or complete #9647 resolution is claimed.

## Python reference

::: chematic.rdkit_compat.GetAdjacencyMatrix

::: chematic.rdkit_compat.GetDistanceMatrix
