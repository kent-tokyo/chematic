# Endpoint-aware alkene path stereo

[RDKit #9719](https://github.com/rdkit/rdkit/issues/9719) reports a cis/trans
classification that depends on the direction of a 1–4 traversal. chematic's
legacy bounds predicate has the same order-sensitive expression: the two
stored stereo-reference atoms are associated with the central bond's begin
and end, so an arbitrary reversed outer pair cannot be compared in that frame.

The additive Rust/Python query accepts the **complete four-atom path** and
normalizes its direction to the central bond before using the existing
predicate. Reversing all four indices preserves the result.

```python
import chematic

mol = chematic.from_smiles("F/C(Cl)=C(Br)/I")
assert mol.rdkit_alkene_path_stereo([0, 1, 3, 5]) == "trans"
assert mol.rdkit_alkene_path_stereo([5, 3, 1, 0]) == "trans"
assert mol.rdkit_alkene_path_stereo([0, 1, 3, 4]) == "cis"
```

Rust exposes
`chematic_3d::rdkit_embed::rdkit_alkene_path_stereo(&mol, [a, b, c, d])`
returning `Result<Option<AlkenePathStereo>, RdkitEmbedError>`, with
`AlkenePathStereo::{Cis, Trans}`. Python returns `"cis"`, `"trans"` or `None`.
These describe the relative side of the **selected outer atoms**, not the
whole bond's CIP E/Z label.

## Input and errors

Indices use the input/prepared molecular view's atom order. All four atoms
must be distinct, in range and consecutively bonded, with a double bond at
the center. An implicit H cannot be named by an atom index: call
`add_hydrogens()` and use its real index first.

Undeclared or unspecified central stereo returns `None`. The query uses only
declared stereo; it infers nothing from 2D/3D coordinates and does not invent
a relation for an unassigned alkene. Invalid integer paths raise `ValueError`
in Python, including negative indices. Wrong lengths/types are rejected by
the binding; integers outside its signed-64-bit extraction range raise
`OverflowError`. Rust returns an error for invalid paths. No graph, stereo
metadata or coordinates are mutated.

## Scope relative to the upstream report

The native legacy predicate's order dependence was reproduced directly on
the original small halogen control: outer pair `(0, 5)` returned different
values when passed as `(5, 0)`. The saved baseline probe and failure log are
under `validation/rdkit_issues/stereo-bounds/` and
`validation/results/stereo-bounds-2026-10-16/`.

This does **not** demonstrate a default chematic bounds-matrix defect. The
current builder already supplies oriented paths in ordinary cases; the
upstream screenshot was not reconstructed, and exploratory molecule controls
did not show a matrix change. The existing `rdkit_bounds_matrix`, its
RDKit-2026.09.1 compatibility predicate and all embedding defaults therefore
retain their behavior. This implements a checked path-query capability
motivated by #9719, not a claim to have fixed the complete upstream builder
or improved conformer quality.

## Validation

Original hand-derived fixtures cover both configurations, reverse SMILES
traversal, all four substituent combinations, explicit hydrogen and undeclared
stereo. Independent RDKit 2025.09.3 reference atoms are combined with an
unordered-set side-parity calculation, a different formulation from the
native endpoint normalization. Every atom permutation of each small fixture
is checked in both path directions: **4,344 permutations and 34,608 decisions,
zero mismatches**. Compressed per-permutation results and both decoded and
compressed hashes are recorded alongside the summary.

Local validation passed seven new Rust tests, ten existing bounds tests and
the pinned bounds-boundary suite, plus 56 Python tests (25 new and 31 existing
renumbering/internal-coordinate tests). These are source-candidate checks,
not released-package evidence, universal stereo parity or a speed benchmark.
No upstream source, screenshot reconstruction or attachments were copied.

```sh
PYTHONPATH=/path/to/candidate-package python \
  validation/rdkit_issues/stereo-bounds/validate.py \
  --output-dir /tmp/alkene-path-validation
```
