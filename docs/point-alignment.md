# Typed point alignment

`chematic.rdkit_compat.rdAlignment.GetAlignmentTransform` aligns paired 3D
points and returns their weighted sum of squared deviations (SSD) and a 4×4
probe-to-reference transform. It addresses the typed point-input item in
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601), verified open on
2026-10-17. This is a new Python compatibility API; it uses NumPy SVD and does
not replace the native Rust `align_coords` API.

```python
from chematic.rdkit_compat import rdAlignment

ssd, matrix = rdAlignment.GetAlignmentTransform(
    refPoints=[[0., 0., 0.], [1., 0., 0.]],
    probePoints=[[3., 2., 1.], [3., 3., 1.]],
    weights=[1., 2.],
)
# For N×3 NumPy row coordinates:
# aligned_probe = probe @ matrix[:3, :3].T + matrix[:3, 3]
```

The API accepts a sequence of coordinate triples, an ordinary numeric NumPy
N×3 array, or a sequence of point objects with readable real `x`, `y`, `z`
attributes (including RDKit `Point3D`). Lists, tuples, strided/read-only arrays
and non-native byte order are supported. Annotations describe the sequence
and structural point alternatives, numeric array input, and a float64 array
result; they do not use a bare `object` input. Generators, strings, booleans,
complex numbers, object arrays and array subclasses such as masked arrays
are rejected. Inputs are copied and never modified.

Points must be paired in the same order and units. There is no atom matching,
chemical validation or conformer generation. The number of points is between
1 and 65,536, checked before copying input arrays. Both clouds must have the
same size. Each row has three finite real values with magnitude ≤1e100.
This is a deliberate supported numeric range, not an arbitrary-precision
geometry guarantee.

Weights are positive real pair weights with the same magnitude bound; `None`
or an empty sequence/ordinary real array means equal weights. Wrong lengths,
zero/negative weights, nonfinite values and weight ratios that underflow after
normalization raise errors. Weight normalization is used only in fitting;
the returned SSD uses the original weights. Scaling all weights by `k` scales
SSD by `k`, while leaving the fitted transform unchanged.

With `reflect=False`, the rotation determinant is +1. With `reflect=True`, an
improper rotation (determinant −1) is required. The option does not mean
"choose whichever reflection gives the smaller SSD". The homogeneous matrix
acts on column vectors, with the bottom row `[0, 0, 0, 1]`. SSD is the weighted
sum of squared residual lengths, not RMSD or mean squared error. Floating-point
roundoff can produce a small positive SSD for an exact fit.

The signature retains `maxIterations=50`, but this SVD backend has no
iterative-eigensolver budget: only that default is accepted, and other values
raise `ValueError`. Non-boolean `reflect` and non-integer `maxIterations`
raise `TypeError`. No supplied option is silently ignored. Degenerate point
clouds are supported; their optimal rotation may be nonunique. Failed SVD or
unrepresentable results raise `ValueError`, without a partial transform.

The return direction, SSD and reflection convention follow the
[official RDKit alignment API](https://www.rdkit.org/docs/source/rdkit.Numerics.rdAlignment.html).
The implementation is independent of RDKit and uses weighted centroids,
scaled covariance, and an SVD with the requested determinant sign. Raw NumPy
SVD output is not compared for exact matrix equality with RDKit's quaternion
solver.

## Validation

Local validation with NumPy 2.5.3 and RDKit 2026.09.1 passed 576 reference
comparisons and 1,152 scaling checks. The maximum absolute SSD difference was
2.39e-10. For nondegenerate examples, the maximum elementwise transform
difference was 3.65e-5; it is retained as a solver-tolerance diagnostic, not
an exact-transform compatibility claim. Sixty-two new tests cover container
forms, weighting, reflections, rank-deficient clouds, numeric/shape errors,
input immutability and the point limit; the selected Python regression suite
passed 224 tests. Type-check evidence includes valid inputs and deliberately
invalid point, weight and option types.

Reproduction commands, source and runtime-fixture hashes are in
`validation/results/point-alignment-2026-10-17/`. This Python/NumPy change used
a reused native extension as an import fixture. No fresh native build, wheel
publication, Rust alignment change or complete RDKit wrapper parity is claimed.
