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

## Align multiple conformers and collect RMS values

`chematic.align_conformers_checked` extends the same solver to coordinate
ensembles and addresses the `AlignMolConformers` RMS-output item in
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601). It accepts explicit
coordinates rather than a molecule holding mutable conformers.

```python
from chematic import align_conformers_checked

rms = [12.5]  # existing results are retained
aligned = align_conformers_checked(
    [
        [[0., 0., 0.], [1., 0., 0.], [0., 2., 0.]],
        [[3., 4., 1.], [3., 5., 1.], [1., 4., 1.]],
    ],
    atom_ids=[0, 1, 2],
    conformer_ids=[0, 1],
    weights=[1., 2., 3.],
    rms_list=rms,
)
# aligned: a new float64 (2, 3, 3) array; rms: [12.5, approximately 0.0]
```

Every frame has the same atom count, atom order and units. The first selected
frame is the reference, and the remaining selected frames are fitted to it in
selection order. Fitting on an atom subset still applies the resulting rigid
transform to **all** atoms of that frame. Reference and unselected frames retain
their coordinates. Frame indices refer to the input array, not external molecule
conformer IDs. None or empty selections mean all atoms/frames. Duplicate and
out-of-range indices raise. Only proper rotations are supported; mirror fitting,
atom correspondence search and conformer generation are outside this API.

The result is a separate array: caller coordinates remain unchanged, including
on failure. `rms_list` is annotated `Optional[list[float]]` and must be an ordinary
list containing finite Python floats. One RMS float per selected probe frame is
appended **only after every fit succeeds**. Tuple/array destinations, list
subclasses and invalid existing entries raise before fitting. There is no bare
`object` output parameter or acceptance of arbitrary `append` methods. A single
selected frame adds no RMS entries, but inputs and weights are still validated.

For selected atoms A and positive weights w, RMS is
`sqrt(sum(w[a] * squared_distance[a] for a in A) / len(A))`. This follows
[RDKit's implementation](https://github.com/rdkit/rdkit/blob/Release_2026_09_1/Code/GraphMol/MolAlign/AlignMolecules.cpp#L420-L478):
it divides by atom count, not the sum of weights. Multiplying all weights by k
retains the fit and multiplies RMS by sqrt(k). Default/empty weights mean ones.
The inherited point and weight validation applies. Ensembles are limited to
1..1024 frames, 1..65536 atoms per frame and one million coordinate points in
total, checked before output allocation. Coordinate magnitudes must be finite
and at most 1e100. Ordinary numeric arrays (including strided, read-only and
non-native-endian arrays) and point-cloud sequences are supported; booleans,
complex/object arrays, generators and ragged shapes are rejected.

On 2026-10-17, 2,592 seeded cases compared 9,072 frames and 4,320 appended RMS
values with RDKit 2026.09.1, varying atom selection/order, frame selection/order,
and weights. RMS agreed within rtol=atol=1e-7 (largest absolute difference
2.30e-10), and every case preserved intra-frame distances within 1e-12. **Output
coordinates are not claimed identical:** the SVD and RDKit iterative solvers
produced a largest component difference of 2.34e-4 in these generated coordinates.
For underdetermined fits, rotations can differ while residuals and distances
agree; one-, two- and three-atom diagnostics check those quantities. These finite
samples establish neither universal numerical parity nor chemistry or speed.

Regression tests cover atomic output on a later solver failure, rejected
inputs, weighted RMS scaling, atom subsets, reference ordering and unchanged
inputs. Reproduce the independent comparison with
`python validation/validate_conformer_alignment.py` with RDKit 2026.09.1
installed. Evidence lives in `validation/results/conformer-alignment-2026-10-17/`.
This extension of the Python/NumPy API reuses the existing native extension as
an import fixture; it is not a fresh native build of the candidate.
