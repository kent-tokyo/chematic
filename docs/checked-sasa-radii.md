# Checked SASA with caller-supplied radii

RDKit [#9601](https://github.com/rdkit/rdkit/issues/9601) and its open
[PR #9611](https://github.com/rdkit/rdkit/pull/9611) include the requirement that
FreeSASA receives one radius per atom. Chematic's existing SASA entry points
use built-in element radii. The new API adds explicit coordinates and per-atom
radii with checked counts, values and computation limits.

```python
import chematic

mol = chematic.from_smiles("CO")
result = mol.sasa_with_radii_checked(
    [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
    [1.0, 1.5],
    probe_radius=0.0,
    sphere_points=8192,
)
print(result["total"], result["per_atom"])
```

Rust exposes `sasa_with_radii_checked(mol, coordinates, radii, probe_radius,
sphere_points) -> Result<CheckedSasa, SasaInputError>`. The result contains
`total` and `per_atom`; the total uses compensated summation.

Coordinates and positive van der Waals radii are in Angstrom; areas are in
square Angstrom. Exactly one coordinate triple and one radius must match each
atom in insertion order, including any explicit H atoms. Implicit hydrogens
are not materialized. Radius inflation adds the nonnegative probe radius
once. There is no coordinate generation, radius classification or default
radius substitution. Empty molecules with empty arrays return zero area.

Nonfinite values, nonpositive radii, invalid counts or unsupported computation
sizes produce errors with no partial result. Python raises `ValueError` for
these inputs, `TypeError` for wrong types (shape errors can also raise
`ValueError`), and propagates iterator exceptions. Both iterables consume at
most one excess item to detect a count mismatch. The point-count argument is
a signed 64-bit Python conversion; values outside it raise `OverflowError`.

## Sampling and numerical boundary

The algorithm samples midpoint Fibonacci directions on inflated spheres and
checks exposure against overlapping neighbors. Fully contained spheres have
zero accessible area; disjoint spheres retain their full area. Tangent
boundaries have zero measure. Equal inflated radii at exactly coincident
centers are refused because assigning per-atom ownership of the shared surface
is ambiguous. Concentric spheres with unequal radii are supported.

At most 2048 atoms and 32..65536 points per sphere are accepted. The neighbor
preparation scans at most the bounded atom-pair set. Partially overlapping
directed neighbor pairs multiplied by point count must not exceed 50 million;
this conservative cap can refuse inputs even when early exits would do less
work. It is not a wall-clock timeout guarantee.

Occlusion uses local coordinate differences normalized by a common radius
and displacement scale. It does not construct absolute translated surface
points, and uses robust norms for finite coordinates whose subtraction may
overflow when spheres are far apart. Areas or totals outside representable
finite float range are errors. Very small nonzero sampled areas that round to
zero are also errors.

Sampling still has finite-resolution error and orientation dependence. This
is a native Shrake-Rupley implementation, not the FreeSASA backend or a
bit-exact RDKit API. It does not provide FreeSASA protein radius classification,
polar/apolar queries, Lee-Richards integration, atom properties or implicit-H
expansion. Existing SASA functions remain unchanged; no WASM binding is added.

The reference API's radius-indexing contract is documented in
[RDKit rdFreeSASA](https://www.rdkit.org/docs/source/rdkit.Chem.rdFreeSASA.html).
The independent probe uses analytic spherical-cap areas and reports RDKit
Lee-Richards values separately, with explicit tolerances instead of claiming
that different integration algorithms return identical numbers.

On the dated macOS arm64 probe, 96 two-sphere geometries, their atom-order
reversals and translations passed the analytic checks (288 native runs).
The largest per-atom error was 0.159% of the corresponding full isolated
sphere area, below the 0.5% threshold. Five scale cases from `1e-100` to
`1e100` also passed. These are bounded geometry controls, not a general
protein or molecular accuracy claim.

RDKit 2026.09.1's default Lee-Richards diagnostic was finite and within its
separate 2.5% full-sphere tolerance in 95 cases, and returned a nonfinite
result in one case. It also accepted a two-atom/one-radius input in an isolated
subprocess on this host; no crash was observed. Both observations are retained
in the evidence rather than treated as an accuracy reference.

Reproduction and dated source evidence are in
`validation/rdkit_issues/sasa_checked/probe.py` and
`validation/results/sasa-checked-2026-10-17/`.
