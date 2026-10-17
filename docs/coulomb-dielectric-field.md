# Checked local-solvation Coulomb field

`chematic.rdkit_compat.rdMIF.CoulombDielectric` evaluates explicit charges
with RDKit's local-solvation dielectric model. It extends the typed
[vacuum Coulomb field](coulomb-interaction-field.md) and partially addresses
the `rdMIF` `charges` / `positions` typing work in
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601).

```python
from chematic.rdkit_compat import rdMIF

field = rdMIF.CoulombDielectric(
    [1.0, -0.3], [(0., 0., 0.), (1., 0., 0.)], epsilon=80., xi=4.
)
energy = field(2., 0., 0., 25.)  # kJ/mol; threshold in squared angstroms
```

`epsilon` is the solvent relative permittivity (default 80); `xi` is the
solute relative permittivity (default 4). Charges are in elementary-charge
units, Cartesian positions and `cutoff` are in angstroms, and `softcoreParam`
is in squared angstroms. Construction copies indexed inputs, including
ordinary numeric NumPy arrays, into an immutable vacuum-field snapshot.
Later mutations of caller data do not change the field.

## Model and numerical contract

The reference is the explicit-charge constructor and scalar call in
[RDKit Release_2026_09_1 MIFDescriptors.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_09_1/Code/GraphMol/MolInteractionFields/MIFDescriptors.cpp).
This is that local empirical model, not a general dielectric boundary solver
or a new physical calibration.

Each site counts other sites within distance strictly less than 4 A. The
probe counts **all** sites with squared distance strictly less than 16 A²,
including zero-charge sites and sites outside the energy threshold.
For a neighbor count `n`, the empirical depth `s(n)` is:

| Neighbors | 0–6 | 7 | 8 | 9 | 10 | 11 | 12+ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Depth (A) | 0 | 0.4 | 0.9 | 1.4 | 1.9 | 2.6 | 4 |

For site `i`, let `r2` be its squared distance from the probe and let
`u = max(r2, cutoff²)` in hard-core mode (`softcoreParam=0`). In soft-core
mode, `u = r2 + softcoreParam`. Only sites with **`r2 < threshold`** contribute:

```text
contrast = (xi - epsilon) / (xi + epsilon)
term_i = charge_i * (1/sqrt(u) + contrast/sqrt(u + 4*s_i*s_probe))
energy = sum(term_i) * 1389.048632870505 * (1/xi) * probeCharge
```

The sum follows input site order. The prefactor is the pinned legacy RDKit
conversion constant. `absVal=True` returns negative absolute **total** energy.
If `epsilon == xi`, contrast is zero and the result is the vacuum field
scaled by `1/xi`. In sparse neighborhoods with zero site depth it reduces
to the vacuum field scaled by `2/(xi+epsilon)`.

Python floating-point evaluation and compiled C++ arithmetic can disagree at
rounded threshold or 4 A neighbor boundaries. Constructor site counting
uses `sqrt(dx*dx + dy*dy + dz*dz) < 4`; probe counting uses the squared
comparison. Exact boundary-inclusion parity across implementations is not
promised. Choose configurations with a distance margin if stable inclusion
matters.

## Input and work limits

- At most **1,024 sites**, checked before indexing/copying either input.
  Site-neighbor counting takes O(N²) time; a query takes O(N) time and local
  scratch memory. No per-query scratch is retained on the field.
- `epsilon` and `xi` must be finite real values in **1..1e6**. Zero,
  negative values, values below one, infinity and NaN are refused.
- Other [vacuum-field input rules](coulomb-interaction-field.md#checked-input-limits)
  apply: matching lengths, exactly three coordinates, bounded finite real
  scalars, strict Boolean `absVal`, positive hard cutoff and supported soft
  core. Its larger vacuum-site cap does not apply to this class.
- Empty fields are allowed and return zero (negative zero with `absVal=True`).
  Type errors raise `TypeError`; length/range/nonfinite violations raise
  `ValueError`; source getter errors propagate.

Molecule/conformer overloads, charge assignment, grid generation and arbitrary
solvent geometry remain outside this profile. This is standard-library Python,
with no new native/Rust/WASM binding or NumPy runtime requirement.

## Local evidence — October 17, 2026

Before the follow-up at PR commit `5473b7de`, `rdMIF` exposed only `Coulomb`.
Baseline main `35b932f9` had no `rdMIF`. The final source passed **284 tests**:
59 new dielectric tests, 63 vacuum tests and 162 existing compatibility tests.
Analytical tests cover every neighbor-table transition, threshold-excluded
zero-charge neighbors, equal permittivities, signs, absolute totals, copied
inputs, concurrent/repeated evaluation and the 1,024-site cap.

RDKit **2026.09.1** comparison: 3,168 constructors, **154,080 ordinary
scalar evaluations** and **72 exact/adjacent 4 A boundary checks** passed
with relative tolerance 2e-12 and absolute tolerance 2e-9 kJ/mol. Maximum
absolute error was `9.313225746154785e-10 kJ/mol`. A separate 85,536
near-threshold diagnostics had zero differences above tolerance; this
synthetic result does not establish universal boundary parity. The earlier
vacuum-only comparison retains its separately documented 108 differences.

Both helpers passed strict mypy in an isolated package; a valid public
consumer passed and four invalid calls produced four expected type errors.
Consumer checking follows imports silently; package-wide strict cleanliness
is not claimed. Source/stub/binding/unsafe/product-name checks and strict
MkDocs passed.

Reproduce with the candidate package and RDKit 2026.09.1 available:

```sh
python3.13 validation/validate_coulomb_dielectric.py \
  --output /tmp/coulomb-dielectric-comparison.json
```

Evidence is in `validation/results/coulomb-dielectric-2026-10-17/`. The existing
native extension was reused for package import and compatibility tests; no
fresh native build, released wheel, speed improvement or upstream fix is
represented. The original vacuum evidence remains pinned to its earlier
commit; this follow-up has its own hashes and commands.
