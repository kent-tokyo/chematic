# Explicit-charge Coulomb interaction field

`chematic.rdkit_compat.rdMIF.Coulomb` constructs a checked vacuum Coulomb
field from caller-supplied charges and Cartesian positions. It addresses the
remaining `rdMIF` `charges` / `positions` typing boundary in
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601), which remains an
upstream tracking issue. This is a partial implementation of that request.

```python
from chematic.rdkit_compat import rdMIF

field = rdMIF.Coulomb([1.0, -1.0], [[0., 0., 0.], [1., 0., 0.]])
energy = field(2.0, 0.0, 0.0, 25.0)  # kJ/mol; squared-distance threshold
assert field.site_count == 2
```

Charges use elementary-charge units; positions and probe coordinates use
angstroms. List, tuple and numeric indexed NumPy inputs are accepted. Each
position must have exactly three real values and the number of positions
must equal the number of charges. Construction copies all scalars into
immutable tuples: later input mutation does not change the field. The helper
uses the Python standard library; NumPy is optional.

## Calculation and compatibility profile

For each site let `r2 = dx*dx + dy*dy + dz*dz`. Only sites satisfying
`r2 < threshold` contribute. **`threshold` is a squared distance in A²**:
this matches the pinned scalar implementation and direct measurements, even
though the upstream Python documentation describes it as a distance.
A site at distance 2 contributes above threshold 4, and is excluded at 4.

The default `softcoreParam=0` uses denominator
`sqrt(max(r2, cutoff*cutoff))`, with default `cutoff=1`. A positive
`softcoreParam` uses `sqrt(r2 + softcoreParam)` and ignores the hard floor.
The ordered sum of `charge / denominator` is multiplied by the probe charge
and RDKit 2026.09.1's legacy unit conversion constant
`1389.048632870505`. This constant is retained for comparison with that
version; it is not a claim to use a newly evaluated physical constant.
`absVal=True` returns the **negative absolute value of the total**.

The reference is
[RDKit Release_2026_09_1 MIFDescriptors.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_09_1/Code/GraphMol/MolInteractionFields/MIFDescriptors.cpp).
Python and compiled C++ can round squared distances differently. Near a
site's threshold boundary this can change whether the site contributes;
bit-for-bit inclusion parity at such boundaries is not promised. The Python
contract uses the ordered floating-point expression above. Choose a threshold
with a margin from site distances when stable inclusion matters.

## Checked input limits

- At most 100,000 sites, checked before indexing or copying; empty inputs are
  allowed. Construction and evaluation each take O(number of sites).
- Charges and probe charge must be finite real numbers in ±1e6. Coordinates
  must be finite in ±1e100. Boolean and complex scalars are refused.
- `threshold` must be finite in 0..1e200. Zero selects no sites.
- `cutoff` must be finite in 0.001..1e100, including in soft-core mode.
  Smaller values are refused instead of reproducing RDKit's tiny-cutoff fallback.
- `softcoreParam` must be zero or finite in 1e-6..1e200 A². Negative values
  and positive values below 1e-6 are refused instead of reproducing undefined
  or silently ignored settings.
- `absVal` must be a Python `bool`. Invalid scalar/sequence types raise
  `TypeError`; numeric bounds, nonfinite values and length errors raise
  `ValueError`. This intentionally tightens the upstream non-sequence
  `ValueError` boundary described in #9601. Source sequence getter errors
  propagate.

This profile has no molecule/conformer overload, automatic charge lookup or
assignment, dielectric model, grid output or batch acceleration. It does not
change the existing native charge calculators. Numerical/range requirements
remain runtime checks; type annotations cannot prove them.

## Local evidence, October 17, 2026

Baseline commit `35b932f9` did not expose `rdMIF`. The candidate passed
225 Python tests (63 new and 162 existing compatibility tests), strict helper
mypy, and public-consumer typing checks: valid usage passed and four invalid
arguments were rejected.

The seeded comparison against RDKit **2026.09.1** covered 864 constructors,
list/tuple/strided NumPy inputs, forward/reversed site orders, probe signs,
hard floors, soft cores, total absolute values and empty inputs. All 120,960
ordinary-point comparisons passed `rel_tol=2e-12`, `abs_tol=2e-9 kJ/mol`;
maximum absolute error was `4.656612873077393e-10 kJ/mol`. A separate 88,128
near-boundary diagnostic evaluations had 108 differences above that tolerance.
Every such difference occurred within eight floating-point ULP of a site
boundary; examples are retained rather than counted as passing parity tests.
These counts are repeated synthetic evaluations, not independent molecules.

Reproduce with the candidate Python package and pinned RDKit on `PYTHONPATH`:

```sh
python3.13 validation/validate_coulomb_field.py \
  --output /tmp/coulomb-comparison.json
```

The checked-in evidence directory is
`validation/results/coulomb-field-2026-10-17/`. The existing native extension
was reused for package import and existing compatibility tests; Coulomb
calculation itself is standard-library Python. No fresh native build, native
binding change, release, upstream fix or performance advantage is claimed.
