# MMFF with caller-supplied partial charges

[Open Babel #2496](https://github.com/openbabel/openbabel/issues/2496) asks to
use externally calculated partial charges in force-field calculations. Its
report concerns GAFF. Chematic adds the transferable input capability to its
prepared **MMFF94/MMFF94s** fields.

```python
import chematic

mol = chematic.from_smiles("[Na+].[Cl-]")
field = mol.mmff_force_field(
    [[0, 0, 0], [4, 0, 0]],
    ignore_interfrag_interactions=False,
    partial_charges=[0.3, -0.8],
)
energy, gradient = field.energy_gradient()
field.set_dielectric_constant_(4.0)
status, energy = field.minimize_(100)
```

The optional keyword `partial_charges` accepts one finite value in units of
elementary charge for each graph atom, in insertion order, including explicit
hydrogens. `None` retains the existing MMFF charge assignment. A charge array
may come from a caller-parsed file or another charge model; this API does not
read or infer a file's atom correspondence. Add explicit H atoms before
assigning the array if the modelling workflow requires them.

Rust exposes `RdkitMmffField::with_terms_and_partial_charges(mol, coords,
variant, non_bonded_thresh, ignore_interfrag, terms, charges)` and the separate
`MmffPartialChargeError`. It checks charge/coordinate counts, finite values,
the finite nonnegative cutoff and representable included-pair coefficients.
Python reports these errors as `ValueError` and wrong types as `TypeError`.
Charge iterables consume at most atom count plus one items, and preserve
iterator exceptions; generators are supported.

## Model and lifecycle

Only electrostatic charge products are replaced. MMFF atom typing, bonded
terms and vdW parameters remain unchanged. Each eligible pair contributes

\[
E_{ij} = 332.0716\,q_iq_j\,s_{ij}/[\epsilon(r_{ij}+0.05)]
\]

in kcal/mol, using Å distances, with `s = 0.75` for 1–4 pairs and 1 otherwise.
1–2/1–3 exclusions, the interfragment option and the construction-time distance
cutoff apply. The convention follows the existing
[MMFF nonbonded expression](https://github.com/rdkit/rdkit/blob/Release_2026_09_1/Code/ForceField/MMFF/Nonbonded.cpp).
Previously zero-charge pairs are considered during custom preparation.

Exact zero suppresses an atom's electrostatics. Explicit small nonzero charges
are retained, including values below the default MMFF near-zero cutoff.
Underflow to zero or overflow of an included charge product/coefficient is
refused before returning a field. These checks do not promise representable
energy/gradient values for every possible geometry or sum of terms.

No net-charge adjustment, normalization or check against formal charge is
performed. The graph and its formal charges/properties are not modified.
Coefficients are owned by the prepared field: changing the original Python
list does not change the field. Coordinate changes, dielectric updates and
minimization use the same supplied charges and prepared pair set. Rebuild the
field to change charges or pair selection.

This is an opt-in hybrid charge model, outside standard MMFF parametrization;
the API does not establish its suitability for any physical system. It adds
no GAFF, charge-file parser, automatic correspondence, UFF electrostatics or
WASM binding. Existing default and stateless MMFF calls retain their behavior.

The independent probe compares electrostatic energy and gradient differences
against a graph-BFS/analytic calculation, and differentiates the full hybrid
energy numerically. An authored methanol control also reproduces Open Babel
3.2.1 GAFF accepting different preassigned atom charges but returning unchanged
energies. Reproduction and dated results are under
`validation/openbabel_issues/mmff-partial-charges/` and
`validation/results/mmff-partial-charges-2026-10-17/`.
