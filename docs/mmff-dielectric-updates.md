# Updating a prepared MMFF dielectric constant

Open Babel [#2811](https://github.com/openbabel/openbabel/issues/2811) reports
that changing the dielectric constant after force-field setup leaves cached
electrostatic terms unchanged. An original Na⁺/Cl⁻ pair at 4 Å reproduces the
problem on Open Babel 3.2.1: the same object reports dielectric constants
1, 4 and 78 but returns the same electrostatic energy. Repeating `Setup` on
that object also leaves the old terms. Fresh instances respect the setting.

chematic's prepared MMFF force field previously offered coordinate updates
but no dielectric-setting API. Rust and Python now support changing a
**constant relative dielectric** on an existing MMFF94 or MMFF94s field.

```python
mol = chematic.from_smiles("CCO").add_hydrogens()
field = mol.mmff_force_field(coords)  # one XYZ per explicit atom, in Å
assert field.dielectric_constant == 1.0

field.set_dielectric_constant_(4.0)
energy, gradient = field.energy_gradient()
field.set_dielectric_constant_(78.0)
status, minimized_energy = field.minimize_(200)
```

Rust exposes `RdkitMmffField::dielectric_constant()` and
`set_dielectric_constant(value)` under `chematic_ff::rdkit_mmff`. Errors use
the new non-exhaustive `MmffDielectricError` type, leaving existing public
error enums intact. Python's setter returns the same object for chaining.

## Numerical contract

For an included atom pair, the electrostatic contribution in kcal/mol is

\[
E_{ij}=\frac{332.0716\,q_iq_j}{\epsilon\,(r_{ij}+0.05)}\,s_{ij},
\]

where distance is in Å, \(s_{ij}=0.75\) for a 1–4 pair and 1 otherwise.
The gradient differentiates this same screened energy. The convention follows
RDKit 2026.09.1's
[MMFF nonbonded implementation](https://github.com/rdkit/rdkit/blob/Release_2026_09_1/Code/ForceField/MMFF/Nonbonded.cpp).

- The default remains 1.0. Default evaluation keeps the existing arithmetic.
- Updates recompute each screened charge product from its **original** product,
  so repeated changes do not compound scaling or accumulate division error.
  Evaluation uses the prepared coefficient; it performs no new per-pair division.
- Atom typing, partial charges, bonded terms, van der Waals parameters and the
  prepared nonbonded pair list stay fixed. Coordinate updates and minimization
  use the current dielectric. Changing dielectric is not a new setup operation.
- Only finite positive constants are accepted. A value that overflows a
  prepared electrostatic coefficient is refused before any state changes.
  Rust returns an error and Python raises `ValueError`; the old constant,
  coordinates, energy and gradient remain available.
- UFF has no electrostatic term in this API. Its dielectric getter and setter
  raise `ValueError`. Empty electrostatic term lists are valid in MMFF.

This API does not add a distance-dependent dielectric, implicit solvent,
periodic electrostatics, solvent calibration or GAFF. The separate native
`Mmff94EnergyModel`, stateless energy/optimization calls, WASM pipeline and
UFF retain their existing defaults and interfaces.

## Validation and reproduction

`validation/openbabel_issues/mmff-dielectric/cases.json` contains seven authored
molecular/ionic cases with seed-42 RDKit 2026.09.1 coordinates and 56 constant-
dielectric energy/gradient references across MMFF94 and MMFF94s. Rust and
Python tests change one prepared object through 1→4→78→0.5→4→1, verify default
restoration, rejected updates, finite-difference gradients, unaffected other
terms, coordinate changes and minimizer reuse. Default MMFF regression tests
remain separate from these new-setting checks.

`scripts/validation/validate_mmff_dielectric.py` rechecks those goldens with pinned RDKit,
then compares fresh seed-7/19/71 conformers and displaced ionic pairs. Each
native field is reused across dielectric changes; independent RDKit reference
fields are constructed with `SetMMFFDielectricModel(1)` and the requested
constant. Ten-step minimizations compare status, energy and coordinates;
the reused native field is also checked against a fresh native field.

The original `openbabel_probe.cpp` compares live updates, repeated `Setup`
and fresh Open Babel instances. Compile it against the installed Open Babel
headers/library and pass the executable through `--openbabel-probe`.
The validator retains all starting coordinates, gradients, minimization
outputs, diagnostic stdout/stderr and source/runtime hashes.

The local source packet is `validation/results/mmff-dielectric-2026-10-16`:

- All 56 checked-in energy/gradient goldens recheck against RDKit 2026.09.1.
- All 252 new energy/gradient comparisons pass. Maximum absolute differences
  are 4.27×10⁻¹⁴ kcal/mol and 4.11×10⁻¹³ kcal/mol/Å respectively.
- All 126 reused-versus-fresh native minimizations agree exactly.
  Against RDKit, all statuses agree and maximum energy difference is
  9.84×10⁻⁹ kcal/mol. Two Na⁺/Cl⁻ cases (seed 7, dielectric 4, both variants)
  exceed the separate 5×10⁻⁶ Å coordinate diagnostic threshold; the maximum
  coordinate difference is 1.51×10⁻⁵ Å. The remaining 124 meet that threshold.
  These two residuals are retained in the summary and raw packet; they do not
  establish exact RDKit minimizer trajectory parity.

These are local source results for the stated fixtures and settings. They do
not establish general solvent accuracy, force-field coverage, speed, package
publication or an upstream Open Babel fix.
