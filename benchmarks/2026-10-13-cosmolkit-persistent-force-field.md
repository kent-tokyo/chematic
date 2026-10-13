# COSMolKit 0.5 feature-gap follow-up: persistent force fields

Date: 2026-10-13
CheMatic base: `origin/main` at `1b761175`
Comparator: COSMolKit 0.5.0 public Python API and source documentation
Host: macOS arm64, CPython 3.13.6

## Why this feature

COSMolKit 0.5 exposes a mutable `MolecularForceField` that retains prepared
MMFF/UFF terms. CheMatic already had public Rust `RdkitMmffField` and
`RdkitUffField` implementations, but its Python energy and minimization methods
rebuilt those fields for every call. The missing public lifecycle was therefore
a bounded API and performance gap rather than missing chemistry.

This change adds Python `Mol.mmff_force_field()` and
`Mol.uff_force_field()`. The returned `MolecularForceField` owns mutable
coordinates and supports energy, gradient, combined energy/gradient, coordinate
updates, and in-place minimization. Existing one-shot methods are unchanged.

The implementation is independent. COSMolKit source code was not copied;
COSMolKit's public API and documentation were used only to identify the product
gap. CheMatic's existing force-field implementations remain the computational
core.

Primary references:

- <https://pypi.org/project/cosmolkit/0.5.0/>
- <https://kit.cosmol.org/python/api>
- <https://github.com/cosmol-studio/COSMolKit>

## Correctness gate

`crates/chematic-py/tests/test_rdkit_forcefield.py` compares persistent MMFF94,
MMFF94s, and UFF energies with the existing one-shot RDKit-compatible APIs at
the same coordinates. It also covers gradients, coordinate replacement,
minimization, invalid lengths, non-finite coordinates, out-of-range atoms, and
invalid variants or thresholds.

The candidate wheel passed all eight focused force-field tests. This establishes
equivalence to CheMatic's existing pinned force-field contract; it is not new
experimental accuracy evidence.

The Python package suite passed 1,193 tests after deselecting the pre-existing
caffeine HBA oracle mismatch (`3` versus RDKit `6`). The unrestricted run had
the same 1,193 passes and that single unrelated failure. The binding crate's
five Rust tests, `cargo clippy -D warnings`, formatting, and diff checks also
passed.

## Exploratory reuse timing

One process performed 1,000 MMFF94 energy calls on explicit-hydrogen
acetanilide coordinates from the existing force-field fixture:

| path | elapsed |
|---|---:|
| one-shot `Mol.rdkit_mmff_energy()` | 55.33 ms |
| retained `MolecularForceField.energy()` | 1.38 ms |

The observed ratio was 39.99x and both paths returned the same energy. This is
a focused demonstration that preparation is actually reused. It is not a
stable speed claim: there are no alternating blocks, confidence interval,
second host, memory result, or published artifact yet.

## Remaining gap

COSMolKit's object also exposes fixed-atom constraints and JavaScript coverage.
CheMatic's first slice does not. Those remain follow-up work; this record does
not claim complete `MolecularForceField` parity.
