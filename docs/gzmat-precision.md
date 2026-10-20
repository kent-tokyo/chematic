# Gaussian-style Z-matrix precision

This local candidate addresses the geometry-conversion request in
[Open Babel #2794](https://github.com/openbabel/openbabel/issues/2794).
It adds a checked Cartesian-to-Z-matrix conversion and an explicitly
configurable precision writer. Existing Cartesian Gaussian input/output APIs
remain compatible. This document describes source validation, not a release
or an upstream issue closure.

## API and units

Rust exports from `chematic_mol`:

- `cartesian_to_zmatrix(elements, positions_angstrom, references)` returns
  `Vec<ZMatrixRow>`.
- `zmatrix_to_cartesian(rows)` reconstructs a canonical Cartesian frame.
- `validate_zmatrix(rows)` checks local row structure, finite values and ranges.
  Reconstruction additionally checks reference planes.
- `write_gzmat(rows, &GzmatWriteOptions)` validates and writes the document,
  including checking reconstruction after requested rounding.

`ZMatrixRow` has an `Element`, zero-based `references: Vec<usize>`, and
`values: Vec<f64>`. The first row has no references; subsequent rows must
have exactly `min(row_index, 3)` distinct earlier references. The reference
sequence defines a distance, angle, then signed torsion. Distances are in
angstroms and angles in degrees. These geometric references do not imply
chemical bonds. Atom ordering is retained; no reference search or reordering
is performed.

Python exposes the same conversion and direct writing from Cartesian arrays:

```python
import chematic

elements = ["C", "N", "O", "F"]
positions = [
    [0.0, 0.0, 0.0],
    [1.234567891234567, 0.0, 0.0],
    [1.1, 1.33333333333333, 0.0],
    [0.1, 0.3, 1.444444444444444],
]
references = [[], [0], [1, 0], [2, 1, 0]]

rows = chematic.cartesian_to_zmatrix(elements, positions, references)
# rows: [{"element": str, "references": list[int], "values": list[float]}, ...]
text = chematic.write_gzmat(
    elements, positions, references,
    charge=0, multiplicity=1, significant_digits=17,
    title="Explicit internal geometry", route="#",
)
```

Python accepts canonical elemental symbols, not isotope names, dummy labels
or wildcard atoms. Reference indices, charge, multiplicity and precision must
be integers; booleans are rejected. Reference iterators are consumed with
row/count bounds. Returned rows are detached dictionaries. Charge and spin
are supplied by the caller, not inferred; physical electron-count/parity
compatibility is not checked.

## Precision contract

`GzmatWriteOptions::new(charge, multiplicity)` defaults to 17 significant
decimal digits, the title `chematic internal coordinates`, and route `#`.
Precision can be set from 1 through 17. Output uses scientific notation for
every internal value, including very small distances, and one distinct
`rN`, `aN` or `dN` variable per coordinate. File references are one-based.
A one-atom document needs no variable section.

At 17 digits, parsing each serialized scalar as `f64` gives the computed
internal value exactly. This is a serialization guarantee; conversion between
Cartesian and internal coordinates still has floating-point error. Internal
coordinates retain shape and handedness but discard absolute translation and
rotation. Lower precision can change shape, even when the frame remains valid.

The route is a single-line Gaussian-style header supplied by the caller. The
`#` default is a template marker; select calculation keywords before running
an actual job. Route semantics, including any unit-changing keywords, are the
caller's responsibility. Length values emitted by this API are angstroms.
Title and route fields are nonempty single lines of at most 1024 bytes;
control characters and section injection are rejected.

## Checked failures and scope

The converter uses finite vector differences, scaled norms and `atan2` for
angles and signed torsions. It rejects coincident references, overflowing
lengths and undefined/near-linear coordinate or reference planes. The minimum
accepted sine is strictly greater than `1e-10`. It does not insert a dummy
atom or choose an arbitrary torsion to hide a singularity.

The writer checks public row data again. It reconstructs the rounded rows
before writing, so rounding that creates a zero/linear angle, collapses a
later reference plane, or overflows coordinates fails. Nondegenerate planar
geometries and zero/180-degree torsions are valid. Linear bond angles are
outside this conservative profile. Work and storage are linear in atom
count, with a hard cap of 100,000 atoms.

This profile contains elements and geometry, not a molecular graph. It does
not serialize isotope, topology, atom-map, stereogroup, periodic, ghost-atom
or dummy-atom metadata. There is no automatic reference selection, general
Gaussian Z-matrix parser, ORCA internal-input writer or quantum job execution.
The existing ORCA adapter continues to preserve raw internal-coordinate input.

Syntax and geometric conventions were implemented independently from the
[Open MOPAC Gaussian Z-matrix specification](https://openmopac.net/Manual/gaussian_z.html)
and checked against the
[ORCA 6.1 coordinate definitions](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html).
No Open Babel GPL implementation or third-party issue attachment was copied.

## Independent validation

The original fixture manifest is
`validation/openbabel_issues/gzmat/precision-geometries.json`.
It records provenance, license, minimization, payload hash, pinned oracles,
settings and tolerances. Seven one- to five-atom fixtures cover signed and
planar geometries, a mirror pair, and non-chain references.

The validator changes atom order (up to 24 permutations), reference pattern,
proper rotation/translation and significant digits (8, 12, 17). RDKit
2026.09.1 independently calculates distances, angles and torsions from the
input Cartesian conformer. Open Babel 3.2.1 reads the resulting GZMAT with
bonding disabled, then emits CJSON with double-precision coordinates. The
validator checks indexed elements, explicit charge/spin, all pair distances,
and nonzero signed tetrahedron volumes; it does not rely on alignment RMSD
that could hide a reflection.

Results are recorded in
`validation/results/openbabel-gzmat-precision-2026-10-15.json`.
The report pins source/fixture/extension/binary hashes and separates exact
scalar serialization from Cartesian reconstruction tolerances. Open Babel
readback is measured; Gaussian and ORCA execution have not been tested.

The run passed **1,242** cases. It independently verified **2,016** nonzero
signed tetrahedron volumes and exact serialization of **2,668** computed
scalars at 17 digits. Maximum pairwise-distance errors were `1.185e-7 Å`
at 8 digits, `1.614e-11 Å` at 12 digits and `8.882e-16 Å` at 17 digits.
The maximum independent RDKit scalar discrepancy was `2.843e-14`
(angstroms or degrees according to the value). These measurements apply to
the seven synthetic fixtures and the stated transforms, not all geometries.

Local Rust validation passed **753 library + 13 new integration tests**;
Python passed **43 new + 38 existing TURBOMOLE tests**. Clippy with warnings
denied, formatting, strict documentation, refactoring budget and binding
inventory checks passed.

An independent synthetic precision probe also reproduces the upstream
limitation: Open Babel's default output and `-x8` output are identical,
with `r2` rounded from `1.234567891234567` to `1.2346`. This is a diagnostic
of the installed 3.2.1 build, not a claim about every future Open Babel version.

## Reproduction

```sh
cargo test -p chematic-mol --lib --test gzmat_precision --offline --quiet
cargo clippy -p chematic-mol -p chematic-py --all-targets --offline --no-deps -- -D warnings
cargo fmt --all -- --check
python -m pytest crates/chematic-py/tests/test_gzmat_precision.py -q
python scripts/validation/validate_gzmat_precision.py
python scripts/check_refactor_budget.py
python scripts/check_binding_surface.py
python scripts/check_product_names.py
mkdocs build --strict
```

The Python checks require the extension built from this source plus the pinned
RDKit and NumPy environment. The validator invokes the installed `obabel`
binary; `--obabel PATH` can select another 3.2.1 binary and records its hash.
Static binding inventory was updated for the two new Python functions. The
Rust integration suite is explicitly included in CI, with no new refactoring
budget allowances.

The candidate preserves the prior isolated source through RDKit meso-group
reduction. Its review patch is
`validation/openbabel_issues/gzmat-precision-2026-10-15.patch`, and evidence is
`validation/results/openbabel-gzmat-precision-2026-10-15.json`. Applying the
patch requires that inherited source context; it is not a substitute for
integrating and validating the earlier candidate changes. Main integration,
commit, push and release are separate steps.
