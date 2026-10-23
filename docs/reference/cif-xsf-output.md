# Checked CIF-to-XSF geometry output

Python `write_xsf_checked` exports an existing `PeriodicStructure` as one
static three-dimensional periodic XSF structure. This supplies an output
path for [Open Babel #2252](https://github.com/openbabel/openbabel/issues/2252):
CIF geometry for XCrySDen. It reuses the native CIF reader and validated
lattice/site types; it adds no crystal parser or symmetry-search algorithm.

```python
from pathlib import Path
from chematic import PeriodicStructure, write_xsf_checked

structure = PeriodicStructure.from_cif(Path("input.cif").read_text())
text = write_xsf_checked(structure, allow_metadata_loss=True)
Path("output.xsf").write_text(text, encoding="ascii")
```

The output contains `CRYSTAL`, `PRIMVEC` and `PRIMCOORD`, atomic numbers and
Cartesian coordinates. Lattice vectors are rows and all lengths are in
angstroms. Directly constructed structures must use the same unit convention.
The supplied basis, atom order and Cartesian positions are retained with
17 significant digits. The writer does not wrap coordinates, reduce a cell,
find symmetry, infer topology or claim a crystallographically primitive basis:
XSF's primitive keywords carry the supplied periodic cell. No conventional
cell basis is generated.

## Geometry and metadata boundaries

`allow_metadata_loss=True` is required because this geometry profile omits
site labels, CIF source/provenance and any POSCAR constraints, velocities or
restart data. The original structure is not modified. Keep source documents
if those fields are needed; this is not a full CIF/POSCAR metadata converter.

Even with that option enabled, the writer always refuses:

- A CIF-sourced structure whose `symmetry_status.is_complete_cell` is false.
  The same guard holds after wrapping/supercell operations preserve that status.
- A site containing multiple species or an occupancy other than exactly 1.0,
  including a zero-occupancy extra species. No disorder is collapsed.
- Empty structures, more than **50,000 sites**, or nonfinite lattice/Cartesian
  output. Finite fractional coordinates can overflow during Cartesian
  conversion; the writer detects this before returning text.

`PeriodicStructure.from_cif(..., expand_symmetry=True)` already expands
operations explicitly present in the source. Its complete-cell flag reflects
P1 or that expansion, not verification against a space-group database.
Space-group names alone are not converted to operations by this writer.
Directly constructed structures have no CIF provenance; their listed sites
are treated as the caller's declared cell.

Wrong structure/option types raise `TypeError`; omitted loss consent,
unsupported geometry and range/numeric failures raise `ValueError`. The site
cap is checked before native site/Cartesian-array getters. The helper uses
existing native getters and their NumPy arrays; NumPy is required for this
workflow. It returns text and performs no filesystem writes.

Only 3D periodic geometry is included. Molecular/slab/polymer modes, forces,
datagrids, animations and conventional-cell inference remain outside this
profile. Generic molecule converters are unchanged.

## Source and local evidence — October 17, 2026

Implemented independently from the
[XCrySDen XSF specification in its distributed documentation](https://web.mit.edu/xcrysden_v1.5.60/www/XCRYSDEN/doc/XSF.html).
The main XCrySDen documentation host timed out; the distributed MIT mirror
was fetched and hashed. No external implementation/example coordinates were
copied. Original synthetic CIF input is in
`validation/openbabel_issues/xsf/ordered-triclinic.cif`.

Baseline `35b932f9` has no XSF writer. Local Open Babel **3.2.1** lists XSF
as **read-only**, and a CIF-to-XSF attempt returns `Invalid output format`.
The candidate passed **198 tests**: 17 new XSF tests, 19 existing crystal tests
and 162 existing compatibility tests. The tests cover basis orientation,
unwrapped positions, input preservation, occupancy/symmetry refusal, POSCAR
extras retention on the source, overflow and the 50,000-site bound.

**ASE 3.26.0** read 60 synthetic cases / 1,360 atom rows with lattice and
Cartesian coordinates matching native getters as binary64 bytes. Cases
include cubic, tilted and left-handed bases, wrapping and supercells.
Additional checks cover the CIF pipeline and all 118 element symbols.
**Open Babel 3.2.1** read 48 cases through its XSF reader; element/order/count
checks passed and XYZ coordinate output differed by at most
`4.984001341412636e-6 A`, within its five-decimal output rounding. Cell-vector
parity was checked with ASE, not with that coordinate-only XYZ output.
XCrySDen GUI execution and rendering have not been tested.

The public type stub also now describes the existing `expand_symmetry` flag
and four read-only symmetry-status properties used by this workflow; runtime
crystal behavior is unchanged. Targeted strict helper mypy passes with
imported-stub diagnostics silenced, valid public usage passes, and four
invalid calls/assignments give four expected errors. A direct strict run
also exposed 81 existing generic-annotation diagnostics in the legacy stub;
whole-package strict cleanliness is not claimed.

Source/stub/binding/unsafe/product-name gates and strict MkDocs passed.
Commands, hashes and logs are in `validation/results/cif-xsf-2026-10-17/`:

```sh
python3.13 validation/validate_cif_xsf.py --output /tmp/cif-xsf-crossread.json
```

The existing native extension was reused for structure operations and tests;
no new native binary was built. This evidence does not establish a package
release, metadata round trip, rendering result or performance advantage.
