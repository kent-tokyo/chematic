# Checked CIF cell and coordinate output

Open Babel [issue #2818](https://github.com/openbabel/openbabel/issues/2818)
is open as checked on 2026-10-14. Its POSCAR fixture has a zero second cell
vector; the report shows NaN cell parameters and fractional positions in the
resulting CIF. This implementation addresses the same invalid-cell failure
class at chematic's plain CIF boundary. It does not add a POSCAR reader.

## APIs

Rust exports `chematic_mol::write_cif_checked(&mol, &coords, cell)` returning
`Result<String, CifError>`, and `UnitCell::validate()`.
Python exposes the same checked operation:

```python
import chematic

mol = chematic.from_smiles("[Si]")
text = chematic.write_cif_checked(
    mol,
    [[0.5, 1.5, 2.5]],  # Cartesian Å, exactly one triple per atom
    {"a": 4.0, "b": 6.0, "c": 8.0,
     "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
)
result = chematic.parse_cif(text)
```

All six cell keys are required in Python. Omitting `cell` writes Cartesian
`_atom_site_Cartn_x/y/z` columns; supplying a cell converts to fractional
`_atom_site_fract_x/y/z`. Output is minimal P1 atomic positions, without bonds,
occupancies, symmetry operations or other molecular metadata.

For explicit molecular connectivity, the separate
[`write_cif_with_bonds_checked` API](cif-molecular-bonds.md) adds a checked
bond-geometry loop. The existing writer and reader contracts above remain
unchanged.

The checked writer rejects:

- Non-finite or non-positive cell lengths, and non-finite angles or angles
  outside the open interval `(0, 180)` degrees.
- Singular, impossible or numerically near-singular cell geometry: the
  dimensionless Gram determinant
  `1 − cos²α − cos²β − cos²γ + 2 cosα cosβ cosγ` must exceed `1e-12`.
  This is a numerical policy; it can reject extremely flat but mathematically
  nonzero cells. The computed volume must also be finite and positive.
- Missing or excess coordinate triples (`InvalidCoordinate`), non-finite
  Cartesian coordinates, and non-finite Cartesian-to-fractional results
  (`InvalidCoordinate`). Invalid cells return `InvalidCellParameters`.

Python maps these errors to `ValueError`. Conversion errors return no CIF text.
Cell parameters and positions use round-trip float formatting, preventing the
legacy four-decimal cell formatting from rounding a small positive length to
zero. Floating-point coordinate transforms still have ordinary numerical error.

## Compatibility

The existing infallible Rust `write_cif` keeps its signature, decimal formatting
and coordinate-padding behavior. It remains unchecked and can emit NaN for an
invalid cell. Callers must select `write_cif_checked` to obtain validation.

The fractional-coordinate plain reader now shares `UnitCell::validate()`.
Consequently it also rejects negative lengths, out-of-range angles and cells
below the Gram-determinant threshold, including cases that previously passed
its volume-only check. Its existing partial-cell defaulting behavior is
unchanged. The periodic CIF and mmCIF adapters have separate validation and
serialization contracts; this change does not certify those paths.

## Local validation

Validated the uncommitted source working tree on 2026-10-14:

- `chematic-mol`: 998 tests passed, zero failures or ignored tests, including
  eight new checked-writer regressions and the five earlier CIF coordinate tests.
- Python: 45 new checked-writer cases and 19 existing crystal cases passed using
  the locally compiled source extension. Cases cover every cell field,
  singular/impossible/near-singular angles, orthogonal/triclinic/small cells,
  missing fields, coordinate counts, non-finite coordinates and conversion overflow.
- Clippy for `chematic-mol` and `chematic-py`, all targets, passed with warnings
  denied. Changed Rust files passed rustfmt, and `git diff --check` passed.

The reproducible probe script is `scripts/validate_openbabel_cif_checked.py`;
its JSON evidence is `validation/results/openbabel-cif-checked-2026-10-14.json`.
It records fixture, source and loaded-extension hashes, 36 invalid-cell
rejections, four local round trips and bounded external probes.

The installed Open Babel reports version 3.2.1. Both the singular POSCAR
conversion and a normal CIF-to-XYZ control exceeded the 10-second timeout;
the XYZ-to-CIF control succeeded. Therefore the external CIF comparison is
inconclusive on this host: these timeouts do not establish an issue-specific
failure or a successful cross-read. The issue's reported NaN behavior and
chematic's checked rejection are distinct evidence. The POSCAR reproduction
is retained in `validation/openbabel_issues/poscar/singular-cell.POSCAR`.

This is local source validation, not a package release or an upstream fix.
