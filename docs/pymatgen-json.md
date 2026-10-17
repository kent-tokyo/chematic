# pymatgen Structure JSON input

[Open Babel #2675](https://github.com/openbabel/openbabel/issues/2675), verified
open on 2026-10-17, requests direct reading of pymatgen/Materials Project JSON.
`chematic.parse_pymatgen_structures` adds this input to Python without installing
pymatgen or Monty. It accepts a bare serialized `Structure`/`IStructure`, or a
Materials Project response shaped as `{"data": [{"structure": ...}, ...]}`.
It reads supplied text; it does not contact Materials Project or require an API
key. A C++ plugin and a Rust reader are outside this change.

```python
from chematic import parse_pymatgen_structures

text = '''{
  "@module": "pymatgen.core.structure",
  "@class": "Structure",
  "lattice": {"matrix": [[4,0,0], [0,5,0], [0,0,6]],
              "pbc": [true,true,true]},
  "charge": 0,
  "sites": [{"species": [{"element":"Fe", "occu":0.6},
                         {"element":"Ni", "occu":0.3}],
             "abc": [-0.25, 0.5, 1.25], "label":"metal"}]
}'''
record = parse_pymatgen_structures(text)[0]
assert record.source_text == text
assert record.projection_losses == ("structure.charge",)
# Native geometry has no structure-charge field, including an explicit zero.
geometry = record.to_periodic_structure(allow_metadata_loss=True)
assert geometry.sites[0].species == [("Fe", 0.6), ("Ni", 0.3)]
assert record.structure_dict()["charge"] == 0
```

## Meaning retained

The matrix rows are lattice vectors in Å; `abc` is fractional, and `xyz` is
Cartesian in Å. Coordinates stay outside the unit cell when supplied that way.
There is no wrapping, inferred symmetry, site expansion, bond perception or
merging of coincident sites. The input is treated as the explicit periodic cell.
The [pymatgen serialization implementation](https://github.com/materialsproject/pymatgen/blob/v2025.6.14/src/pymatgen/core/structure.py)
and [API documentation](https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure.Structure.from_dict)
are the reference for the input layout.

Each immutable `PymatgenStructureRecord` contains:

- `source_text`: the complete original decoded UTF-8 document, retaining JSON
  spelling, whitespace, response metadata and all records. Bytes must decode
  as UTF-8; this is a text record rather than a byte-for-byte container.
- `record_index`: the zero-based `data` entry index, or zero for a bare structure.
- `structure_dict()`: a fresh ordinary JSON-value dictionary for this structure,
  including charge, oxidation states/spin, properties and unknown fields.
  Dictionary values follow JSON/Python numeric decoding; original numeric
  spelling remains in `source_text`. Serialized classes are never imported or
  reconstructed through Monty.
- `projection_losses`: paths of information absent from the native geometry.
  `response` summarizes response-envelope metadata; `data[i].material_id`, for
  example, identifies per-entry metadata. The entire response remains in source.

The native `PeriodicStructure` projection keeps the lattice matrix, site/species
order, labels, fractional coordinates and all element/occupancy entries. A mixed
site stays mixed, vacancies remain partial occupancies, and separate Fe2+/Fe3+
entries remain separate Fe occupancy entries. Oxidation states, spin, charge,
site/structure properties and unknown metadata remain in the record; they are
not native geometry attributes. By default `to_periodic_structure()` refuses
when `projection_losses` is nonempty. Opting in does not remove the source or
metadata from the record. Common `Structure.as_dict()` output includes an
explicit `charge: 0.0`, which also requires this opt-in.

## Checked input profile

Only full 3D PBC is supported. Missing `pbc` means `[true,true,true]`; partial
PBC, nonboolean flags, molecules, derived structure classes such as slabs and
dummy elements raise. Optional class/module tags must identify `Structure` or
`IStructure` in `pymatgen.core.structure` when non-null. Unknown metadata is
retained, but it cannot override native geometric validation.

`abc` and a 3×3 lattice matrix are required. If `xyz` is also present, it must
agree with `abc @ matrix` within absolute 1e-7 Å and relative 1e-8 tolerance.
Optional redundant `a/b/c`, `alpha/beta/gamma` and `volume` must agree with the
matrix within absolute/relative 1e-7 tolerance (lengths in Å, angles in degrees,
volume in Å³). Inconsistent dual representations raise instead of selecting a
silent alternative. The matrix is authoritative within these rounding limits.

Native lattice and occupancy validation applies: singular/poorly conditioned
cells, invalid element symbols, negative occupancies and overfull sites raise.
Coordinates and geometric numeric fields exclude booleans and require finite
magnitudes ≤1e100. JSON duplicate keys, NaN/Infinity and float64 overflow or
nonzero-to-zero underflow are rejected, including in otherwise opaque metadata.
Malformed responses and a bad later record fail the entire call; no record is
silently skipped or returned as an unrelated first structure.

Limits are checked before constructing native objects: 16 MiB UTF-8 input,
64 JSON nesting levels, 1..256 records, 1..50000 sites per structure and 50000
sites in total, 1..32 species entries per site and 100000 species in total.
The text parser itself first checks the input size/depth; collection sizes are
checked after JSON decoding. Record dictionaries are separate copies, and
record attributes are frozen. This is not a lazy file/stream reader.

## Validation

The 2026-10-17 regression run passed 258 tests: 77 new input/source/metadata/loss
cases, 19 existing crystal cases and 162 compatibility cases. The fixed-version
pymatgen 2025.6.14 oracle serialized 384 synthetic structures with 3936 sites and
5448 species entries, including triclinic/left-handed cells, out-of-cell
coordinates, disorder, oxidation states and site properties. It also checked
384 two-record response envelopes and one additional 118-element structure.

Matrices, fractional coordinates, labels and element/occupancy entries matched
exactly in that corpus. The largest Cartesian component difference was
7.11e-15 Å (absolute tolerance 1e-12 Å). Re-reading the record's dictionary with
pymatgen retained its original `as_dict()` content, including oxidation states
and site properties. These are generated compatibility cases, not a live
Materials Project request or an assertion about all pymatgen subclasses.

On this host Open Babel 3.2.1 did not recognize the `pymatgen` input format.
That establishes a missing input route, not a numerical defect in another
Open Babel reader. Reproduce the oracle with
`python validation/validate_pymatgen_json.py` in an environment with pymatgen
2025.6.14. Evidence lives in `validation/results/pymatgen-json-2026-10-17/`.

The parser is Python-only, using existing native crystal constructors. An
existing extension was reused as an import/constructor fixture; this run does
not establish a fresh native build, wheel publication or performance result.
