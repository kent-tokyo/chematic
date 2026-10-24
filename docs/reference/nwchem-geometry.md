# NWChem Cartesian GEOMETRY blocks

The direct NWChem APIs read, edit, and write **one molecular Cartesian
`GEOMETRY` block**, rather than a complete calculation deck. They preserve the
represented center fields and header controls without perceiving bonds,
adding atoms, changing tags, centering coordinates, or converting units.

## Represented meaning

Each center retains its case-sensitive tag and three coordinates. Repeated
tags are allowed: tags can associate several centers with the same basis or
ECP settings. Element-prefixed tags, `bq` dummy-center tags, and `x` dummy-center
tags are accepted; `Xe` tags remain unmodified. The API stores tags directly
and does not project them into an element-only molecular graph.

Optional center fields remain distinct from absent values:

| Field | Meaning |
|---|---|
| `velocity_au` | Three velocities in atomic units, independent of coordinate units |
| `nuclear_charge` | Explicit charge of the nucleus/center, **not** a molecular formal charge or atomic partial charge |
| `mass_amu` | Explicit positive mass in atomic mass units |
| `nucleus` | `point` or `finite` nuclear model; center override is separate from the geometry default |

Header fields retain a geometry name, explicit coordinate-unit alias, custom
angstrom-to-AU factor, center/nocenter, autosym/noautosym and its optional
tolerance, autoz/noautoz, print/print xyz/noprint, bqbq, and default nuclear
model. Optional controls stay absent when omitted. This preserves the runtime
defaults rather than filling them in with library guesses.

The accepted unit spellings are `angstroms`, `an`, `au`, `atomic`, `bohr`, `nm`,
`nanometers`, `pm`, and `picometers`. Their spelling is retained in lowercase;
an omitted unit remains `null` and means the documented NWChem angstrom
default. Coordinates retain the declared unit. No scaling is performed, even
when a custom conversion factor is present.

Numbers are finite f64 values. Output uses seventeen decimal places in
scientific notation to retain their values, including subnormal values.
Fortran D exponents are accepted. Overflow, nonfinite values and nonzero
numbers that underflow to zero refuse. Mass, conversion factor and explicit
symmetry tolerance must be positive.

## APIs

```python
import chematic

record = chematic.parse_nwchem_geometry("""geometry water units bohr nocenter
O_ether 0.12345678901234566 0 0 charge 7.5 mass 17.999159612
H1 0 0 1.8
end
""")
record["centers"][0]["mass_amu"] = 18.01010101010101
text = chematic.write_nwchem_geometry(record)
assert chematic.parse_nwchem_geometry(text) == record
```

Rust exposes `NwchemGeometry`, `NwchemCenter`, `NwchemLimits`,
`parse_nwchem_geometry`, `parse_nwchem_geometry_with_limits`,
`write_nwchem_geometry`, and `nwchem_geometry_from_json`. Python accepts text
or bytes and returns a dictionary. WASM/Node exposes `nwchem_geometry_json`
and `nwchem_geometry_from_json`. JSON/dictionary writers refuse unknown fields
at all levels and validate the current edited record.

```sh
chematic nwchem-geometry -i geometry.nw -o normalized.nw
chematic nwchem-geometry -i geometry.nw --inspect -o geometry.json
chematic nwchem-geometry -i geometry.json --from-json -o edited.nw
```

Input comments retain their ordered text in `comments`; canonical output moves
them before the block. Inline comment positions, whitespace, keyword order,
keyword case and nuclear-model aliases are normalized. This is an editable
record API, not a byte-identical source route.

## Boundaries

Full input decks containing START/RESTART, TITLE, CHARGE, BASIS, SCF/DFT, TASK,
or any other outer directive refuse. This prevents a geometry-only extraction
from silently dropping calculation settings. In particular, the top-level
CHARGE directive is separate from a center's nuclear charge and is unsupported
here. NWChem output files (`nwo`) are a different format and remain unsupported.

This initial profile also refuses periodic lattices, explicit SYMMETRY groups,
Z-matrices/internal-coordinate edits, ADJUST, LOAD, INCLUDE/EOF, multiple
geometry blocks, quoted/escaped tokens, continuation lines, semicolon-separated
statements, unknown attributes, and duplicate/conflicting options. Names/tags
use unquoted ASCII alphanumerics, underscore, dot or hyphen. Tags are limited to
16 bytes and must have a recognized element or dummy prefix; reserved
directives cannot be tags. Geometry names are limited to 64 bytes and cannot
collide with geometry options.

Default limits are 16 MiB input/output, 1023 bytes per physical input line, and
100000 centers. JSON input is limited to 64 MiB. WASM further limits text input
to 1000000 bytes, centers to 10000, and JSON/output to 16 MiB. The writer uses
the current structured record, rather than replaying stale input after edits.

## Validation and comparison

The original packet in `validation/nwchem_geometry/cases.json` contains 25
valid block cases and 57 refusal cases, shared by Rust, Python, WASM and CLI.
Tests cover omitted/explicit controls, all nine unit aliases, custom factors,
dummy centers, repeated/case-sensitive tags, center overrides, isotope masses,
velocities, comments and numerical limits. Native boundary tests retain the
bits of the smallest positive subnormal, smallest normal, maximum finite f64,
and negative zero. WASM tests exercise the 10000/10001-center boundary.

Open Babel 3.2.1 advertises `nw` as **write-only**, so it cannot serve as a
same-format round-trip comparator. The committed comparison instead measures
**six original indexed XYZ geometries exported to NWChem Cartesian input**.
The chematic JSON adapter is independently checked against the same XYZ tags
and angstrom coordinates before execution. Open Babel uses
`-ixyz -onw -ab -l1`. Both tools produced observable geometry blocks in 6/6
cases; exact indexed coordinates survived in chematic 6/6 and Open Babel 1/6.
The remaining Open Babel outputs rounded coordinates to five decimal places.
This is a small precision diagnostic, not a representative chemistry benchmark.

Only the tags and Cartesian coordinates/units are compared. The observer does
not validate the surrounding Open Babel job header, and attributes absent
from XYZ are not scored. It does not claim that Open Babel lost nuclear
charges, masses or velocities that were never supplied in that input.
Full commands, per-case output/error receipts and source/input/binary hashes
are in `validation/nwchem_geometry/openbabel-3.2.1-comparison.json`.

```sh
python3 scripts/validation/check_nwchem_geometry.py \
  --cli target/debug/chematic --obabel /path/to/obabel \
  --output /tmp/nwchem-comparison.json
```

The implementation follows the public manual. Local validation builds the
real Python and WASM bindings and exercises parsing/output; it does not execute
NWChem calculations. No NWChem executable was available for this local run.

## Primary references

- [NWChem geometry keywords and units](https://nwchemgit.github.io/Keywords-for-the-GEOMETRY-directive.html)
- [NWChem Cartesian centers and optional attributes](https://nwchemgit.github.io/Cartesian-coordinate-input.html)
- [NWChem input syntax](https://nwchemgit.github.io/Getting-Started.html)
- [Open Babel NWChem input format and write-only direction](https://openbabel.org/docs/FileFormats/NWChem_input_format.html)
