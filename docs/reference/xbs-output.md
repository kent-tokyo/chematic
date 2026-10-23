# Checked XBS scene input and output

Python `write_xbs` emits the static `atom`, `spec` and `bonds` records accepted
by the XBS ball-and-stick scene format. It is inspired by
[Open Babel #2576](https://github.com/openbabel/openbabel/issues/2576), which
asks for XBS-style output. Open Babel's `bs` is a **different Ball and Stick
format**, despite sharing the `.bs` extension. The new API explicitly names
XBS and does not redefine `bs` in any generic converter.

```python
from pathlib import Path
from chematic import XbsSpecies, XbsBondRule, write_xbs

text = write_xbs(
    ["O", "H", "H"],
    [(0., 0., 0.), (0.92, 0.13, 0.), (-0.2, 0.88, 0.07)],
    species=[XbsSpecies("O", 0.8, 0.4), XbsSpecies("H", 0.5, 1.0)],
    bonds=[XbsBondRule("O", "H", 0., 1.3, 0.1, 0.7)],
)
Path("water.bs").write_text(text, encoding="ascii")
```

`XbsSpecies(label, radius, gray)` defines a species drawing radius and gray
value, from black 0 to white 1. `XbsBondRule(species1, species2, min_distance,
max_distance, radius, gray)` declares a species-pair distance rule for drawing.
These rules describe a scene; they do not encode chemical bond orders.
All coordinates, radii and distance limits share the caller's chosen unit,
normally angstroms. XBS records carry no unit declaration or conversion.

The writer preserves atom order and emits 17 significant digits, sufficient
to round-trip finite binary64 values through a binary64 reader. Native viewer
precision and rendering results are separate and have not been measured.
Species styles and rules are emitted in caller order. No palette, default
radius, bond rule or molecular topology is inferred. Extra unused styles
are permitted. The function returns a complete string and writes no file.

## Checked profile

- Python sequences such as lists and tuples are accepted; generators and
  NumPy arrays must first be materialized as ordinary lists/tuples.
- Between 1 and 10,000 atoms; exactly three coordinates per atom, finite in
  ±1e6. Labels and positions must have equal lengths.
- Up to 128 species styles and 1,024 bond rules. Size checks precede copying
  or indexing. The writer uses bounded indexed access to style/rule inputs.
- Labels are 1–8 ASCII letters/digits/underscores, beginning with a letter.
  They are scene species names, not validated chemical element symbols.
  Wildcards, whitespace and executable viewer commands are refused.
- Each used label needs exactly one style. Undefined rule species and
  duplicate undirected pair rules (including reversed pairs) are refused.
- Radii must be positive and at most 1e6. Gray values must be in [0, 1].
  Distance limits satisfy `0 <= min_distance < max_distance <= 1e6`.
  All numeric fields must be finite real values; booleans are refused.
- Style/rule dataclasses are frozen. Invalid types raise `TypeError`;
  invalid values, dimensions, references or bounds raise `ValueError`.
  Source sequence errors propagate. Inputs are not modified.

This is a standard-library Python writer. RGB/named colors, wildcard rules,
view commands, periodic cell metadata, movies, GUI execution and
chemical graph conversion are outside this profile. It has no Rust/WASM
binding and does not add a format alias to `convert_format`.

## Source-preserving input — October 19, 2026

`read_xbs(text)` and `XbsScene(text)` read the same static gray-style subset.
The author-supplied [XBS README](https://github.com/Montmorency/XBS#readme)
documents `atom`, `spec`, `bonds` and whole-line `*` comments. This reader extends
the existing output capability motivated by Open Babel #2576; it is not a fix
to Open Babel's separate `bs` dialect.

```python
from chematic import XbsScene, read_xbs

text = "* Original abstract site\r\natom Site1 -0 0 1.25\r\nspec Site1 0.5 0.75\r\n"
scene = read_xbs(text.encode("ascii"))
assert isinstance(scene, XbsScene)
assert scene.labels == ("Site1",)
assert scene.positions == ((-0.0, 0.0, 1.25),)
assert scene.atom_lines == (2,)
assert scene.to_xbs() == text
normalized = scene.to_xbs(preserve_source=False)
assert read_xbs(normalized).labels == scene.labels
```

The immutable result retains `source_text`, ordered `labels` and `positions`,
ordered `species` styles and `bonds` drawing rules, and their 1-based physical
`atom_lines`, `species_lines`, and `bond_lines`. Nested containers are tuples;
styles and rules are frozen. Construction itself validates input. The compact
representation reports counts instead of printing the complete source.

`to_xbs()` returns exact decoded source, preserving comments, whitespace, CRLF,
decimal spelling and record interleaving. `to_xbs(preserve_source=False)` uses
the existing writer to emit LF, grouped records and 17 significant digits;
this deliberately omits comments and original numeric spelling. Parsed binary64
values, including negative zero and subnormal coordinates, remain unchanged.
ASCII byte input can be reconstructed exactly from `source_text.encode("ascii")`.

Reader-specific limits and behavior:

- Input is `str` or immutable ASCII `bytes`, at most 4 MiB, 20,000 physical
  lines and 4,096 bytes per physical line. LF/CRLF and space/tab separation
  are accepted. Bare CR, other controls, non-ASCII and unsupported records fail.
- Atoms/styles/rules retain the writer's 10,000/128/1,024 limits, label syntax,
  numeric ranges, required styles and unique undirected pair rules. Unused styles
  are allowed. Records may interleave and refer to styles declared later.
- Each physical line contains one complete record, a blank line or a whole-line
  `*` comment. Inline comments, semicolon commands, named/RGB colors, wildcard
  rules, movies and viewer commands are refused instead of ignored.
- Numbers use decimal syntax with optional exponent. Nonfinite numbers and
  nonzero decimal literals that round to zero in binary64 are rejected. Other
  decimal values use normal binary64 rounding; exact decimal arithmetic is not
  promised. Coordinate, radius and distance units are retained without inference.
- Wrong input types raise `TypeError`; malformed/unsupported data raise
  `ValueError`. Record errors include the physical source line; byte/total-line
  caps are document errors. There is no partial scene return or filesystem I/O.

The reader is an explicit scene-data API. Species labels need not be chemical
elements. Distance-based drawing rules do not encode molecular connectivity,
bond orders, charges or periodicity; no graph or viewer is created.

Local follow-up validation passed **271 tests** (62 new reader cases, 47 existing
writer cases and 162 compatibility cases), including positive/negative resource
boundaries. Strict mypy passed for the two XBS helpers and a public consumer;
four invalid consumer calls produced four expected errors. The package root's
existing generic-annotation errors were excluded from this scoped check.

Independent cross-reading with **hofmann 0.21.1** covers 144 authored scenes /
8,016 atom rows in scientific/general decimal notation, both atom orders and
LF/CRLF. Original and normalized text both match binary64 coordinates, styles
and rules: 288 ordinary cross-reads plus one 10,000-atom scene. Original source
text and physical atom-line correspondence are checked independently too.
No GUI/rendering, bond perception, fresh native build or published wheel is
claimed. See `validation/results/xbs-input-2026-10-19/` and reproduce with:

```sh
python3.13 validation/validate_xbs_input.py --output /tmp/xbs-input-crossread.json
```

## Sources and local validation — October 17, 2026

Implemented independently from the author-supplied
[XBS manual distributed by Debian](https://manpages.debian.org/bookworm/xbs/xbs.1.en.html),
with an external parser check through
[hofmann's documented XBS reader](https://hofmann.readthedocs.io/en/latest/xbs-format.html).
No viewer implementation or manual example coordinates were copied.
The fixtures are original synthetic project data.

The exact Python files at baseline `35b932f9` have no `write_xbs`,
`XbsSpecies` or `XbsBondRule`. The candidate passed **209 tests**: 47 new
XBS regressions plus 162 existing compatibility tests. Strict helper mypy
passed; the public typed consumer passed, and four invalid consumer calls
produced the four expected errors. The legacy package was not checked in
strict mode as a whole.

A seeded cross-read through **hofmann 0.21.1** passed 96 scenes / 5,344 atom
rows, 288 style checks and 96 rule checks, plus a separate 10,000-atom scene.
Coordinates matched binary64 bytes, including negative zero and subnormal
values. No rendering, bond detection or original XBS GUI run is claimed.
Open Babel **3.2.1** produced its different `bs` dialect on the original
three-site fixture, and its format listing had no `xbs` identifier.
That is a format distinction, not a newly reproduced defect in its `bs` writer.

Evidence and commands are in `validation/results/xbs-output-2026-10-17/`.
Reproduce with the candidate package, hofmann 0.21.1 and NumPy available:

```sh
python3.13 validation/validate_xbs_output.py --output /tmp/xbs-crossread.json
```

The local native extension was reused for package import and existing tests;
this addition contains no native calls. Source gates and strict MkDocs
passed. No fresh native build, wheel/release publication, performance result
or upstream Open Babel change is implied.
