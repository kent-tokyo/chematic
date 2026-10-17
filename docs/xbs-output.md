# Checked XBS scene output

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
view commands, periodic cell metadata, movies, reading, GUI execution and
chemical graph conversion are outside this profile. It has no Rust/WASM
binding and does not add a format alias to `convert_format`.

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
