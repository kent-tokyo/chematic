# PDBQT coordinates and partial charges

The Rust PDBQT parser already returns coordinates and partial charges. The
Python `from_pdbqt` and WASM `mol_from_pdbqt` readers previously exposed only
atom elements. Python `convert_format` also discarded the arrays: PDBQT input
needed caller-supplied coordinates for coordinate output, and rewriting with
those coordinates defaulted charges to zero.

The additive readers expose numeric data in source ATOM/HETATM order:

```python
import chematic

mol, coords, partial_charges = chematic.from_pdbqt_with_data(source)
# coords: [[x, y, z], ...] in Å; partial_charges: [q, ...] in e
rewritten = chematic.convert_format(source, "pdbqt", "pdbqt")
xyz = chematic.convert_format(source, "pdbqt", "xyz")
```

```javascript
const { atomic_numbers, coords, partial_charges } =
  JSON.parse(wasm.pdbqt_atom_data_json(source));
```

These readers retain parsed finite binary64 values, including signed zero and
charges wider than the usual fixed field when separated from the type by
whitespace. They neither round nor recalculate those values. Partial charges
stay separate from formal charges; no formal charge is inferred from them.
Python accepts `str` or UTF-8 `bytes`; invalid UTF-8 raises `UnicodeDecodeError`.
Malformed atom records fail the entire read. An atom-free input returns empty
arrays, matching the existing Rust reader. WASM throws on invalid atom data and
bounds input to 1,000,000 bytes and atom records to 10,000 before building the
output. Native parser line budgets and the WASM 16 MiB output cap also apply.

The Python bridge now takes both arrays from PDBQT input unless individually
overridden. Its PDBQT output uses the [checked writer](pdbqt-checked-output.md):
exact per-atom counts, finite values and representable field widths are required.
Empty molecules, unsupported default types and overflowing overrides fail.
Coordinates and charges are written with three decimals, so the array reader
is needed when additional charge precision matters. Other input formats keep
their existing PDBQT writer behavior. PDB/XYZ and other bridge target formats
do not retain PDBQT partial charges.

## Scope

The [AutoDock format specification](https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf)
provides the coordinate/charge field layout. Widened charge tokens are an
explicit chematic reader extension, rather than a fixed-column compatibility
claim. The new APIs add access to existing Rust parsing; they do not expand its
supported atom-type vocabulary.

No bonds are perceived. The returned molecule is an element-only representation
with zero formal charges; it is not suitable evidence of original bond order,
implicit hydrogen count, aromaticity, isotope or stereochemistry. AutoDock
aliases are mapped to elements. Source atom serial/name/residue, source AutoDock
types, remarks and torsion trees are not retained. In particular a PDBQT rewrite
uses synthetic identities and heuristic types, and does not preserve ligand
flexibility or establish docking readiness. Retain the source separately when
those fields matter. Existing `from_pdbqt` and `mol_from_pdbqt` stay available.

## Reproduction

The six original MIT fixtures in `validation/pdbqt_atom_data/cases.json` specify
indexed element, coordinate and charge expectations, units, tolerances,
provenance and normalized SHA-256 hashes. They cover type aliases, explicit H,
charged metal/halide records, extended charge precision, signed zero, coordinate
width boundaries, CRLF and a missing final newline. These are development
fixtures used during implementation, not an independent holdout.

The baseline Python artifact is attested to main
`a83ed2653037a7b17033e03ee22501a1ddda8bb9` in the validation index. Observed
outcomes, including errors and complete probe stdout/stderr, are in
`validation/results/pdbqt-atom-data-{baseline,candidate}-2026-10-22.json`.

| Contract on six inputs | Baseline | Candidate |
| --- | --- | --- |
| Public Python numeric reader | Unavailable | Exact arrays 6/6 |
| PDBQT rewrite without coordinate overrides | 0/6, coordinates required | 6/6, charges within 0.000500000001 e |
| Rewrite with supplied coordinates | Charges defaulted to zero | Source charges retained to output precision |
| Node/WASM numeric reader | Unavailable | Exact arrays 6/6 |

A separate C++ probe linked to Open Babel **3.2.1** reads stored charges with
automatic bonding disabled (`b` input option); it does not run ligand preparation.
Both readers reproduce all indexed elements/coordinates/charges on the five
fixed-column fixtures. Open Babel refuses the widened-charge extension fixture;
Python and WASM return its source values. This comparison does not establish
superiority on standard PDBQT, other formats or docking chemistry.

```sh
# First build a Python extension and a Node/WASM package from this source.
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_pdbqt_atom_data.py \
  --output /tmp/pdbqt-atom-data.json
CHEMATIC_WASM_MODULE=/path/to/chematic_wasm.js \
  node crates/chematic-wasm/tests/pdbqt_atom_data.test.mjs

# Optional external read comparison; needs Open Babel development headers/libs.
c++ -std=c++17 $(pkg-config --cflags openbabel-3) \
  validation/openbabel_issues/pdbqt-checked/openbabel_read_probe.cpp \
  $(pkg-config --libs openbabel-3) -o /tmp/pdbqt-ob-probe
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_pdbqt_atom_data.py \
  --openbabel-probe /tmp/pdbqt-ob-probe --output /tmp/pdbqt-atom-data-paired.json
```

The Python regression suite also checks coordinate-only and charge-only
overrides, malformed numeric/type records, invalid UTF-8, count mismatch and
overflow. Real Node/WASM tests exercise the 10,000/10,001 atom boundary and the
input byte limit. The validation index records local artifacts, logs and source
hashes; passing source checks does not establish a published-package release.
