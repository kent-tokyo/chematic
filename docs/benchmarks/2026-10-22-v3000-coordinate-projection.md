# V3000 coordinates through CJSON projection

Python `convert_format(source, "v3000", "cjson")` and WASM
`convert_common_format(source, "v3000", "cjson")` previously read only the graph
and wrote CJSON without coordinates. They now carry the source atom-row x/y/z
values into the existing CJSON writer. This is a cross-format projection rather
than immutable source copying.

```python
mol, name, xyz = chematic.from_mol_v3000_with_coords_3d(source)
output = chematic.convert_format(source, "v3000", "cjson")
```

Rust exports `parse_mol_v3000_with_coords_3d`, returning
`(Molecule, MolMetadata, Vec<(f64, f64, f64)>)`. The new coordinate API is bounded
to 16 MiB input and retains the native atom/bond/count limits. WASM common
conversion additionally enforces its existing 1,000,000-byte/10,000-atom caps.
The parser's source z buffer is returned alongside the existing report; legacy
readers discard that buffer without copying it. Existing 2D API shapes and
stereo-perception behavior stay unchanged.

Source coordinates are finite binary64 values in Angstrom. They retain signed
zero, values below perception's flatness tolerance and planes offset from z=0.
Dimension headers still affect stereo interpretation; they do not zero source
coordinate components. A missing perception conformer does not imply missing
source coordinates. Decimal spelling, original dimension declarations and
metadata are not preserved in CJSON.

Python common conversion also supplies these coordinates to PDB, XYZ and PDBQT
writers unless the caller passes `coords`, including an explicit empty list.
Those writers keep their existing precision and metadata limitations. No
partial charge is inferred from a V3000 formal charge. Other common-bridge
outputs retain their existing behavior.

This does not verify cross-format stereochemistry, implicit hydrogen semantics,
query/polymer meaning or opaque attributes. The coordinate API inherits the
native V3000 parser profile; it is not a general CTfile validator. The new
path does not replace source retention or checked conversion reports.

## Development packet and observations

Eight original MIT synthetic fixtures cover tiny z, high precision, offset
planes with formal charges, sparse source IDs and carbon-13, flat signed zero,
attribute-tail continuation, a 2D header with nonzero z and four-atom noncoplanar
geometry. LF, CRLF and no-final-newline give **24** inputs. This is a development
packet used during implementation, not an independent holdout or an upstream
issue reproduction.

The oracle declares atomic numbers, source-row coordinates as decimal strings,
formal charges, isotope masses and bonds between zero-based row indices. The
independent JSON observer checks indexed elements/charges/bond endpoints/orders
separately from coordinate presence, exact binary64 components and the existing
local `atoms.isotopes` extension. Absent empty bond sections normalize to zero
bonds; that formatting choice is not a loss. This is a narrow indexed-graph
check, not graph isomorphism, atom-order-independent equivalence or stereo.

| Route | JSON object | Indexed graph | Coordinates present | Exact coordinate bits | Local isotope array |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baseline Python common | 24/24 | 24/24 | 0/24 | 0/24 | 24/24 |
| Baseline WASM common | 24/24 | 24/24 | 0/24 | 0/24 | 24/24 |
| Candidate Python common | 24/24 | 24/24 | 24/24 | 24/24 | 24/24 |
| Candidate WASM common | 24/24 | 24/24 | 24/24 | 24/24 | 24/24 |
| Open Babel 3.2.1 `-imol -ocjson` | 15/24 | 15/24 | 15/24 | 15/24 | 12/24 |
| Open Babel 3.2.1 `-imol -ocjson -xs` | 15/24 | 15/24 | 15/24 | 15/24 | 12/24 |

Open Babel retains every measured coordinate bit in all 15 emitted objects;
there is no rounding loss among those objects. The other nine have empty stdout
with exit code zero: three flat-2D variants, three 2D-header/nonzero-z variants
and three attribute-tail-continuation variants. The two 2D families report no
3D coordinates; the continuation family reports zero molecules. Three emitted
sparse-ID/carbon-13 variants omit `atoms.isotopes`. The local isotope array is
not declared in the linked CJSON schema; this is an extension-level observation,
not a schema-compliance advantage. Full output is retained for inspection.

The continuation fixture uses `-` after the atom-map token, then `CHG=0` on the
next V30 line. BIOVIA's atom-block example places optional attributes on a
continued line. A diagnostic unfolding of this fixture succeeds in Open Babel;
its terminal result is recorded separately from the scored 24 inputs. No
coordinate generation, hydrogen addition, normalization or raw-copy mode is
used. `-xs` selects single-geometry CJSON output, not a precision override.

Baseline main package source is `b9412a06` with its attested native extension;
WASM is the attested `c713b21a` artifact whose unrelated CJSON same-format route
does not change the V3000-to-CJSON branch. Relevant native/Python sources are
identical across those bases. Candidate source starts at main `d4208b9d`; its
intervening CIP guard does not change this V3000 path. The previous immutable
source-retention PR is not a dependency.

The result establishes better preservation and usable conversion on this small
local workflow. It does not establish toolkit-wide superiority, broad format
coverage, performance, release availability or complete CTfile/CJSON semantics.

## Reproduction and evidence

With built Python and Node-target WASM packages:

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_v3000_coordinate_projection.py \
  --wasm-module /path/to/chematic_wasm.js --obabel /path/to/obabel \
  --output /tmp/v3000-coordinate-projection.json
```

`--baseline` records prior outputs without requiring candidate success.
The subprocess runner bounds both pipes and elapsed time. Native tests check
source coordinates, legacy 2D coordinates, graph fields and CJSON reparse;
Python and real Node/WASM tests exercise the common routes, malformed input and
bounds. Python tests also cover explicit coordinate precedence and downstream
coordinate writers. Observer negative controls cover precision, signed zero,
z, element order, charges, bonds, isotopes, missing coordinates, duplicate keys
and nonfinite JSON. CI discovers the native/Python/Node suites and runs observer
unit tests.

The packet is `validation/v3000_coordinate_projection/cases.json`. Compressed
complete baseline/candidate results, terminal logs, artifact attestation and
source hashes are indexed in
`validation/results/v3000-coordinate-projection-validation-2026-10-22.json`.

References: [BIOVIA CTfile Formats 2020, general syntax and atom block](https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf),
[Open Chemistry CJSON schema, Angstrom coordinate arrays](https://github.com/OpenChemistry/chemicaljson/blob/main/cjson.schema).
