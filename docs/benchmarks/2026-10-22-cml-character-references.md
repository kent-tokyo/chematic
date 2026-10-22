# CML character references and typed-field preservation, 2026-10-22

CML XML character references previously changed the parsed chemistry. A formal
charge `&#49;` silently became zero, an isotope reference became unspecified,
and an encoded coordinate became zero. Encoded elements, IDs and bond orders
could fail to read. The shared attribute decoder also decoded `&amp;lt;` twice
into `<`, changing literal CDXML attribute values.

The candidate decodes the five predefined entities and valid decimal/hex XML
character references once, according to the reference and character rules in
[XML 1.0](https://www.w3.org/TR/xml/#sec-references). The decoder scans the source
in linear passes and does not rescan decoded ampersands. It does not expand
custom DTD entities or make the lightweight readers full XML validators.

CML readers also reject malformed or out-of-core-range formal charges,
isotopes and supplied 2D coordinates. Missing fields retain their existing
defaults. The `isotope` alias and `isotopeNumber` must agree numerically when
both exist, and atom IDs must be unique after reference decoding. These checks
apply to both normal and structurally strict readers. No public API signature
or error-enum variant is added.

## Same-input development packet

The packet uses an authored nitrogen-15 ammonium atom attached to one explicit
deuterium atom. Literal input and six equivalent lexical variants test charge,
isotope, element, coordinate, atom ID and bond order. Seven additional inputs
contain malformed fields, conflicting aliases, duplicate decoded IDs or an
invalid/literal reference in a numeric field. Coordinates are supplied 2D
Angstrom values, not optimized geometry.

Both engines receive identical input bytes on CML → CML. Python ElementTree
decodes the XML independently, then compares authored indexed elements,
charges, isotopes, total attached H counts, coordinates and explicit bonds.
Values use exact integers and dyadic coordinates with no tolerance. Each
engine has a terminal row for every input; stdout, stderr and successful
output text are retained.

| Terminal outcome, all 14 inputs retained | Before repair | Fixed source | Open Babel 3.2.1 default |
|---|---:|---:|---:|
| Specified packet fields preserved | 1 | 7 | 7 |
| Output with field loss | 3 | 0 | 0 |
| Read/write failure on valid input | 3 | 0 | 0 |
| Explicit rejection | 0 | 7 | 1 |
| Output on ambiguous/invalid source | 7 | 0 | 6 |

The seven valid inputs match Open Babel's preservation coverage. The repair
also prevents output on the seven declared invalid/ambiguous inputs; rejection
is not successful conversion and is not added to preservation coverage.
These are original development cases including the motivating failures,
not a sealed holdout, published-package proof, full CML conformance or evidence
of toolkit-wide superiority. Uninterpreted metadata, stereo and 3D coordinates
are outside scoring. The legacy four-decimal CML writer is unchanged; this
packet deliberately uses exactly representable coordinates within that profile.

The production baseline is `00fc7457`. The
[fixture recipe and manifest](https://github.com/kent-tokyo/chematic/blob/main/validation/cml_character_references/cases.json),
[baseline](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cml-character-refs-baseline-2026-10-22.json),
[candidate](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cml-character-refs-candidate-2026-10-22.json)
and [validation index](https://github.com/kent-tokyo/chematic/blob/main/validation/results/cml-character-refs-validation-2026-10-22.json)
record normalized input hashes, expected fields, engine identities and outcomes.
The recipe is original project data under `MIT OR Apache-2.0` and does not use
third-party or sealed corpus inputs.

## Reproduce

```sh
cargo test --locked -p chematic-mol --test cml_character_references
cargo build --locked -p chematic-mol --example file_io_semantic_probe
python3 -m scripts.validation.check_cml_character_references \
  --probe target/debug/examples/file_io_semantic_probe \
  --obabel /opt/homebrew/bin/obabel --output /tmp/cml-references-new.json
python3 -m pytest -q scripts/tests/test_cml_character_references.py
```

Results require a new output path and pin Open Babel 3.2.1. Additional Rust
regressions cover XML character bounds, single-pass literal entities, i8/u16
representation bounds, alias agreement, normal/strict reader parity and CDXML
attribute inspection with unchanged exact source. Python/WASM run the same
fourteen inputs through both public CML readers; their graph-only APIs do not
claim to retain source coordinates.

The broad WASM cross-binding suite still fails its macrocycle QED assertion
on this host, in debug and release candidates. The same assertion fails with
the recorded pre-CML artifact, whose original hash was verified. New CML,
format-parity, CDXML-document and ChemDraw suites pass. No QED code or expected
value was changed; the wider failure remains an explicit validation limitation.
