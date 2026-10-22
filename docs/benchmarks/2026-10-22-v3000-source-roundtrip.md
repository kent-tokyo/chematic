# V3000 source retention through ordinary conversion

The Python/WASM common V3000 → V3000 paths previously reconstructed a molecule
without retaining input coordinates or metadata. WASM's explicit
`roundtrip_mol_v3000_block` retained selected metadata but used the 2D writer:
z was zeroed, source IDs normalized, groups reordered and x/y rounded to four
decimals. Small nonzero z values could not be recovered from the perception
conformer either, because its zero tolerance serves a different purpose.

The same-format routes now validate and return the exact UTF-8 source:

```python
retained = chematic.convert_format(source, "v3000", "mol_v3000")
assert retained == source
assert chematic.roundtrip_mol_v3000_source(source.encode("utf-8")) == source
```

```javascript
const retained = wasm.roundtrip_mol_v3000_block(source);
if (retained !== source) throw new Error("source changed");
// convert_common_format(source, "v3000", "v3000") uses the same path.
```

Rust exports `roundtrip_mol_v3000_source`,
`roundtrip_mol_v3000_source_with_limits` and `V3000SourceLimits`. The source is
immutable; original x/y/z tokens, IDs, atom maps, SGROUP/COLLECTION ordering,
attributes, header text and line endings survive without reserialization.
Opaque chemistry/stereo is retained without editing or verification. This is
not a new editable molecule model, general CTfile validator or stereochemical
fidelity claim. Existing graph readers and all regenerated MOL writers keep
their behavior, including their precision and 2D/3D boundaries. The
[BIOVIA CTfile Formats 2020](https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf)
describes V3000's free-form records and continuation syntax; these APIs accept
a narrower local profile through the existing native parser.

## Compatibility and bounds

The explicit WASM round-trip API now returns source IDs, original formatting
and original group order; it no longer promises the normalized text produced
by the old 2D writer. Python same-format conversion rejects coordinate, charge,
name and comment overrides, including explicit empty arrays. Default arguments
remain accepted. Other source/target routes keep their existing behavior.

The retained input must contain one complete CTAB and terminate with `M  END`.
Trailing SDF fields/delimiters, extra CTABs or V30 records outside that CTAB,
dangling continuations, duplicate/zero atom or bond IDs, nonfinite coordinates
and malformed native core structures fail without output. Unknown fields are
not interpreted; acceptance is not full compliance or consistency verification
of extensions, query/polymer semantics or coordinate-derived stereo. The
existing parser's continuation/projection limitations remain applicable.

Default Rust/Python bounds are 16 MiB input, 100,000 atoms, 1,000,000 physical
lines and 1 MiB per physical line. Byte/line and source-ID/atom-count checks
precede graph parsing; existing native atom/bond caps also apply. WASM tightens
input to 1,000,000 bytes and atoms to 10,000. Output is the checked source and
therefore fits those input bounds. Python bytes must decode as UTF-8.

## Development packet and comparison

Six original MIT one/two-atom inputs cover tiny z, full precision 3D, a plane
offset from z=0, sparse IDs with groups, radical/valence attributes and a
token-boundary continuation. LF/CRLF/no-final-newline give **18** variants.
`validation/v3000_source_roundtrip/cases.json` records provenance, license,
minimization, units, oracle and normalized SHA-256. These fixtures were used
during implementation and are not an independent holdout.

Baseline Python uses an attested native extension based on main `b9412a06`
with main package source. Baseline WASM is attested commit `c713b21a`; that
commit's unrelated CJSON source route does not change the V3000 paths under
test. The final candidate incorporates the subsequent main CIP guard
`d4208b9d`, which does not modify this I/O path. Artifact/source hashes and
integration checks are in the validation index.

An independent limited token observer compares source-ordered element labels,
raw x/y/z binary64 bits (including signed zero) and IDs. This is a coordinate
retention metric, not graph isomorphism or chemistry interpretation. Output
without an observable atom table is counted separately. Formatting alone does
not count as coordinate loss; exact source is a distinct metric.

| Route | Atom table observed | Exact source | Exact coordinate bits |
| --- | ---: | ---: | ---: |
| Baseline Python common conversion | 18/18 | 0/18 | 0/18 |
| Baseline WASM common conversion | 18/18 | 0/18 | 0/18 |
| Baseline WASM explicit round trip | 18/18 | 0/18 | 0/18 |
| Candidate, each of these three routes | 18/18 | 18/18 | 18/18 |
| Open Babel 3.2.1 `-imol -omol -x3` | 17/18 | 0/18 | 3/18 |
| Open Babel 3.2.1 `-imol -omol -x3 -xa` | 17/18 | 0/18 | 3/18 |
| Open Babel raw-copy SDF carrier control | 18/18 | 12/18 | 18/18 |

Open Babel's ordinary V3000 output rounds coordinates; the CRLF continuation
input produces no output despite exit code zero. Full terminal outcomes are
retained, so that no-output case is distinguishable from rounding. The raw-copy
control uses a seekable `.sdf` file with explicit `$$$$` framing and
`-isdf source.sdf -ocopy -aT`. Its expected exact bytes include that carrier
framing. Title-only reading does not claim graph validation. All coordinate
values survive in that control; six CRLF variants differ in source bytes.
Thus the ordinary-writer comparison does not establish absence of a raw-copy
workflow in Open Babel or toolkit-wide superiority.

An initial development fixture included a nonstandard atom `LABEL` extension
that Open Babel did not emit. It was removed while retaining standard RAD/VAL
attributes. The earlier packet and both complete outcomes remain in
`v3000-source-roundtrip-development-v1-2026-10-22.json.gz`; this revision is not
silently treated as an independent packet. The final results above use the
revised fixture. This work establishes the local immutable-source workflow,
not cross-format/editable preservation, general CTfile coverage, performance,
a release, or interpreted SGROUP/stereo fidelity.

## Reproduction

With a built Python package and a Node-target WASM module:

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_v3000_source_roundtrip.py \
  --wasm-module /path/to/chematic_wasm.js --obabel /path/to/obabel \
  --output /tmp/v3000-source-roundtrip.json
```

`--baseline` records older artifacts without requiring candidate success.
Subprocess time/output are bounded. The observer tests detect lost z, rounded
coordinates and signed-zero loss separately from harmless headers/ID changes.
Native tests cover single-block boundaries, malformed input and inclusive
custom limits. Real Node/WASM tests include the 10,000/10,001-atom boundary.
Python tests cover both public routes, bytes/UTF-8 and override refusal. CI
runs the observer tests and discovers the new binding/native tests.

Complete per-input results are retained in
`validation/results/v3000-source-roundtrip-{baseline,candidate}-2026-10-22.json.gz`;
the validation index records decompressed hashes, sources, artifacts and logs.
