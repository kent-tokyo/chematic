# Ordinary CML source retention, 2026-10-24

This exposed development comparison tests **ordinary CML → CML conversion**,
not the report API that already retained source. It reuses the original
[six-document packet](../cml_source_preservation/cases.json) with LF, CRLF and
no-final-newline forms: 18 evaluations, not 18 independent chemical structures.
The cases cover 2D, simultaneous 2D/3D, high-precision 3D, repeated properties /
foreign extensions, namespace prefixes and comments/processing instructions.

| Path | Exact source | Selected fields retained | Graph mismatch | Other selected-field loss |
| --- | ---: | ---: | ---: | ---: |
| Candidate Python | 18/18 | 18/18 | 0/18 | 0/18 |
| Candidate WASM | 18/18 | 18/18 | 0/18 | 0/18 |
| Candidate CLI | 18/18 | 18/18 | 0/18 | 0/18 |
| Earlier local Python/WASM/CLI, each | 0/18 | 0/18 | 18/18 | 0/18 |
| Open Babel 3.2.1 default | 0/18 | 9/18 | 0/18 | 9/18 |
| Open Babel 3.2.1 `-xp -xm` | 0/18 | 9/18 | 0/18 | 9/18 |

The earlier local snapshot is commit `6b4afbb9`, built for the CJSON-output
work, not a released or separately attested current-main baseline. Its CML
self-conversion still uses the old graph bridge. All 18 outputs in that lane
are well-formed XML but lose the source bond: the old graph reader does not
read the inline `<bondArray><bond .../></bondArray>` used by these fixtures.
It also omits coordinates and metadata. The selected-graph observer identifies
the graph mismatch, rather than calling the output absent or preserved.

All Open Babel evaluations return one molecular XML output. Three source cases
× three newline forms lose selected fields in both modes: dual-dimensional
input loses x2/y2 when selecting 3D; the precise-3D case changes coordinates by
more than the absolute 1e-12 Å tolerance; the opaque-properties case loses a
foreign atom attribute and the ordered property/extension payload. Added
metadata outside source subtrees is not counted as a source-information loss.
Molecule root ID/title, implicit hydrogens, chemical/schema validity and stereo
interpretation are outside the selected-field scorer. Comments, processing
instructions, prefix/number spellings and newline preservation are exact-text
diagnostics, not additional chemical-structure successes.

The independent standard-library XML observer uses source atom identity,
element/charge/isotope and one single C–O bond, coordinate component presence
and values, source atom attributes, and expanded property/extension subtrees.
No toolkit is the oracle. Six existing observer negative controls plus three
new classification controls detect removed dimensions/attributes, reordered
properties, changed graphs, malformed XML and empty output. Child processes
are bounded to 10 seconds and 1 MiB combined stdout/stderr.

Every input, output, error and artifact/source identity is retained in the
[baseline report](../results/cml-common-retention-baseline-2026-10-24.json.gz) and
[candidate report](../results/cml-common-retention-candidate-2026-10-24.json.gz).
The new ordinary conversion returns raw source directly through the same
bounded XML validator used by the existing report API. It retains unknown
content opaquely and makes no claim that its chemistry or references were
validated. This result does not establish cross-format fidelity, editing
support, holdout performance or Open Babel-wide superiority.

Reproduce with fresh builds:

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_cml_common_retention.py \
  --cli /path/to/chematic \
  --wasm-module /path/to/chematic_wasm.js \
  --output /tmp/cml-common-candidate.json.gz
```

Malformed XML/prolog/framing, unsupported nested molecule documents, source
override rejection, UTF-8 handling, atom/byte/depth/output boundaries, report
parity and existing-file safety are tested separately. They are not added to
the paired preservation counts. Python's common wrapper remains str-only;
the direct retention API additionally accepts UTF-8 bytes.
