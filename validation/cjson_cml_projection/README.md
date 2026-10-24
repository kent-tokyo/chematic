# Checked CJSON → CML projection, 2026-10-24

This exposed development packet measures graph, Cartesian and name retention
in a cross-format conversion. Eight original common fixtures × LF/CRLF/compact
give 24 evaluations, not 24 independent molecules. A ninth fixture × three
layouts is a separate compatibility lane: CJSON order 4 is chematic's local
aromatic alias, while other implementations can interpret it differently.
No sealed holdout, chemical-validity or toolkit-wide superiority claim is made.

| Common lane, 24 evaluations | Indexed graph | Isotopes | Coordinate bits / absence | Name requirement | All selected fields |
| --- | ---: | ---: | ---: | ---: | ---: |
| Candidate Python | 24 | 24 | 24 | 24 | 24 |
| Candidate WASM | 24 | 24 | 24 | 24 | 24 |
| Candidate CLI | 24 | 24 | 24 | 24 | 24 |
| Earlier local snapshots, each | 24 | 24 | 3 | 9 | 0 |
| Open Babel 3.2.1 default | 15 | 18 | 15 | 9 | 3 |
| Open Babel 3.2.1 `-xp -xm` | 15 | 18 | 15 | 9 | 3 |

The name requirement includes nine variants with no declared source name.
Open Babel does not retain the declared name on any of the other 15 variants
(three of those also have no molecular output). This avoids presenting the
nine no-name controls as nine preserved declared names.

The baseline uses the local `c693dfb7` builds from the preceding CML retention
work, not a published release or separately attested current-main baseline.
Its CJSON → CML route still writes graph-only CML, omitting XYZ and source names.
Artifact hashes identify the precise Python, Node WASM and CLI binaries used.

All Open Babel common evaluations return exit 0. The three absent-geometry
variants emit no molecular output and report `Coordinates are not specified` /
`0 molecules converted`; they are counted as missing output, not preservation.
Other outputs are XML, with six-decimal coordinates. Precision and signed-zero/
tiny-coordinate differences are observed without alignment or tolerance. In
this local reader profile, charge/isotope arrays and declared source names are
not all retained by Open Babel. This is measured field coverage, not a claim
that every reader must implement chematic's extensions. The `-xp -xm` mode does
not change these results. Generated metadata outside the scored fields is
not penalized. Empty-name presence is scored explicitly; absent source names
do not fail merely because a comparator adds a title.

The compatibility lane retains all selected fields for candidate Python/WASM/
CLI on 3/3 evaluations. Open Babel's different bond-order-4 interpretation is
recorded separately and is not used to rank common-format correctness.

The independent Python ElementTree observer resolves emitted atom IDs in
source order, normalizes only undirected bond endpoint orientation, and checks
elements, formal charges, isotopes, bond orders, optional finite binary64 XYZ
and source name. It accepts either CML `<name>` text or a consistent molecule
`title`. No toolkit supplies the expected answers. Missing output is distinct
from malformed XML and field mismatch. Five observer controls catch changed
bonds, missing coordinate components, signed-zero/rounding loss, name loss and
invalid/absent XML. Every child process has a 10-second / 1 MiB combined-output
budget.

The [baseline report](../results/cjson-cml-projection-baseline-2026-10-24.json.gz)
and [candidate report](../results/cjson-cml-projection-candidate-2026-10-24.json.gz)
retain every input, output, error, command, artifact/source identity and per-field
observation. Reproduce against fresh Python, CLI and Node WASM builds:

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_cjson_cml_projection.py \
  --cli /path/to/chematic \
  --wasm-module /path/to/chematic_wasm.js \
  --output /tmp/cjson-cml-candidate.json.gz
```

Ten refused-input fixtures separately cover partial charges, stereo/unit-cell/
vendor/2D extensions, conflicting or non-XML names, empty graphs, duplicate JSON
keys and incomplete XYZ. These are API regression tests, not preservation wins
or paired Open Babel measurements. Tests also cover overrides, Unicode, byte/
atom/output limits, finite numeric extremes, and CLI output-file safety.

The projector reparses and verifies the emitted selected fields. It does not
verify source chemistry, inferred H counts or coordinate-derived stereo. CML
aromatic atom flags and total H declarations are derived by the existing writer/
reader; they are not explicit CJSON source fields. Original JSON bytes and
opaque metadata are not silently treated as preserved: opaque fields cause
refusal. Exact same-format retention and report-based conversion keep their
existing contracts.
