# Checked CML to CJSON projection

CML-to-CJSON conversion now verifies a selected graph/Cartesian/name profile
before returning output. Previously the ordinary Python and WASM route retained
3D coordinates but omitted names and unrepresented attributes. The CLI used a
separate graph bridge that also omitted geometry. All three routes now use one
native projector and refuse payload outside the profile.

```python
import json
import chematic

source = '''<molecule xmlns="http://www.xml-cml.org/schema">
<name>water &amp; oxygen</name>
<atom id="oxygen" elementType="O" x3="-0.0" y3="1.25" z3="2.5"/>
</molecule>'''
output = chematic.project_cml_to_cjson(source)  # UTF-8 str or bytes
assert json.loads(output)["name"] == "water & oxygen"
assert chematic.convert_format(source, "cml", "cjson") == output
```

Rust exposes `project_cml_to_cjson` and `project_cml_to_cjson_with_limits`.
WASM exposes `project_cml_to_cjson` and routes `convert_common_format` through
it. CLI `convert --input-format cml --output-format cjson` uses the same check.
The report API's supported routes are unchanged: it does not yet offer a
CML-to-CJSON report. Python's ordinary conversion wrapper remains str-only;
the direct projector also accepts UTF-8 bytes. Coordinate, charge, name and
comment overrides are refused on this source-preserving projection.

## Accepted fields and limits

- Exactly one nonempty molecule, either a molecule root or one direct molecule
  in a `cml` wrapper. Namespaced prefixes and the legacy unnamespaced form are
  accepted; foreign and mixed-namespace payload is refused.
- Source-order atom elements, formal charges and isotope masses; unique atom IDs
  resolve forward or reversed bond references. Single, double and triple bonds
  retain their indexed endpoints and orders. XML IDs are document-local
  references: molecule/atom/bond ID spellings are not retained in CJSON.
- Complete, finite x3/y3/z3 values on every atom, or no XYZ on any atom. Values
  are interpreted as Cartesian angstroms. Numeric comparison uses binary64 bits,
  including signed zero; no geometry is generated when it is absent.
- An optional molecule `title` and/or direct text-only `name`. Repeated names
  must agree after XML decoding. Presence and decoded value, including an
  explicit empty name, CR/LF/tab and Unicode, survive as root CJSON `name`.
- `hydrogenCount` may repeat a count derivable by the local CJSON reader. After
  reparsing the candidate, the non-node H count must match for every atom.
  A disambiguating count, including one changed by explicit H/isotope nodes,
  is refused. This checks local inference, not an independent chemistry model.

2D coordinates (including documents with both 2D and 3D), packed/fractional
coordinates, aromatic/query/zero bonds, stereo declarations, spin multiplicity,
units, properties, crystal data and unknown attributes/elements are refused.
Zero isotope labels are refused. Comments, processing instructions, whitespace,
CDATA versus text spelling and namespace prefixes are representation details;
they are not copied. The selected fields are not an exact-source contract.

The existing bounded graph/3D reader is unchanged and remains available when a
caller deliberately wants only that subset. Use the separate same-format source
retention contract when original XML, opaque annotations or ID spellings matter.
The new common-route refusals are intentional changes from earlier acceptance
that silently omitted fields. They are not counted as preserved conversions.

Byte, physical/normalized line, atom/bond and XML node/depth bounds apply before
projection. Output size is checked before returning. Candidate CJSON is reparsed
and checked for indexed graph, isotope labels, effective local H counts, name and
XYZ presence/bits. Chemical validity, full CML schema validity and stereo inferred
from geometry are not certified. The
[Open Babel CML documentation](https://openbabel.org/docs/FileFormats/Chemical_Markup_Language.html)
lists a broader set of interpreted fields; this profile makes no full-format
coverage claim. CJSON field context is documented by
[OpenChemistry](https://github.com/OpenChemistry/chemicaljson).

## Recorded comparison

`validation/cml_cjson_projection/cases.json` contains 11 original synthetic
positive fixtures and 25 refusal fixtures, each evaluated in compact, LF and
CRLF layout: 33 positive and 75 outside-profile evaluations. They are not 108
independent molecules or an unseen holdout. The fixtures target names, precision,
signed zero, missing geometry, charges/isotopes, H-count agreement, XML prefixes
and forward references. An independent ElementTree source observer and typed
CJSON observer compare authored selected fields. Native tests additionally check
H-count refusal, resource bounds and numeric extremes; those extremes are not
included in the paired molecular packet.

| Route/build | Indexed graph | Isotopes | XYZ bits or absence | Name presence/value | All four |
|---|---:|---:|---:|---:|---:|
| Candidate Python, WASM, CLI, each | 33/33 | 33/33 | 33/33 | 33/33 | 33/33 |
| Earlier local `ed050821`, Python/WASM, each | 33/33 | 33/33 | 33/33 | 15/33 | 15/33 |
| Earlier local `ed050821`, CLI | 20/33 | 32/33 | 6/33 | 15/33 | 0/33 |
| Open Babel 3.2.1 default | 16/33 | 14/33 | 16/33 | 6/33 | 6/33 |
| Open Babel 3.2.1 `-xs` | 16/33 | 14/33 | 16/33 | 6/33 | 6/33 |

There are 18 evaluations with declared names. Candidate retains all 18; the
previous local routes retain zero and Open Babel retains six. The previous
routes' 15 name matches are no-name controls. This metric distinguishes absent
from explicitly empty names. Comparator generated fields, such as partial
charges, are outside the four-field score; a full score is not a claim that no
other output fields were added.

The measured Open Babel installation emits no output in 17 positive evaluations:
11 CRLF inputs, four compact/LF inputs lacking 3D geometry, and two compact/LF
empty-name inputs terminated by SIGSEGV. Its process returned exit 0 for the
other 15 no-output inputs. These are observations of this macOS executable, not
claims about every Open Babel build. The prefixed wrapper also records XML parser
errors for the geometry-free case. Candidate explicitly refuses all 75
outside-profile evaluations without output; Open Babel emits nothing for 27.
Acceptance outside this local profile does not establish a CML-spec violation,
and safe refusal is not a successful preservation result.

Raw records retain all input/output/error text, command arguments, process bounds
(10 seconds and combined stdout/stderr 1 MiB), comparator version, runtime hashes,
fixture/runner hashes and candidate source hashes. The baseline uses preserved
local PR #1186 runtime snapshots, not a released or main build; its CLI hash was
matched to that PR's recorded candidate artifact. Candidate binaries were built
from the listed working-tree source hashes on parent `255edf73`, before commit.
No timing or toolkit-wide superiority conclusion follows from this packet.

- `validation/results/cml-cjson-projection-baseline-2026-10-24.json.gz`
- `validation/results/cml-cjson-projection-candidate-2026-10-24.json.gz`

Reproduce with a fresh Python extension on `PYTHONPATH`, a fresh Node-target WASM
module, and a CLI built from the same source:

```sh
python3 scripts/validation/check_cml_cjson_projection.py \
  --wasm-module /path/to/chematic_wasm.js \
  --cli /path/to/chematic \
  --output /tmp/cml-cjson-candidate.json.gz
python3 -m unittest scripts/tests/test_cml_cjson_projection.py
```

The dated CML 3D and namespace evidence remains frozen. Their current runners
expect explicit refusals for the 2D/annotation fixtures that the older routes
accepted by dropping fields; their historical positive counts are not current
common-route acceptance claims.
