# V2000 source coordinates: bounded local candidate

This packet measures V2000 MOL/SDF → CJSON on an authored synthetic corpus.
It is not a release, full-format fidelity result or performance comparison.
The implementation base is `38b4c92c430d9b3bfe70e47e97293b9017ced990`.

## Observations

The 10 original profiles are expanded across LF/CRLF/no-final-newline and
MOL/one-record-SDF framing: 60 routes. Source triples are independently decoded
from the fixed-width fields, then output CJSON is observed for indexed elements,
formal charges, isotopes, bond orders and exact binary64 coordinate bits. The
observer is shared with the separately tested V3000 packet.

| Artifact | Coordinate arrays | Exact coordinate bits | All declared fields |
| --- | ---: | ---: | ---: |
| Previous local Python/WASM packages | 0/60 | 0/60 | 0/60 |
| Candidate Python | 60/60 | 60/60 | 60/60 |
| Candidate real Node/WASM | 60/60 | 60/60 | 60/60 |
| Open Babel 3.2.1 default CJSON | 54/60 | 54/60 | 48/60 |
| Open Babel 3.2.1 CJSON with `-xs` | 54/60 | 54/60 | 48/60 |

Open Babel's six tiny-z/2D-header cases emit no CJSON and report that no 3D
coordinates exist. These are explicit refusals, not silent wrong outputs.
The six isotope-ion variants emit valid CJSON without retaining source isotope
labels. Other inferred Open Babel metadata/partial charges are outside this
contract. `atoms.isotopes` is chematic's existing local CJSON profile; this is
a source-fidelity check, not a claim that all ChemicalJSON implementations
promise that field. The 10 profiles are not 60 independent chemical structures.

The fixtures include normal four-decimal records and deliberate wider-precision
or scientific tokens within 10-byte coordinate fields. The latter are accepted
extension tests. RDKit 2026.09.1 was inspected separately: it refuses two such
profiles and reads one wider-precision profile differently. It is not treated
as a definitive coordinate oracle for those extensions. Source-column decoding
and exact input bytes define this coordinate-retention contract.

## Artifacts and reproduction

`baseline.json.gz` and `candidate.json.gz` retain every command, stdout/stderr,
resource/exit status, source hash, input hash, observation and tool binary hash.
Gzip uses mtime 0. The candidate's head is the pre-commit base; its dirty tested
source files and Python/WASM binaries are identified separately by hashes.

Baseline packages came from PR #1177 (`80fe1a10943cb49f74e0b24d4e26d514da0c85ac`)
and contain that PR's unrelated PDBQT feature. Its Python conversion wrapper,
native V2000 reader and WASM common-format source match the main base; the
baseline is not represented as a published main package. Baseline report source
identities are omitted rather than attributing candidate files to that build.

```sh
PYTHONPATH=/path/to/candidate-python python3 scripts/validation/check_v2000_coordinate_projection.py \
  --wasm-module /path/to/node/chematic_wasm.js \
  --obabel /opt/homebrew/bin/obabel \
  --output /tmp/v2000-coordinate-candidate.json
```

The old package can be measured with the same command and `--baseline`.
Comparators receive the same V2000 bytes and CJSON target. `-xs` is recorded as
a separate Open Babel lane to avoid conflating default conformer handling.

Validation passed: 822 native molecular-library tests; 3 new native contracts;
124 Python V2000/V3000/bridge tests and 59 related CML/SDF tests; real Node/WASM
V2000, V3000, SDF-document and MOL2→CJSON suites. The new native/Python/WASM
contracts reject all 18 malformed coordinate fields (3 axes × 6 tokens).
Formatting, Clippy with warnings denied, stubs, API compatibility/inventory,
cross-binding manifest, refactor budget, script layout, integration-test discovery
and strict MkDocs build passed. The docs build used tracked files plus the new
page, excluding unrelated ignored/private documentation.
