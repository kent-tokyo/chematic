# Namespace-aware CML projection: local source evidence

The candidate is based on `38b4c92c430d9b3bfe70e47e97293b9017ced990` and fixes
the Cartesian/CJSON projection. It is not a published release or a full CML,
stereo, metadata or performance comparison.

## Same-input, same-target observations

The authored packet has 20 positive profiles × compact/LF/CRLF = 60 positive
routes, and 28 negative profiles × three layouts = 84 refusal routes. These
are layout variants of synthetic inputs, not 144 independent molecules.
Python ElementTree independently verifies source expanded names and graph fields;
a separate CJSON observer checks indexed elements, formal charges, isotope
annotations, bond orders and exact coordinate bits in the output.

| Artifact | All declared fields, valid inputs | Graph + coordinates | No output on valid inputs | No output on negative inputs |
| --- | ---: | ---: | ---: | ---: |
| Previous local Python / WASM | 21/60 | 21/60 | 27/60 | 51/84 |
| Candidate Python | 60/60 | 60/60 | 0/60 | 84/84 |
| Candidate real Node/WASM | 60/60 | 60/60 | 0/60 | 84/84 |
| Open Babel 3.2.1 CJSON | 28/60 | 30/60 | 20/60 | 42/84 |
| Open Babel 3.2.1 CJSON `-xs` | 28/60 | 30/60 | 20/60 | 42/84 |

Every positive candidate result passes every declared-field observation. Every
negative candidate result raises an explicit error before returning output.
Both Open Babel lanes show the following positive-input outcomes:

- 28 outputs match every declared field.
- 10 outputs change the selected graph and coordinates: foreign atoms,
  foreign ancestors, namespace-prefix rebinding, an empty-namespace atom and
  formula composition are flattened into molecular graph atoms in compact/LF
  layouts (five profiles × two layouts).
- Two labelled-ion outputs omit isotope annotations while retaining graph and
  coordinates. This uses chematic's existing local CJSON isotope profile, not
  a claim that all ChemicalJSON implementations promise that field.
- 20 routes emit no CJSON. 18 are CRLF layouts; the other two are compact/LF
  empty-namespace notes containing a separate metadata molecule. Exit codes,
  exact stderr and empty stdout are retained. Some say only `0 molecules
  converted`, so these are not classified as explicit typed refusals.

CRLF behavior is separately reproduced with both a real file and stdin in
`crlf-file-reproduction.json`; the exact input is `default-crlf.cml.gz`. The
compact/LF-only subset still shows 10 selected-graph differences and two isotope
omissions, so namespace contamination is separate from newline handling.
Some otherwise matching wrapper outputs have parser diagnostics on stderr;
those diagnostics are preserved rather than being called clean conversions.
Acceptance of a negative input means acceptance outside this local profile,
not automatically a CML-specification violation.

## Reproduce and provenance

`baseline.json.gz` and `candidate.json.gz` use deterministic gzip mtime 0 and
retain every source string, expected field, command, output, error/resource/exit
status, observation and binary/source identity. The report `head` is the current
checkout base at execution, with dirty candidate sources identified separately
by hashes; it is not the provenance commit of the baseline package.

Baseline packages were built during PR #1178 at
`bc4148244902045597cfb7d1834d76a0d560cf79`; that PR changes the unrelated V2000
coordinate route. Its CML reader/projection is the main-base implementation.
The baseline is not presented as a published main package.

```sh
PYTHONPATH=/path/to/candidate-python python3 scripts/validation/check_cml_namespace_projection.py \
  --wasm-module /path/to/node/chematic_wasm.js \
  --obabel /opt/homebrew/bin/obabel \
  --output /tmp/cml-namespace-candidate.json
```

Use the previous package with `--baseline` to observe the old behavior.
Both comparators receive identical CML bytes and the CJSON target. No coordinate
generation is requested. Partial-charge inference, opaque metadata and stereo
are outside the observation contract.

Validation passed: 822 native library tests and 32 related CML integration tests;
217 Python CML/bridge/SDF tests; new namespace and existing CML-3D Node/WASM
suites; the seven independent CJSON-observer tests. Clippy with warnings denied,
formatting, stub/API compatibility/inventory, cross-binding manifest, refactor
budget, script layout and integration discovery pass. Strict MkDocs builds use
tracked documentation plus the new page, excluding unrelated ignored/private
documentation.
