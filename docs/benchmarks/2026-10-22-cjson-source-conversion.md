# CJSON source retention in the common conversion bridge

The ordinary Python and WASM conversion APIs previously reconstructed CJSON
from a projected molecule, even when both format arguments were `cjson`.
This discarded fields already supported by the immutable record APIs. The
Python route kept 3D coordinates but dropped charge datasets and metadata; the
WASM topology bridge also dropped geometry.

Same-format CJSON now uses the existing checked source record path:

```python
retained = chematic.convert_format(source, ".CJSON", "cjson")
assert retained == source
```

```javascript
const retained = wasm.convert_common_format(source, "cjson", ".CJSON");
if (retained !== source) throw new Error("source changed");
```

Both APIs validate the [local immutable CJSON profile](../format-capabilities.md#immutable-chemicaljson-records)
before returning the exact UTF-8 source. Coordinates, declared charge datasets,
names, isotopes, formal charges, opaque extension fields, large integer tokens,
numeric spelling, whitespace and line endings survive because no molecule is
serialized back over the input. This is source retention, not verification or
editing of opaque chemistry. There is no new native parser, schema or export.
The [Chemical JSON project](https://github.com/OpenChemistry/chemicaljson)
describes arrays mapped to atoms/bonds and ongoing representation extensions;
this local profile is intentionally narrower than arbitrary CJSON documents.

## Compatibility and limits

Python rejects coordinate, charge, name and comment overrides on this path,
including explicit empty coordinate/charge arrays. This prevents an override
from being silently ignored while the source is returned. The existing defaults
(`None`, `None`, `"LIG"`, `""`) are accepted. Use an editable molecular API when
editing is intended, or the checked conversion report for documented cross-format
losses. Other target formats retain their existing behavior.

The common same-format bridge now rejects duplicate decoded keys, unsupported
version declarations, inconsistent arrays, invalid/self/duplicate bonds,
out-of-range atom indices and unsupported bond orders according to the existing
immutable profile. These inputs could previously be accepted and silently
normalized by the permissive graph reader. The legacy `from_cjson` and
`mol_from_cjson` readers keep their behavior. Unknown fields are retained only
opaquely: coordinate stereo, vendor references and their consistency are not
interpreted or verified.

Python checks its 16 MiB UTF-8 input cap before source validation. WASM checks
1,000,000 input bytes and delegates to its bounded record reader with 10,000
atoms and existing JSON parser budgets. No validation is
bypassed by the same-format return. The output is the input and therefore also
fits the existing WASM output budget.

## Paired development packet

`validation/cjson_source_conversion/cases.json` contains six original MIT
inputs, each with two atoms, provenance, minimization, expected contract and
SHA-256. They isolate a declared charge dataset, isotope/formal-charge data,
multiple charge methods, alternative coordinates, opaque references/Unicode
and high-precision numbers. LF, CRLF and compact/escaped JSON give **18**
variants. They were used during implementation, not sealed as an independent
holdout.

The baseline is main `a83ed2653037a7b17033e03ee22501a1ddda8bb9`.
The previously attested baseline Python extension was reused after verifying
its hash and matching the package source to that commit. A fresh baseline WASM
build was made before changing production sources. Candidate artifacts were
built from the scoped patch and are identified by hash in the validation index.
Before PR creation, the candidate was rebased onto main `b9412a06` (the
ETKDG conformer-pooling change, which does not modify the CJSON paths) and both
bindings were rebuilt and revalidated. The pinned baseline remains unchanged.

| Engine/route | Successful conversions | Exact source | All source JSON fields retained |
| --- | ---: | ---: | ---: |
| Baseline Python bridge | 18/18 | 0/18 | 0/18 |
| Candidate Python bridge | 18/18 | 18/18 | 18/18 |
| Baseline Node/WASM bridge | 18/18 | 0/18 | 0/18 |
| Candidate Node/WASM bridge | 18/18 | 18/18 | 18/18 |
| Open Babel 3.2.1 default | 18/18 | 0/18 | 0/18 |
| Open Babel 3.2.1 `-xm` | 18/18 | 0/18 | 0/18 |
| Open Babel 3.2.1 `-xv` | 18/18 | 0/18 | 0/18 |

The commands are `obabel -icjson -ocjson`, with the listed output option when
applicable. Every invocation completed and emitted valid JSON. The independent
observer requires each source object field and value but allows additional
candidate fields. It checks array order/length and numeric values, including
signed zero and integers above JavaScript's exact integer range. Formatting
alone does not count as JSON-field loss. A second metric compares the complete
source text. No graph-isomorphism or chemical interpretation is inferred from
these JSON observations.

Representative changed source paths in Open Babel default output were
`/partialCharges/USER_CHARGES`, `/atoms/formalCharges`, `/atoms/isotopes`,
`/atoms/coords/2d`, `/atoms/coords/3dSets`, `/atoms/vendorRefs` and
`/properties/name`. The precision input also changed a coordinate from
`-1.4175225513503917` to `-1.4175225513503916`. Full per-variant outputs,
changed paths, exit codes and stdout/stderr are retained in both paired JSON
reports. Empty explicit bond arrays are counted as source fields; their absence
alone is not a claim of graph loss.

These measurements establish the local same-format source-retention contract.
They do not establish general CJSON compatibility, cross-format preservation,
editing support, chemistry correctness, toolkit-wide superiority or a published
release. The candidate retains extensions without needing to understand them;
Open Babel's reconstruction workflow has a different scope.

## Reproduction

Build the Python extension and a Node-target WASM module, then run:

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_cjson_source_conversion.py \
  --wasm-module /path/to/chematic_wasm.js --obabel /path/to/obabel \
  --output /tmp/cjson-source-conversion.json
```

`--baseline` records older artifacts without requiring candidate retention. All
subprocesses have bounded time/output. The unit-tested observer distinguishes
missing fields, charge precision, atom order, signed zero and large integer loss
from harmless added fields. CI discovers the Python and real Node/WASM contract
tests; the observer tests run in the Python CI job. The Node contract also checks
WASM input/atom limits and malformed typed fields/duplicate decoded keys.

Evidence:

- `validation/results/cjson-source-conversion-baseline-2026-10-22.json`
- `validation/results/cjson-source-conversion-candidate-2026-10-22.json`
- `validation/results/cjson-source-conversion-validation-2026-10-22.json`
- Compressed build/test logs under `validation/results/cjson-source-conversion-2026-10-22/`
