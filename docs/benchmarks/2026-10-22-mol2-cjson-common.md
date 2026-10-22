# MOL2 selected fields through ordinary CJSON conversion

Python `convert_format(source, "mol2", "cjson")` and WASM
`convert_common_format(source, "mol2", "cjson")` previously used the graph-only
bridge. Source xyz, partial charges and molecule name were omitted even though
the existing checked conversion report already retained and verified them.
These ordinary routes now share a selected-field projection over that report.

```python
output = chematic.convert_format(source, "mol2", "cjson")
# Equivalent selected-field API; accepts UTF-8 bytes too:
output = chematic.project_mol2_to_cjson(source.encode("utf-8"))
```

Rust exports `project_mol2_to_cjson` and
`project_mol2_to_cjson_with_limits`. They use the existing `allow_loss` report
internally, then require verified atom/bond correspondence and bitwise xyz
agreement. Any reported loss or uncertainty for parsed name, declared partial
charge values/method or projected formal charges rejects the projection. The
source method label is retained without recomputing or verifying that method.
Partial charges are never rounded into formal charge or generated here.

This is **selected-field preservation**, not the report's `strict` policy:
source text, atom names/types/IDs, substructure ownership, opaque sections,
implicit-hydrogen semantics and stereochemistry are still lost or unverified.
`convert_with_report(..., policy="strict")` still refuses these losses.
`allow_loss` remains available to inspect them and explicitly obtain a lossy
candidate. Same-format retention remains the source-byte workflow.

## Charge assignment and compatibility

A declared non-`NO_CHARGES` method requires a numeric value for every atom.
Missing one or all values refuses the ordinary projection. Nonzero numeric
columns without an assigned method also refuse. The full report now explicitly
marks a declared dataset with no source values as `partial_charges: unverified`;
previously the partial-charge observation was absent in that case.

`NO_CHARGES` with absent columns means no assigned dataset. Numeric zeros under
that declaration are unassigned placeholders, including signed zeros, and are
omitted without creating `partialCharges`. The report identifies the omitted
raw columns as `unassigned_zero_columns: dropped`. This keeps the existing
SMILES → generated MOL2 → CJSON workflow usable. This omission does not claim
bitwise retention of the unassigned raw placeholder columns. A declared
`USER_CHARGES` zero dataset, by contrast, survives with its original sign bits.

Python rejects explicit coordinates/charges, a nondefault name or a nonempty
comment on this route instead of replacing or silently ignoring source fields.
This tightens previous acceptance: callers needing intentional edits can use
the format-specific APIs, and callers accepting data loss can inspect the
explicit report. Default Rust/Python input/output limits are each 16 MiB with
10,000 atoms; other shared MOL2 line/count/bond limits apply. WASM tightens input
to 1,000,000 bytes and keeps the same atom/output limits. Native callers may
choose stricter limits; output is bounded before it is returned.

## Development packet

`validation/mol2_cjson_common/cases.json` contains 16 original synthetic MIT
fixtures in LF, CRLF and no-final-newline layouts: 48 inputs. Nine selected-field
cases (27 inputs) cover complete declared charges, GASTEIGER label retention,
negative partial charge on N.4 with projected formal charge +1, a UTF-8 name,
absent charges, tiny values, signed zero, opaque metadata and unassigned zero
placeholders. The UTF-8 name and vendor section are local-profile extensions;
they are not evidence of strict Tripos ASCII conformance.

Four ambiguous-charge fixtures (12 inputs) require explicit loss acceptance.
Three other fixtures (9 inputs) cover nonfinite values, duplicate source IDs and
multiple molecules outside the one-record API's domain. A multi-record MOL2
file is not intrinsically invalid for a batch converter. None of these fixtures
is a holdout or upstream issue reproduction.

All scored two-atom sources deliberately use sparse atom IDs 50/7 and bond ID
31. The oracle maps their declared bond endpoints by source ID into row indices.
This is the native reader's local profile. The original Tripos document says
atom IDs are reassigned on SYBYL reading and distinguishes those IDs from bond
endpoint fields; this packet does not establish general MOL2/SYBYL equivalence
for sparse numbering. Separate sequential-ID diagnostics are retained below.

The independent JSON observer separates indexed elements/formal charges/bonds,
xyz binary64 bits, parsed name and partial-charge values/absence from method
label/absence. It accepts absent empty bond sections as zero bonds, but the
scored sources declare one bond. No chemical equivalence, stereochemistry or
hydrogen assignment is verified.

| Route, selected-field inputs | JSON object | Indexed graph | xyz bits | Name | Charge values/absence | Method/absence |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Baseline Python common | 27/27 | 27/27 | 0/27 | 0/27 | 6/27 | 6/27 |
| Baseline WASM common | 27/27 | 27/27 | 0/27 | 0/27 | 6/27 | 6/27 |
| Candidate Python common | 27/27 | 27/27 | 27/27 | 27/27 | 27/27 | 27/27 |
| Candidate WASM common | 27/27 | 27/27 | 27/27 | 27/27 | 27/27 | 27/27 |
| Open Babel 3.2.1 `-imol2 -ocjson` | 18/27 | 0/27 | 18/27 | 18/27 | 14/27 | 14/27 |
| Open Babel 3.2.1 with `-xs` | 18/27 | 0/27 | 18/27 | 18/27 | 14/27 | 14/27 |

Open Babel emits no bonds for the scored sparse-ID sources. Its nine CRLF
variants have empty stdout/exit zero; the other 18 emit coordinate/name values
unchanged. All 14 emitted assigned-charge datasets retain measured numeric
bits and labels. Four emitted unassigned cases instead contain a `NO_CHARGES`
dataset of zeros. These distinctions, not a single aggregate ranking, explain
the table. Sequential-ID diagnostic copies retain the bond in Open Babel;
the zero graph score must not be generalized to ordinary sequential MOL2.
The sequential-ID `NO_CHARGES`/absent-column control emits nonzero values
`[0.21748719006663936, -0.21748719006663936]` under `partialCharges.NO_CHARGES`.
The source declares no assigned dataset; the candidate emits none. The control
records observed output, without attributing a particular charge algorithm.

The candidate refuses all 12 ambiguous and all 9 one-record-invalid/out-of-domain
inputs. Baseline common bridges return no output for 0/12 and 6/9 respectively.
Open Babel returns no text for 4/12 and 5/9; some emitted text is malformed JSON,
so nonempty text is not counted as a successful molecule conversion here. On
LF incomplete `USER_CHARGES`, for example, it emits a zero for the missing value.
The complete terminal outputs make that behavior distinguishable from a refusal.
No charge calculation, geometry generation, normalization or raw-copy options
are requested. `-xs` is the single-geometry CJSON option, not a precision knob.

Baseline Python uses main package `b9412a06` and its attested extension; WASM
uses attested `c713b21a`, whose unrelated CJSON-source branch does not change
MOL2 → CJSON common conversion. Candidate starts at main `d4208b9d`. The earlier
V3000/CJSON/CML/PDBQT PRs are not dependencies. Source/artifact identity checks
and complete results are recorded in the validation index.

This establishes the bounded ordinary route and clearer refusal behavior on
a small local development packet. It does not establish toolkit-wide
superiority, broad format coverage, speed, a release or full MOL2 semantics.

## Reproduction and references

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_mol2_cjson_common.py \
  --wasm-module /path/to/chematic_wasm.js --obabel /path/to/obabel \
  --output /tmp/mol2-cjson-common.json
```

`--baseline` records older artifacts without requiring candidate success.
Subprocess time and both output pipes are bounded. Native tests check source
fields, explicit report policy and inclusive input/atom/bond/output bounds.
Python and real Node/WASM tests cover common routes, refusals, own-writer zero
placeholders and runtime limits; Node covers 10,000/10,001 atoms. Observer
negative controls separate changed graph, coordinates, method, charge values,
name, invented datasets and malformed JSON. CI discovers the suites and runs
the observer tests.

Full baseline/candidate results, sequential-ID diagnostics, tool help and
compressed logs are indexed by SHA-256 in
`validation/results/mol2-cjson-common-validation-2026-10-22.json`.

References: [Tripos MOL2 format, original SYBYL document mirrored by DockRMSD](https://zhanggroup.org/DockRMSD/mol2.pdf),
[Open Chemistry CJSON description](https://github.com/OpenChemistry/chemicaljson/blob/main/chemicaljson.md)
and [CJSON schema](https://github.com/OpenChemistry/chemicaljson/blob/main/cjson.schema).
