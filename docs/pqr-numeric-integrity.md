# PQR atom-field integrity

`write_pqr_checked` writes supplied atom records without rounding coordinates,
partial charges or radii. It uses 17 significant digits and verifies that each
decimal token reads back to the same finite binary64 value, including signed
zero and subnormal values. No values are calculated or rescaled.

The layout follows the [APBS PQR field specification](https://apbs.readthedocs.io/en/latest/formats/pqr.html):
whitespace-separated atom identity, optional chain, residue identity, x/y/z,
partial charge and radius. Coordinates and radii use Ångström, charges use
electrons. The output has no trailing element column and does not impose PDB's
fixed coordinate widths.

## APIs and supported records

- Rust: `write_pqr_checked(&[PqrAtomRecord])`,
  `write_pqr_checked_with_limits(&[PqrAtomRecord], &PqrWriteLimits)` and
  `PqrWriteError`. Existing `PqrError` variants are unchanged.
- Python: `write_pqr_checked(atoms, *, max_output_bytes=67_108_864)` accepts a
  bounded iterable of dictionaries.
- WASM: `write_pqr_json_checked(records_json)` accepts a strict JSON array.

```python
parsed = chematic.parse_pqr(text)
output = chematic.write_pqr_checked(parsed["atoms"])
reparsed = chematic.parse_pqr(output)
```

Atom records retain their order, ATOM/HETATM kind, signed 64-bit serial and
residue numbers, names, optional chain and insertion code. Names and chains
must be nonempty printable ASCII tokens of at most 64 bytes; insertion codes
are optional ASCII letters. Radius must be finite and nonnegative; coordinates
and charge must be finite. Empty input and unrepresentable records are refused
without partial output.

PQR has no authoritative element field. `infer_element` is a documented name
heuristic, not chemical perception. The checked writer refuses a supplied
element that would change under this heuristic after reparse. It neither
calculates charges nor assigns formal charges/bonds, verifies charge models,
validates duplicate serials or infers chemistry from radii.

Python and WASM require `group_pdb`, `serial`, `atom_name`, `res_name`, `res_seq`,
`x`, `y`, `z`, `charge` and `radius`; missing values are not filled with zeros.
`chain_id`, `icode` and `element` are optional. Omitted/null chain and insertion
code mean absence; omitted/null element invokes the same inference heuristic.
Unknown fields and boolean numeric values are refused. WASM also refuses
duplicate JSON fields, including escaped-key aliases and repeated null fields.
Python receives dictionaries, which cannot represent duplicate keys.

This API writes the supplied **atom fields**, not an original document. Source
comments, other record types, layout, numeric spelling and metadata outside
`PqrAtomRecord` are not represented. PQR parsing has its existing atom-projection
contract. No new common-format route or CLI command is added.

The legacy Rust/Python `write_pqr` and WASM `write_pqr_json` remain available
with their existing defaults and output: coordinates round to three decimal
places, charge/radius to four. Use the checked APIs when that rounding would
alter supplied parameters.

## Resource bounds

| Surface | Atom limit | Input limit | Output limit |
| --- | ---: | --- | ---: |
| Rust default | 100,000 | supplied slice | 64 MiB |
| Python | 100,000 | bounded iterable | 1 byte–64 MiB, configurable |
| WASM | 10,000 | 1,000,000 UTF-8 JSON bytes | 16 MiB |

Rust callers can configure `PqrWriteLimits`. Output accounting includes spaces,
newlines and all atom fields. Python stops at the 100,001st row of a longer
iterable, and invalid output options are rejected before consuming geometry.
WASM checks JSON byte length before decoding. These are operation budgets,
not a process memory cap. The JSON byte cap may bind before its atom cap.

For WASM, the contract concerns the values received in `records_json`.
Re-stringifying parsed JS numbers changes negative zero and can round integer
IDs beyond JS's safe integer range. Tests pass raw JSON when those distinctions
matter. The PQR output itself retains received signed 64-bit IDs and f64 bits.

## Development measurements

Original fixtures live in `validation/pqr-checked-2026-10-22/cases.json`.
The paired writer lane supplies the same atom arrays to checked and legacy
APIs, then uses an independent Python field observer. Neither engine is the
oracle. All five numerical fields and ordered identity metadata are checked.

| Writer lane, 12 common atom arrays | Checked Python/WASM | Legacy Python/WASM |
| --- | ---: | ---: |
| All numeric fields bit-identical | 12/12 | 3/12 |
| Ordered identity metadata retained | 12/12 | 12/12 |
| Maximum coordinate absolute error (Å) | 0 | 0.00045678901234680325 |
| Maximum charge absolute error (electrons) | 0 | 0.000043210987654335176 |
| Maximum radius absolute error (Å) | 0 | 0.000045678901234635916 |

The legacy control is the old API in the same candidate runtime. Its native
function body is verified equal to main `e9e578e2`; this is not a separately
rebuilt main executable or a published-release result. Eight synthetic numeric
boundary controls and six metadata controls are separate candidate lanes, not
additional Open Babel failures. Native properties run 4,096 deterministic bit
pattern trials (skipping nonfinite samples) plus explicit numerical boundaries.

### Open Babel read/write pipeline lanes

Local Open Babel 3.2.1 runs `-ipqr -opqr -ab` on original authored PQR inputs.
These are **read/write pipeline** observations; they do not isolate writer
rounding or hold identical internal `OBMol` state.

1. Twelve APBS-style whitespace inputs omit a trailing element column. Their
   source fields match authored expectations before execution. Open Babel's
   resulting atom records are unreadable by the independent APBS-field observer
   in all twelve cases. Raw stdout/stderr are archived; unreadability is not
   scored as a per-field loss, safe refusal or invalid-input ranking.
2. Four separate controls supply fixed coordinate slots and the trailing element
   field used by Open Babel's input dialect. These inputs also independently
   match the authored values. Open Babel emits readable records, but 0/4 retain
   all five numeric fields bit for bit. The checked Python/WASM APIs retain 4/4
   when supplied the same authored atom values; the legacy writer retains 0/4.

The dialect distinction is grounded in the pinned
[Open Babel 3.2.1 PQR implementation](https://github.com/openbabel/openbabel/blob/openbabel-3-2-1/src/formats/pqrformat.cpp).
The current chematic PQR reader does not accept a trailing element extension;
the extension controls are not claimed as chematic file-input coverage. No
toolkit-wide superiority, APBS execution or independently held-out chemistry
claim follows from these development controls.

Python passes 26 positive record controls and 40 refusal controls; actual
Node/WASM passes the same positives/refusals plus three duplicate-key JSON
controls. Additional API tests cover nonfinite fields, bounded iterators,
100,000/100,001 atoms, exact output bytes and the WASM UTF-8 input cap.
Six independent observer tests run in CI and reject altered fields, signed-zero
normalization, metadata changes, truncated/extra records, nonfinite output and
contradictory success/refusal envelopes.

## Evidence and reproduction

`validation/results/pqr-checked-2026-10-22.json` and its `.packet.json.gz`
record source/fixture/runtime hashes, worker outputs, raw Open Babel outputs,
per-field observations and process limits. Subprocesses are capped at 30 seconds
and 1 MiB combined stdout/stderr; no process memory cap is imposed.

```sh
cargo test -p chematic-mol --test pqr_checked
python3 -m pytest crates/chematic-py/tests/test_pqr_checked.py
node crates/chematic-wasm/tests/pqr_checked.test.mjs
python3 -m unittest discover -s scripts/tests -p test_pqr_checked_validation.py
python3 scripts/validation/validate_pqr_checked.py \
  --python-root /path/to/built-package-root \
  --wasm-module /path/to/chematic_wasm.js
```
