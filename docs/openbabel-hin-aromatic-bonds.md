# HIN aromatic bonds and immutable source records

The new bounded HyperChem HIN reader addresses the input side of
[Open Babel #3000](https://github.com/openbabel/openbabel/issues/3000), checked
open on 2026-10-15. On an original synthetic six-carbon ring, Open Babel 3.2.1
writes all six HIN `a` bonds as MDL type **5** (single-or-double query).
chematic maps `a` to `BondOrder::Aromatic` and writes MDL type **4**.
This is a bond-representation check, not a general chemistry or speed claim.

## Supported profile

- Exactly one `mol <positive-id> "title"` / `endmol <same-id>` record.
- `atom id name element type flags charge x y z degree [neighbor bond]...`.
- Reciprocal, matching `s`, `d`, `t`, `a` connectivity. Source IDs are mapped
  explicitly to zero-based atom-row indices; nonconsecutive IDs are accepted
  by the local profile and tested separately from the paired reproduction.
- Finite partial charges and Cartesian coordinates in Å; names/type/flag
  tokens are retained independently of the graph. Partial charges never become
  formal charges. Numeric spelling, comments and CRLF/LF survive exact-source
  retention; the typed view preserves parsed f64 values.
- Default limits: 1 MiB UTF-8 input, 4,096 atoms, 16,384 unique bonds.
  Rust callers can provide `HinParseLimits`. WASM additionally retains its
  shared 1,000,000-byte input and 16 MiB JSON-output limits.

Empty/multiple molecules, dangling/self/duplicate references, mismatched
reciprocity, unknown bond tokens, nonfinite numbers, malformed counts and
unsupported directives reject explicitly. The accepted comment syntax is a
whole line beginning with `;`. Residue/system/force-field directives and
inline annotations are outside this profile.

The graph projection uses neutral formal-charge defaults and the core's
ordinary implicit-H query rules. It does not interpret atom types/flags,
assign stereo from coordinates, add atoms or attach coordinates/partial
charges to `Molecule`. The source record retains these fields. There is no
writer for an edited graph and no new generic conversion-bridge/CLI route.

## APIs

Rust: `parse_hin_record`, `parse_hin_record_with_limits`, immutable `HinRecord`
accessors, `HinRecord::to_molecule`, `HinRecord::to_json`, `write_hin_record`.
Existing MOL/conformer writers can consume the projected graph and separately
retained coordinates; cross-format metadata preservation is not implied.

```python
import json
import chematic

record = json.loads(chematic.parse_hin_record_json(text))
assert chematic.roundtrip_hin_record(text) == text
mol = chematic.from_hin(text)
block = mol.to_mol_block()  # graph output; source coordinates/charges remain in record
```

WASM exposes `parse_hin_record_json` and `roundtrip_hin_record` with the same
`chematic.hin-record.v1` JSON schema. Both bindings call the same Rust reader.

## Validation and reproduction

`validation/openbabel_issues/hin/aromatic-ring.hin` is original synthetic
project data, not an upstream attachment. The manifest records its license,
provenance, SHA-256 and expected literal bond semantics. The paired fixture
uses sequential IDs so that aromatic bond handling is isolated from source-ID
handling differences.

```sh
cargo test -p chematic-mol --lib --test hin_contract --offline --quiet
python3 -m pytest -q crates/chematic-py/tests/test_hin.py
node crates/chematic-wasm/tests/hin.test.mjs
python3 scripts/validation/validate_openbabel_hin.py
```

Python needs the extension built from this source. Node defaults to `pkg-node`
and accepts `CHEMATIC_WASM_MODULE` for an explicitly built module. The paired
runner records source/binary/runtime hashes and raw outputs under
`validation/results/openbabel-hin-aromatic-2026-10-15.*`.

Local verification passes 769 existing MOL unit tests and 11 HIN contract tests,
56 Python tests (16 HIN plus existing MOL2/CML checks), and the Node/WASM HIN
contract. The same-input run observes `[4,4,4,4,4,4]` for chematic and
`[5,5,5,5,5,5]` for Open Babel 3.2.1. The extension and WASM modules are local
source builds, not published packages. The first Python build and a later Rust
build encountered disk exhaustion; task-owned unused caches were removed and
the checks rerun. No upstream Open Babel code or issue was changed.
