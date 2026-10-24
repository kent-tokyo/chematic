# V2000 property precedence

`cases.json` is an original synthetic MIT packet: five controls and nine
reproductions with a spec-authored oracle. No RDKit/Open Babel output supplies
the expected fields. Sparse `M  CHG` or `M  ISO` supersedes the complete old
atom-block channel, including atoms not listed in the properties.

See [the contract, compatibility difference and limits](../../docs/reference/v2000-property-precedence.md).
The baseline is `342d1d5f`; its affected reader is identical to the candidate's
parent `506a5219`. Keep the baseline Python/Node WASM artifacts separate from
candidate rebuilds; use the Cargo.lock-pinned `wasm-bindgen` (0.2.129 here).

With RDKit 2025.09.3, Node and Open Babel 3.2.1 at `/opt/homebrew/bin/obabel`:

```sh
PYTHONPATH=/path/to/baseline-python python3 scripts/validation/check_v2000_property_precedence.py \
  --wasm-module /path/to/baseline-wasm/chematic_wasm.js \
  --baseline --runtime-revision 342d1d5f9ad09aa6ade4e3079a902cba856209e2 \
  --output validation/results/v2000-property-precedence-baseline-2026-10-24.json.gz
PYTHONPATH=/path/to/candidate-python python3 scripts/validation/check_v2000_property_precedence.py \
  --wasm-module /path/to/candidate-wasm/chematic_wasm.js \
  --output validation/results/v2000-property-precedence-candidate-2026-10-24.json.gz
python3 -m unittest scripts.tests.test_v2000_property_precedence -v
```

Reports retain inputs, outputs, expected/observed fields, process outcomes and
source/runtime/tool hashes. Candidate `head` and `runtime_revision=working_tree`
describe a pre-commit source measurement; verify its reader hash against Git.
External subprocess limits are 10 seconds and 1 MiB of combined output.
Empty, failed or unreadable output cannot pass.

SMILES observation matches the packet's uniquely identified elements, charges,
isotopes, maps and bonds while allowing atom-order changes. Repeated elements
are explicitly outside that observer. RDKit's direct-MOL source diagnostic is
separate because its isotope precedence can differ from the MDL oracle. Open
Babel uses documented name/atom-class SMILES options and passes all 14 cases.
No H/stereo/geometry/source-byte/performance or general chemistry-validity claim
is made. Native-only empty property controls are legacy compatibility tests,
not scored format-validity successes.
