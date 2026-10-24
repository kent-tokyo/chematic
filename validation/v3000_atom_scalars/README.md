# V3000 atom scalars

`cases.json` contains original, synthetic MIT-licensed single-atom V3000 blocks:
8 ordinary scalar cases, 2 local quoted `LABEL` compatibility controls and 25
invalid/unsupported scalar cases. The authored four-field oracle is independent
of the runtime parsers. No external molecule corpus was copied.

See [the contract and measured limits](../../docs/reference/v3000-atom-scalar-integrity.md).
The fixed baseline is `342d1d5f9ad09aa6ade4e3079a902cba856209e2`.

Build separate Python and Node WASM runtimes for baseline and candidate; use the
`wasm-bindgen` version pinned by `Cargo.lock` (0.2.129 in this run). Keep the
runtimes in distinct directories and invalidate Cargo's shared target cache when
switching source trees. In this measurement the native and WASM regression
tests confirmed candidate interpretation before scoring it.

With RDKit 2025.09.3, Open Babel 3.2.1 at `/opt/homebrew/bin/obabel`, and Node:

```sh
PYTHONPATH=/path/to/baseline-python python3 scripts/validation/check_v3000_atom_scalars.py \
  --wasm-module /path/to/baseline-wasm/chematic_wasm.js \
  --baseline --runtime-revision 342d1d5f9ad09aa6ade4e3079a902cba856209e2 \
  --output validation/results/v3000-atom-scalars-baseline-2026-10-24.json.gz
PYTHONPATH=/path/to/candidate-python python3 scripts/validation/check_v3000_atom_scalars.py \
  --wasm-module /path/to/candidate-wasm/chematic_wasm.js \
  --output validation/results/v3000-atom-scalars-candidate-2026-10-24.json.gz
python3 -m unittest scripts.tests.test_v3000_atom_scalars -v
```

The reports retain inputs, expected fields, outputs, terminal outcomes and hashes.
Candidate `runtime_revision=working_tree` and `head` identify a pre-commit source
measurement; verify `production_sources` hashes against the committed files.
The reader observes emitted SMILES, not the original V3000 dialect. Refusal,
exact source retention, ordinary field preservation and quoted compatibility
are separate observations. An empty output, signal or timeout never counts as
preservation. No speed or broad-format claim follows from this packet.
