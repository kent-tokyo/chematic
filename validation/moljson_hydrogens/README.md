# Local atom-record MolJSON hydrogen authority

`cases.json` is original MIT synthetic data: 13 achiral inputs and 9 invalid
present H values. The expected non-node H arrays and reference SMILES are
authored independently of the reader. Zero H, an omitted field, partial/full
explicit H nodes, charges, isotopes and aromatic N-H have distinct controls.
The u8 storage endpoints (0/255) are tested separately in Rust; 255 is not a
chemical-validity or external-reader claim.

Build the Python and Node/WASM bindings from the current checkout, then run
with Python able to import both that fresh chematic package and RDKit:

```sh
PYTHONPATH=/path/to/fresh/python/package \
python scripts/validation/check_moljson_hydrogens.py \
  --wasm-module /path/to/fresh/pkg-node/chematic_wasm.js \
  --output /tmp/moljson-h-candidate.json.gz
cargo test --locked -p chematic-mol --lib --test moljson_hydrogens
cargo test --locked -p chematic-cli
python -m pytest crates/chematic-py/tests/test_moljson_hydrogens.py
node crates/chematic-wasm/tests/moljson_hydrogens.test.mjs /path/to/fresh/pkg-node/chematic_wasm.js
python -m unittest scripts.tests.test_moljson_hydrogens
```

The paired runs use identical packet and observation code. The older runtime
is the saved Python/Node build from `cfbce7d66a096f0a44cf46541bc239af989cad75`
(PR #1190); its production `moljson.rs` is identical to this change's main
parent `aeae3f53deb34579be60c054cdf22cd686d9254d`. This is a previous local
runtime, not an installed release. Pass `--baseline` and
`--baseline-runtime-revision <SHA>` for that run. The revision is operator
provenance; the runner records runtime hashes rather than proving arbitrary
binaries were built from a supplied SHA.

Reports are saved as `../results/moljson-hydrogens-{baseline,candidate}-2026-10-24.json.gz`.
They record source, packet, harness, extension, WASM, RDKit parser and Open Babel
binary hashes, runtime versions, source inputs, candidate outputs, Node commands
and results, individual observations and refusal vectors. Gzip timestamps are
zero. CLI subprocesses use 10-second/1-MiB combined-output limits. RDKit observes
SMILES with sanitization on and H removal off. The observer compares canonical
SMILES, formula, node count and radical-electron count; it does not independently
parse the local MolJSON dialect. Six unit controls exercise damaged output,
count type/shape, graph changes and unreadable graph handling through the actual
observation function; the diagnostic uses the real RDKit reader.

The baseline retains 5/13 graphs in each binding and refuses 0/9 invalid H
values. Candidate: 13/13 and 9/9. Repeated round trips are regression tests, not
extra corpus samples. Open Babel has no advertised MolJSON format in the saved
inventory and is not assigned an artificial loss score. The current oxpig
schema is a different dialect, as documented in the reference page.
