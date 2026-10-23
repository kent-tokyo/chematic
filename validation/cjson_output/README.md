# Checked ChemicalJSON output, 2026-10-24

This development record measures **CML → CJSON** for eight original synthetic
two-atom graph records with complete source Cartesian geometry. These are
representation tests, not a realistic chemistry corpus or a toolkit ranking.

The three positive records cover single, double and triple bonds. The five
negative records use chematic's CML reader profile for zero and query bonds:
`0`, `any`, `S/D`, `S/A`, `D/A`. The query tokens are local reader conventions;
this experiment does not establish that every CML implementation must support
them. The failure being addressed is chematic emitting a different bond with
no indication that its input model was changed.

| Path | Positive graph and XYZ retained | Unsupported records emitted as ordinary bonds | Explicit refusal |
| --- | ---: | ---: | ---: |
| Python checked common bridge | 3/3 | 0/5 | 5/5 |
| WASM checked common bridge | 3/3 | 0/5 | 5/5 |
| Python legacy `Mol.to_cjson` | 3/3 | 5/5 | 0/5 |
| Open Babel 3.2.1 default | 3/3 | 5/5 | 0/5 |
| Open Babel 3.2.1 `-xs` | 3/3 | 5/5 | 0/5 |

The legacy path is called on the **same candidate build** after parsing the
same CML input, to isolate the old writer's behavior. It is not a separately
built main/released baseline. Open Babel returned exit 0 and `1 molecule
converted` for all eight records in both modes. It emitted order 1 for `0`,
`any`, `S/D`, `S/A`, and order 2 for `D/A`; no output field retained the source
query/zero distinction. The complete stdout/stderr and outputs are recorded.
These observations establish a difference in refusal behavior for this profile,
not a claim of CML specification nonconformance or general superiority.

The independent JSON observer checks source-indexed elements, formal charges,
isotopes, bond endpoints/orders and binary64 coordinate bits. Expected values
are literal fixture data, not another toolkit's output. Missing output,
malformed output, changed output and refusal are counted separately. Processes
are bounded to 10 seconds and 1 MiB combined output. The report retains all
inputs, outputs, errors, and hashes/identities of runtime artifacts and sources:
[raw compressed report](../results/cjson-output-2026-10-24.json.gz).

Reproduce against fresh Python and Node WASM builds:

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_cjson_output.py \
  --wasm-module /path/to/chematic_wasm.js \
  --output /tmp/cjson-output.json.gz
```

Native regression tests additionally cover every coerced bond order, including
dative direction, quadruple/aromatic alias ambiguity and Up/Down stereo;
atom properties, enhanced stereo groups, stashed aromatic bond directions,
coordinate row counts/nonfinite values, self bonds outside the reader profile,
and exact supported graph/XYZ roundtrips.
Python and real Node WASM tests cover common-bridge refusals and supported
output, plus eight SMILES with unsupported properties. CLI integration tests
also require nonzero exit and empty stdout for those eight SMILES, and retain
the supported CCO topology. This does not add CLI source-coordinate retention.
Those SMILES tests are
regressions, not part of the paired Open Babel table.

The checked writer preserves the **already parsed in-memory profile**. It
cannot recover source metadata or stereo ignored earlier by a reader.
In particular, CML output from chematic explicitly declares hydrogen counts;
reading that CML and requesting CJSON now refuses to discard them. This is an
intentional compatibility change in the common bridges. The old lossy Rust
`write_cjson` and Python `Mol.to_cjson` remain available. Exact CJSON retention
and the separate MOL2 selected-field projection retain their existing contracts.
