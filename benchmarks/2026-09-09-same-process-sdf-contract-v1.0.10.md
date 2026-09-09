# v1.0.10 same-process SDF semantic contract

The schematic Python extension and RDKit were run in the same CPython process
against the identical `streaming.sdf` bytes. Across 20 repetitions each engine
returned 40 records with zero failures, and canonical SMILES matched in every
repetition. Two malformed SDF controls produced zero accepted records and one
failure in both engines. The schematic side uses the file-backed batched API,
whose progress manifest exposes rejected records, so the comparison does not
silently discard failure information.

The recorded elapsed times are context only. They are not a speed ranking:
the APIs use different SDF readers and the harness is a semantic contract
gate.

Reproduction:

```text
PYTHONPATH=crates/chematic-py/python:/Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages \
  python3 scripts/check_same_process_sdf_contract.py \
  --output benchmarks/2026-09-09-same-process-sdf-contract-v1.0.10.json
```

Machine-readable evidence: [`2026-09-09-same-process-sdf-contract-v1.0.10.json`](2026-09-09-same-process-sdf-contract-v1.0.10.json).
