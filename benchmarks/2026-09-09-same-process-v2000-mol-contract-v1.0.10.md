# v1.0.10 same-process V2000 MOL semantic contract

The two V2000 MOL blocks in `streaming.sdf` were split at the SDF delimiter
and parsed in one CPython process by schematic's `from_mol_block` and RDKit's
`MolFromMolBlock`. Across 20 repetitions, both engines returned 40 records
with zero failures and zero canonical-SMILES mismatches. Two malformed controls
were rejected once by each engine.

Elapsed times are context only and are not a speed ranking. The gate compares
materialized MOL-block semantics, not identical internal implementations.

Reproduction:

```text
PYTHONPATH=crates/chematic-py/python:/Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages \
  python3 scripts/check_same_process_mol_contract.py \
  --output benchmarks/2026-09-09-same-process-v2000-mol-contract-v1.0.10.json
```

Machine-readable evidence: [`2026-09-09-same-process-v2000-mol-contract-v1.0.10.json`](2026-09-09-same-process-v2000-mol-contract-v1.0.10.json).
