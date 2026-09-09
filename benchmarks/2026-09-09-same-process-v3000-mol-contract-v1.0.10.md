# v1.0.10 same-process V3000 MOL semantic contract

The V3000 ethanol fixture was parsed in one CPython process by schematic's
`from_mol_v3000` and RDKit's `MolFromMolBlock`. Across 20 repetitions, both
engines returned 20 records with zero failures and zero heavy-atom/ring-count
signature mismatches. Two malformed V3000 controls were rejected once by each
engine.

Elapsed times are context only and are not a speed ranking. The gate compares
materialized V3000 MOL semantics and keeps parser/process boundaries explicit.

Reproduction:

```text
PYTHONPATH=crates/chematic-py/python:/Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages \
  python3 scripts/check_same_process_v3000_contract.py \
  --output benchmarks/2026-09-09-same-process-v3000-mol-contract-v1.0.10.json
```

Machine-readable evidence: [`2026-09-09-same-process-v3000-mol-contract-v1.0.10.json`](2026-09-09-same-process-v3000-mol-contract-v1.0.10.json).
