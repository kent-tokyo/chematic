# v1.0.10 same-process MOL2 semantic contract

The ethanol MOL2 fixture was parsed in one CPython process by schematic's
`from_mol2` and RDKit's `MolFromMol2Block`. Across 20 repetitions, both
engines returned 20 records with zero failures and zero heavy-atom/ring-count
signature mismatches. Two malformed controls were rejected once by each
engine.

MOL2 charge inference and canonical SMILES are deliberately not compared:
the fixture leaves charge interpretation to the reader and the engines report
different but structurally bounded representations. Elapsed times are context
only and are not a speed ranking.

Reproduction:

```text
PYTHONPATH=crates/chematic-py/python:/Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages \
  python3 scripts/check_same_process_mol2_contract.py \
  --output benchmarks/2026-09-09-same-process-mol2-contract-v1.0.10.json
```

Machine-readable evidence: [`2026-09-09-same-process-mol2-contract-v1.0.10.json`](2026-09-09-same-process-mol2-contract-v1.0.10.json).
