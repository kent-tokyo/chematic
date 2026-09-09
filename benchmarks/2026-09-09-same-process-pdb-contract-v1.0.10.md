# v1.0.10 same-process PDB semantic contract

The minimal PDB fixture was parsed in one CPython process by schematic's
`from_pdb` and RDKit's `MolFromPDBBlock`. Across 20 repetitions, both engines
returned 20 records with zero failures and zero mismatches for total atoms,
heavy atoms, and finite coordinates.

The two malformed controls expose a deliberate parser-boundary difference:
schematic's PDB reader is lenient and returns a partial record, while RDKit
rejects them. This is recorded as an open malformed-parity gap; it is not
treated as a compatibility pass. Bond parity is also not claimed.

Reproduction:

```text
PYTHONPATH=crates/chematic-py/python:/Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages \
  python3 scripts/check_same_process_pdb_contract.py \
  --output benchmarks/2026-09-09-same-process-pdb-contract-v1.0.10.json
```

Machine-readable evidence: [`2026-09-09-same-process-pdb-contract-v1.0.10.json`](2026-09-09-same-process-pdb-contract-v1.0.10.json).
