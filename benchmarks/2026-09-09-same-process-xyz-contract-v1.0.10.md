# v1.0.10 same-process XYZ semantic contract

The two XYZ frames in `streaming.xyz` were parsed in one CPython process by
schematic's `from_xyz` and RDKit's `MolFromXYZBlock`. Across 20 repetitions,
both engines returned 40 frames with zero failures and zero total/heavy-atom
count or coordinate-finiteness mismatches. Two malformed controls were
rejected once by each engine.

XYZ bond perception is deliberately not compared: the contract covers frame
and coordinate semantics only. Elapsed times are context and are not a speed
ranking.

Reproduction:

```text
PYTHONPATH=crates/chematic-py/python:/Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages \
  python3 scripts/check_same_process_xyz_contract.py \
  --output benchmarks/2026-09-09-same-process-xyz-contract-v1.0.10.json
```

Machine-readable evidence: [`2026-09-09-same-process-xyz-contract-v1.0.10.json`](2026-09-09-same-process-xyz-contract-v1.0.10.json).
