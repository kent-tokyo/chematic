# v1.0.10 same-process Extended XYZ semantic contract

The two Extended XYZ frames in `streaming.extxyz` were parsed in one CPython
process. schematic used `from_extxyz`; RDKit received a coordinate-only XYZ
projection because it has no Extended XYZ reader in this lane. Across 20
repetitions, both sides returned 40 frames with zero failures and zero
frame-signature mismatches for atom count, heavy-atom count, and finite
coordinates. Two malformed coordinate/frame controls were rejected once by
both.

Extended XYZ metadata and per-atom properties are checked only on the
schematic side; property parity is not claimed. Elapsed times are context only
and are not a speed ranking.

Reproduction:

```text
PYTHONPATH=crates/chematic-py/python:/Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages \
  python3 scripts/check_same_process_extxyz_contract.py \
  --output benchmarks/2026-09-09-same-process-extxyz-contract-v1.0.10.json
```

Machine-readable evidence: [`2026-09-09-same-process-extxyz-contract-v1.0.10.json`](2026-09-09-same-process-extxyz-contract-v1.0.10.json).
