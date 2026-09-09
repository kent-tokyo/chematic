# Same-process cross-engine contract bundle (v1.0.10)

This is a local semantic-contract rerun against the current v1.0.10 Python
extension and the installed RDKit 2025.09.3. It is not a throughput ranking:
the timing fields in the per-format JSON files are retained only as context.

## Result

Eight format-specific runners completed with 20 repetitions each. SDF,
V2000 MOL, V3000 MOL, MOL2, XYZ, Extended XYZ, and CDXML had zero valid-record,
failure-count, and applicable structural-signature mismatches. PDB had zero
valid-input signature mismatches. PDB and CDXML malformed controls are expected
boundary differences: the default schematic readers are intentionally lenient,
while RDKit rejects those controls; the mismatches are recorded rather than
normalized away.

The machine-readable reports are:

- [`2026-09-10-same-process-sdf-contract-v1.0.10.json`](2026-09-10-same-process-sdf-contract-v1.0.10.json)
- [`2026-09-10-same-process-v2000-mol-contract-v1.0.10.json`](2026-09-10-same-process-v2000-mol-contract-v1.0.10.json)
- [`2026-09-10-same-process-v3000-mol-contract-v1.0.10.json`](2026-09-10-same-process-v3000-mol-contract-v1.0.10.json)
- [`2026-09-10-same-process-mol2-contract-v1.0.10.json`](2026-09-10-same-process-mol2-contract-v1.0.10.json)
- [`2026-09-10-same-process-xyz-contract-v1.0.10.json`](2026-09-10-same-process-xyz-contract-v1.0.10.json)
- [`2026-09-10-same-process-extxyz-contract-v1.0.10.json`](2026-09-10-same-process-extxyz-contract-v1.0.10.json)
- [`2026-09-10-same-process-pdb-contract-v1.0.10.json`](2026-09-10-same-process-pdb-contract-v1.0.10.json)
- [`2026-09-10-same-process-cdxml-contract-v1.0.10.json`](2026-09-10-same-process-cdxml-contract-v1.0.10.json)

## Reproduction

The source-extension location is intentionally supplied through `PYTHONPATH`
so the normal user-site installation is not modified:

```sh
PYTHONPATH=/private/tmp/chematic-py-target-reaction-map-v3 \
  python3 scripts/check_same_process_sdf_contract.py --repeats 20 \
  --output benchmarks/2026-09-10-same-process-sdf-contract-v1.0.10.json
```

Repeat the same command with the seven other format-specific runners and output
paths listed above. The checked-in bundle gate then verifies all eight reports:

```sh
python3 scripts/check_same_process_contracts.py
```

The gate requires current-version metadata, malformed-case reporting, and an
explicit non-ranking timing boundary. It does not claim Extended XYZ property
parity because RDKit's comparison lane intentionally projects to coordinates.
Broader malformed-corpus coverage, generated WASM artifact parity, and
throughput-equivalent benchmarking remain open roadmap items.
