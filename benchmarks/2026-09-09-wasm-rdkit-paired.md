# Paired Node/WASM comparison: schematic vs RDKit.js

This is the same-corpus paired follow-up to the isolated-process gate. Both
WASM implementations are loaded and measured in one Node.js v24.5.0 process
on macOS arm64, using 1,000 rows and 20 warm-up rows. The machine-readable
result is [`2026-09-09-wasm-rdkit-paired.json`](2026-09-09-wasm-rdkit-paired.json).

| Lane | schematic p50 | RDKit.js p50 |
|---|---:|---:|
| Parse (parse + free) | 0.0091 ms | 0.6277 ms |
| SMILES write (parse + write) | 0.0386 ms | 0.7519 ms |
| RDKit-compatible ECFP4 / Morgan (parse + fingerprint) | 0.5614 ms | 1.0035 ms |

The RDKit lane reports RDKit 2025.03.4. The RDKit-compatible fingerprint
matched exactly for **1,000/1,000** rows. The fingerprint definitions and
bit-packing conversion are explicit in the JSON record.

This paired lane removes the separate-process startup boundary from the
timed operations, but it does not make the implementations chemically
identical. Peak RSS is intentionally not interpreted per lane here because
the two runtimes share one process; use the isolated-process report for that
axis. Browser-engine and future RDKit package measurements remain separate
follow-up gates.

Reproduce with:

```sh
node scripts/bench_wasm_vs_rdkit.mjs --same-process \
  --rdkit-package /tmp/chematic-rdkit-gate/node_modules/@rdkit/rdkit \
  --rows 1000 --warmup 20 \
  --output /tmp/chematic-wasm-rdkit-paired.json
```
