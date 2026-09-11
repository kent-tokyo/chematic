# Official RDKit.js comparison — v1.0.12

This is a same-process Node/WASM comparison between the regenerated v1.0.12
`crates/chematic-wasm/pkg-node` artifact and
`@rdkit/rdkit@2025.3.4-1.0.0` on macOS arm64. It uses 1,000 rows from
`scripts/descriptor_census_corpus.smi`, 20 warm-up rows, and Morgan radius 2
with 2,048 bits. The machine-readable result is
[`competitive-benchmark-rdkitjs-2026-09-11-v1.0.12.json`](../validation/results/competitive-benchmark-rdkitjs-2026-09-11-v1.0.12.json).

| Metric | chematic | official RDKit.js |
|---|---:|---:|
| WASM raw | 3,934,713 bytes | 6,914,823 bytes |
| WASM gzip-9 | 1,434,786 bytes | 2,045,412 bytes |
| Initialization | 0.043 ms | 41.545917 ms |
| Parse p50 | 0.003375 ms/mol | 0.233792 ms/mol |
| SMILES write p50 | 0.018542 ms/mol | 0.321708 ms/mol |
| ECFP4/Morgan p50 | 0.225437 ms/mol | 0.391333 ms/mol |
| Peak RSS | 109,772,800 bytes | 145,866,752 bytes |
| Exact fingerprint matches | 1,000 / 1,000 | 1,000 / 1,000 |

Reproduction:

```sh
npm install --prefix /private/tmp/chematic-rdkitjs-v1.0.12 \
  --ignore-scripts --no-package-lock @rdkit/rdkit@2025.3.4-1.0.0
node scripts/bench_wasm_vs_rdkit.mjs \
  --rdkit-package /private/tmp/chematic-rdkitjs-v1.0.12/node_modules/@rdkit/rdkit \
  --schematic-dir crates/chematic-wasm/pkg-node \
  --output /private/tmp/chematic-rdkitjs-v1.0.12/compare.json \
  --rows 1000 --warmup 20 --same-process
```

This is a Node same-process result, not browser-engine evidence. It does not
establish universal chemistry compatibility or superiority outside the named
corpus, operations, runtime, and hardware.
