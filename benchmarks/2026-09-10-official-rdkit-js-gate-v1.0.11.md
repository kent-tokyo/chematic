# Official RDKit.js comparison gate — 2026-09-10

This is a same-process Node/WASM comparison between the v1.0.11 workspace
artifact and `@rdkit/rdkit@2025.3.4-1.0.0` (RDKit core version 2025.03.4). It uses the first 1,000 rows of
`scripts/descriptor_census_corpus.smi` (SHA-256
`d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`), 20
warm-up rows, Morgan radius 2, and 2,048 bits.

| Metric | chematic | official RDKit.js |
|---|---:|---:|
| WASM raw | 3,726,056 bytes | 6,914,823 bytes |
| WASM gzip-9 | 1,357,406 bytes | 2,045,412 bytes |
| Initialization | 3.767 ms | 31.762 ms |
| Parse p50 | 0.003167 ms/mol | 0.237437 ms/mol |
| SMILES write p50 | 0.017312 ms/mol | 0.320313 ms/mol |
| ECFP4/Morgan p50 | 0.228542 ms/mol | 0.401084 ms/mol |
| Peak RSS | 107,151,360 bytes | 142,557,184 bytes |
| Exact fingerprint matches | 1,000 / 1,000 | 1,000 / 1,000 |

The raw artifact SHA-256 values are:

- chematic WASM: `26e78553c5dd929b22b130bc1f0a700689256876099b13d4df63c1a5edb27468`
- RDKit WASM: `e0967d44fed59e44a2d07bfc08f3a86821c54a14c95bb7e6e392fe867333b8c8`
- RDKit JS: `58d3c996ade7e0b0137d4f9363ff6c204b545689428fa0898eeaed579ef788d9`
- RDKit package metadata: `c9029cad92384253cd6f885818ec5a357c045fdbda9ec196dc772e7db6df98fb`
- RDKit TypeScript declarations: `2744d98da7f92fdd77d7dbae6434a9346ad381e5c5c41096030311d01b2882be`

This is a Node same-process result, not browser-engine evidence. It does not
establish universal chemistry compatibility, and it does not compare APIs that
are not present in both lanes. The machine-readable result and reproduction
command are in
`validation/results/competitive-benchmark-rdkitjs-2026-09-10-v1.0.11.json`:

```sh
node scripts/bench_wasm_vs_rdkit.mjs \
  --rdkit-package /private/tmp/chematic-rdkitjs-gate/node_modules/@rdkit/rdkit \
  --output validation/results/competitive-benchmark-rdkitjs-2026-09-10-v1.0.11.json \
  --rows 1000 --warmup 20 --same-process
```
