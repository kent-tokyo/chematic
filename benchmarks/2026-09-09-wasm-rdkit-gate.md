# Node/WASM comparison gate: schematic vs RDKit.js

Measured 2026-09-09 on macOS arm64 with Node.js v24.5.0. This is a
same-corpus implementation comparison between the generated v1.0.10 schematic
WASM artifact and the installed `@rdkit/rdkit` npm package exposing RDKit
2025.03.4. It is not a comparison against the not-yet-measured future release
described by the RDKit packaging work.

Machine-readable output: [`2026-09-09-wasm-rdkit-gate.json`](2026-09-09-wasm-rdkit-gate.json).

## Protocol

- 1,000 rows from `scripts/descriptor_census_corpus.smi`, SHA-256
  `d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`;
- 20 warm-up rows, then 1,000 timed rows per isolated Node.js process;
- parse and write lanes parse each SMILES inside the timed operation;
- fingerprint lane uses schematic `rdkit_ecfp4_bitvec` and RDKit.js Morgan
  radius 2, 2,048-bit binary output;
- peak RSS is the per-process `resourceUsage().maxRSS`, not a browser-tab
  memory measurement.

## Results

| Lane | schematic p50 | RDKit.js p50 |
|---|---:|---:|
| Initialization | 7.827 ms | 81.025 ms |
| Parse | 0.0071 ms | 0.5608 ms |
| SMILES write | 0.0316 ms | 0.7689 ms |
| RDKit-compatible ECFP4/Morgan | 0.5314 ms | 0.9760 ms |

The reported run used 98.3 MB peak RSS for schematic and 133.1 MB for RDKit.js.
The WASM artifact was 3,726,056 bytes raw / 1,357,406 bytes gzip for schematic,
versus 6,914,823 bytes raw / 2,045,412 bytes gzip for RDKit MinimalLib.

Fingerprint parity was **1,000/1,000 exact matches** for the named
RDKit-compatible operation. This is a representation check, not a general
chemistry-accuracy claim.

## Reproduction

```sh
npm install --prefix /tmp/chematic-rdkit-gate --ignore-scripts --no-package-lock \
  @rdkit/rdkit@2025.3.4-1.0.0
node scripts/bench_wasm_vs_rdkit.mjs \
  --rdkit-package /tmp/chematic-rdkit-gate/node_modules/@rdkit/rdkit \
  --rows 1000 --warmup 20 \
  --output /tmp/chematic-wasm-rdkit-gate.json
```

The gate intentionally keeps native Python/RDKit measurements separate from
these Node/WASM measurements. A larger corpus, browser-engine runs, and the
new official release candidate remain follow-up gates.
