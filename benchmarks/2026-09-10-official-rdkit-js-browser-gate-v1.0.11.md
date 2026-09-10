# Official RDKit.js browser gate — 2026-09-10

This is a real Playwright Chromium comparison of the v1.0.11 browser artifact
and `@rdkit/rdkit@2025.3.4-1.0.0`. It uses the first 1,000 rows of
`scripts/descriptor_census_corpus.smi` (SHA-256
`d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`), 20
warm-up rows, Morgan radius 2, and 2,048 bits.

| Metric | chematic | official RDKit.js |
|---|---:|---:|
| WASM raw | 3,726,056 bytes | 6,914,823 bytes |
| WASM gzip-9 | 1,357,406 bytes | 2,045,412 bytes |
| Initialization | 66.900 ms | 163.800 ms |
| Parse p50 | below 0.1 ms/mol | 0.2 ms/mol |
| SMILES write p50 | below 0.1 ms/mol | 0.3 ms/mol |
| ECFP4/Morgan p50 | 0.2 ms/mol | 0.4 ms/mol |
| Exact fingerprint matches | 1,000 / 1,000 | 1,000 / 1,000 |

The browser timing resolution rounds very small values to zero in the raw
report; the mean values remain available in the machine-readable result. This
is browser-engine evidence only and must not be mixed with the Node process
peak-RSS lane. The artifact and package digests are retained in
`validation/results/competitive-browser-rdkitjs-2026-09-10-v1.0.11.json`.

Reproduction:

```sh
node scripts/bench_wasm_vs_rdkit.mjs \
  --rdkit-package /private/tmp/chematic-rdkitjs-gate/node_modules/@rdkit/rdkit \
  --output /private/tmp/browser-html-generation-placeholder.json \
  --html-output /private/tmp/chematic-browser-gate-1000.html \
  --rows 1000 --warmup 20 --init-timeout-ms 10000
python3 scripts/bench_browser_wasm_vs_rdkit_playwright.py \
  --html /private/tmp/chematic-browser-gate-1000.html \
  --output validation/results/competitive-browser-rdkitjs-2026-09-10-v1.0.11.json \
  --schematic-dir demo/pkg \
  --rdkit-package /private/tmp/chematic-rdkitjs-gate/node_modules/@rdkit/rdkit \
  --browser /Users/k_tanabe/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell \
  --corpus scripts/descriptor_census_corpus.smi --rows 1000 --timeout-ms 60000
```
