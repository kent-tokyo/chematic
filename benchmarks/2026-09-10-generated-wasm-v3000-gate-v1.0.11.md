# Generated WASM V3000 metadata gate — 2026-09-10

The source-built v1.0.11 Web artifact was generated with `wasm-pack 0.13.1`
and exercised in Playwright Chromium `HeadlessChrome/153.0.8010.12` against
the pinned `@rdkit/rdkit@2025.3.4-1.0.0` package. The same 1,000-row corpus
and timing lane as the official browser gate were used.

The generated artifact exported `roundtrip_mol_v3000_block` and passed the
V3000 metadata probe: an opaque `SUP` SGROUP line survived, and output placed
`SGROUP` before `COLLECTION`. Fingerprint parity remained 1,000/1,000. The
generated-artifact machine result, including all artifact and package digests,
is [here](../validation/results/competitive-browser-rdkitjs-generated-2026-09-10-v1.0.11.json).

This is generated-artifact evidence for the local Web build, not a published
package or cross-platform result. The checked-in `demo/pkg` artifact remains
separate until the release workflow intentionally regenerates it.

Reproduction:

```sh
wasm-pack build crates/chematic-wasm --target web \
  --out-dir /private/tmp/chematic-wasm-v1.0.11-browser-gate --mode no-install
node scripts/bench_browser_wasm_vs_rdkit.mjs \
  --rdkit-package /private/tmp/chematic-rdkitjs-gate/node_modules/@rdkit/rdkit \
  --output /private/tmp/ignored-browser-result.json \
  --html-output /private/tmp/chematic-browser-generated-gate-1000.html \
  --rows 1000 --warmup 20 --init-timeout-ms 10000
python3 scripts/bench_browser_wasm_vs_rdkit_playwright.py \
  --html /private/tmp/chematic-browser-generated-gate-1000.html \
  --output validation/results/competitive-browser-rdkitjs-generated-2026-09-10-v1.0.11.json \
  --schematic-dir /private/tmp/chematic-wasm-v1.0.11-browser-gate \
  --rdkit-package /private/tmp/chematic-rdkitjs-gate/node_modules/@rdkit/rdkit \
  --browser /Users/k_tanabe/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell \
  --corpus scripts/descriptor_census_corpus.smi --rows 1000 --timeout-ms 60000
```
