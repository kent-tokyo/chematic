# Official RDKit.js isolated browser comparison — 2026-09-16

This is a fresh-process comparison in Chromium, Firefox, and WebKit of a locally
rebuilt current source artifact (workspace manifest version `1.0.15`) and the
public `@rdkit/rdkit@2026.03.6` package. It is **not** evidence about the already
published `v1.0.15` tag. It supersedes neither the older v1.0.11 browser gate nor the separate
fingerprint-compatibility evidence: it fixes their same-page measurement
limitation for this one environment and operation contract.

## Conditions

- Browser: Google Chrome on macOS arm64; each arm used a fresh Chromium process
  and cache-free local routes.
- Corpus: first 1,000 valid rows of `scripts/descriptor_census_corpus.smi`,
  SHA-256 `d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`.
- Repetitions: the initial 1k/5k records used three fresh processes per arm;
  the fixed 10k records below use 20. Arm order alternates and operation
  timings use 20 warm-up rows. Cold initialization has no warm-up.
- Operations: parse; parse + canonical-SMILES write; parse + radius-2,
  2,048-bit fingerprint. They are explicitly parse-inclusive, not prepared-object
  timings. All 1,000 input rows were accepted in both arms.

| Metric (p50 of per-process means unless stated) | chematic current source (manifest 1.0.15) | RDKit.js 2026.03.6 |
|---|---:|---:|
| WASM raw | 4,000,896 bytes | 7,333,095 bytes |
| WASM gzip-9 | 1,459,253 bytes | 2,379,975 bytes |
| Init p50 / p95 | 12.3 / 15.99 ms | 50.8 / 52.51 ms |
| Parse | 0.00690 ms/mol | 0.15460 ms/mol |
| Parse + canonical write | 0.09790 ms/mol | 0.18880 ms/mol |
| Parse + radius-2 FP | 0.31380 ms/mol | 0.17030 ms/mol |

The 1,000-row raw record is
`validation/results/competitive-browser-rdkitjs-isolated-v1.0.15-2026-09-16.json`.
It includes the package version, npm registry tarball URL and SRI integrity,
package/lock/artifact hashes, individual process results, and browser-exposed
JS-heap snapshots.

The same source-pinned corpus contains 5,000 usable rows, so a second run used
all 5,000 rather than duplicating entries to manufacture a 10,000-row result.
It kept the same three alternating fresh-process repetitions and 20-row warmup.

| Metric (p50 of per-process means unless stated) | chematic current source (manifest 1.0.15) | RDKit.js 2026.03.6 |
|---|---:|---:|
| Init p50 / p95 | 11.6 / 14.21 ms | 49.5 / 51.39 ms |
| Parse | 0.00436 ms/mol | 0.14614 ms/mol |
| Parse + canonical write | 0.09992 ms/mol | 0.17256 ms/mol |
| Parse + radius-2 FP | 0.24332 ms/mol | 0.15622 ms/mol |

The 5,000-row record is
`validation/results/competitive-browser-rdkitjs-isolated-5k-v1.0.15-2026-09-16.json`.
Both records use the same artifact hashes; timing differences between rows are
normal fresh-process variation, not a claim of a trend.

The companion direct bit-vector gate compares this same 5,000-row corpus in
one browser page: `rdkit_ecfp4_bitvec` versus RDKit.js Morgan (radius 2, 2,048
bits) is **5,000/5,000 exact**. Its raw record is
`validation/results/competitive-browser-rdkitjs-ecfp4-parity-5k-v1.0.15-2026-09-16.json`.
It covers that fingerprint configuration and corpus only.

The follow-up uses a **fixed 10,000-row common-parse corpus**, not duplicated
5k rows. `scripts/build_browser_rdkit_comparison_corpus.py` pins three exposed
source hashes into a 12,000-row candidate set; `scripts/qualify_browser_rdkit_common_corpus.py`
then retains the first 10,000 inputs accepted by both browser-resident engines.
The resulting corpus hash is
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
It is explicitly **not sealed** and must not be used for accuracy tuning or a
generalization claim.

| Metric (p50 of per-process means unless stated) | chematic current source (manifest 1.0.15) | RDKit.js 2026.03.6 |
|---|---:|---:|
| WASM raw / gzip-9 | 4,005,280 / 1,460,499 bytes | 7,333,095 / 2,379,975 bytes |
| Init p50 / p95 | 11.7 / 12.51 ms | 34.9 / 35.44 ms |
| Parse | 0.00222 ms/mol | 0.07218 ms/mol |
| Parse + canonical write | 0.05682 ms/mol | 0.10086 ms/mol |
| Parse + radius-2 FP (9,999 supported rows) | 0.145885 ms/mol | 0.094089 ms/mol |

The isolated raw record is
`validation/results/competitive-browser-rdkitjs-isolated-10k-v1.0.15-2026-09-16.json`.
The direct configured-bit gate is **9,999/9,999 exact on the declared supported
domain**. The former fused porphyrin-like aromatic residual was fixed by matching
RDKit's one-shared-bond, 24-bond ring-neighbour rule and is now a Rust/browser
regression. One Fe(II), degree-10 graph with two anionic carbon ligands remains
a typed `unsupported RDKit coordination sanitization` refusal; it remains in the
fixed corpus and is recorded rather than silently excluded. `OCl(=O)(=O)=O` was
a third residual until the RDKit-Morgan-only
hypervalent-halogen normalization was added; it is now a native and browser
regression. The parity raw record is
`validation/results/competitive-browser-rdkitjs-ecfp4-parity-10k-v1.0.15-2026-09-16.json`.

The Chromium arm was then extended from three to **20** fresh-process repetitions,
without changing the artifact hashes, corpus, operation contract, 20-row warm-up,
or one typed fingerprint refusal. The p50/p95 values below are Chromium-only;
Firefox and WebKit have their own independent 20-repetition records below and
are never pooled with it.

| Metric (p50 / p95 of per-process means) | chematic current source | RDKit.js 2026.03.6 |
|---|---:|---:|
| Init | 11.80 / 16.69 ms | 45.60 / 71.88 ms |
| Parse | 0.002670 / 0.003698 ms/mol | 0.087055 / 0.115387 ms/mol |
| Parse + canonical write | 0.068035 / 0.091582 ms/mol | 0.121025 / 0.153014 ms/mol |
| Parse + radius-2 FP (9,999 supported rows) | 0.171217 / 0.218389 ms/mol | 0.113076 / 0.142165 ms/mol |

The raw record is
`validation/results/competitive-browser-rdkitjs-isolated-10k-chromium-20rep-v1.0.15-2026-09-16.json`.
It strengthens the local timing estimate but still does not provide a
cross-host result, a download-plus-ready measurement, unique physical memory,
or a general fingerprint-performance claim.

## Published npm package confirmation (Chromium)

The preceding rows use a current-source rebuild and therefore do not prove the
contents of the already published npm package. A separate Chromium lane
downloaded `@kent-tokyo/chematic@1.0.15` from npm, retained the tarball SHA-256
`ddbe24150f6c4c130bab7c1e1791aea5e2407afda79d60fd432bd7cb07bb3314`, extracted
that package without rebuilding it, and repeated the same fixed 10k operation
against the pinned public RDKit.js package. Both arms have 20 fresh browser
processes, 10,000 accepted parse/write rows, and 9,999 fingerprint rows; the
same Fe(II) input remains an explicit typed fingerprint refusal.

| Metric (p50 / p95 of per-process means) | chematic npm 1.0.15 | RDKit.js 2026.03.6 |
|---|---:|---:|
| WASM raw / gzip-9 | 4,009,079 / 1,461,625 bytes | 7,333,095 / 2,379,975 bytes |
| Init | 13.95 / 17.71 ms | 58.40 / 68.28 ms |
| Parse | 0.003645 / 0.004018 ms/mol | 0.108980 / 0.122483 ms/mol |
| Parse + canonical write | 0.083500 / 0.096262 ms/mol | 0.148200 / 0.177089 ms/mol |
| Parse + radius-2 FP (9,999 supported rows) | 0.315347 / 0.366077 ms/mol | 0.134228 / 0.165473 ms/mol |

The raw record is
`validation/results/competitive-browser-rdkitjs-published-npm-10k-chromium-20rep-v1.0.15-2026-09-16.json`.
It records package metadata and hashes for both arms. The npm package does not
export `wasm_linear_memory_bytes`, so its linear-memory field is `unavailable`;
the matching RDKit.js field is also unavailable.

### Cold local download-to-ready (Chromium)

The published-package Chromium lane was rerun with a separately recorded
navigation-start-to-ready measure. Each arm uses a fresh browser process and
cache-disabled localhost routes; the timer includes HTML, JS glue/module, and
WASM asset requests plus initialization, and ends before the parse/write/FP
workload. It is therefore a reproducible **cold local asset** measure, not an
internet/CDN download claim or a general page-load benchmark.

| Metric (20 fresh processes) | chematic npm 1.0.15 | RDKit.js 2026.03.6 |
|---|---:|---:|
| Download-to-ready p50 / p95 | 42.40 / 51.62 ms | 81.20 / 98.495 ms |

The p95 ratio is 52.4% in this declared environment. This satisfies the
plan's 75% numerical target for the narrow local measure, but no paired
confidence interval, remote-network result, cross-host result, or universal
startup claim follows from it. The raw record is
`validation/results/competitive-browser-rdkitjs-published-npm-download-ready-10k-chromium-20rep-v1.0.15-2026-09-16.json`.

### Published-package process-tree RSS (Chromium)

The same tarball was also measured in three fresh Chromium CDP launches. Peak
RSS was sampled every 50 ms over each owned browser process tree; its median
was 2,189,262,848 bytes for chematic and 2,272,215,040 bytes for RDKit.js.
This is neither unique resident memory nor a memory-budget pass: the `ps` sum
can double-count shared pages, has only three repetitions, and does not cover
Firefox, WebKit, or another host. The raw record is
`validation/results/competitive-browser-rdkitjs-published-npm-10k-chromium-rss-v1.0.15-2026-09-16.json`.

The same exact npm tarball, 10k corpus, warm-up, and 20 fresh-process design
was also run in Firefox and WebKit. Browser results are not pooled.

| Browser / metric (p50 / p95 of per-process means) | chematic npm 1.0.15 | RDKit.js 2026.03.6 |
|---|---:|---:|
| Firefox init | 172.00 / 200.15 ms | 326.50 / 380.25 ms |
| Firefox parse | 0.020500 / 0.022575 ms/mol | 0.554200 / 0.655055 ms/mol |
| Firefox parse + canonical write | 0.601800 / 0.666435 ms/mol | 0.787100 / 0.944935 ms/mol |
| Firefox parse + radius-2 FP (9,999 supported rows) | 2.284678 / 2.501260 ms/mol | 0.753225 / 0.859701 ms/mol |
| WebKit init | 19.00 / 26.35 ms | 40.00 / 54.25 ms |
| WebKit parse | 0.002900 / 0.003430 ms/mol | 0.108500 / 0.121405 ms/mol |
| WebKit parse + canonical write | 0.071850 / 0.090630 ms/mol | 0.136150 / 0.164880 ms/mol |
| WebKit parse + radius-2 FP (9,999 supported rows) | 0.256576 / 0.314631 ms/mol | 0.123162 / 0.153860 ms/mol |

Their raw records are
`validation/results/competitive-browser-rdkitjs-published-npm-10k-firefox-20rep-v1.0.15-2026-09-16.json`
and
`validation/results/competitive-browser-rdkitjs-published-npm-10k-webkit-20rep-v1.0.15-2026-09-16.json`.
These published-package lanes neither establish end-to-end download time nor
peak/unique memory, and do not justify broad fingerprint superiority.

In the same three isolated Chrome runs, the browser-exposed `usedJSHeapSize`
after the workload had a median of 6,924,779 bytes for chematic and 25,429,694
bytes for RDKit.js (before-load medians: 2,132,108 and 2,352,144 bytes). This is
a snapshot, not peak memory. chematic's exported Wasm page count measured
1,507,328 bytes after initialization and 1,900,544 bytes after the workload.
RDKit.js MinimalLib does not expose its `WebAssembly.Memory`, so its linear
memory is `not_measured`; WebKit/Firefox do not expose the same JS heap metric,
and their process RSS remains `not_measured`. A separate Chrome CDP lane launched
each arm under an owned parent PID and sampled its complete browser process tree
every 50 ms. Across three 10k runs, peak process-tree RSS had a median of
2,178,007,040 bytes for chematic and 2,316,894,208 bytes for RDKit.js. The raw
record is `validation/results/competitive-browser-rdkitjs-isolated-10k-rss-v1.0.15-2026-09-16.json`.
This is summed `ps` RSS (shared pages may be counted by multiple processes), not
unique physical memory or a cross-browser result.

Firefox and WebKit then reran the same fixed corpus with **20** fresh-process
repetitions per arm. Their direct bit gates also produce **9,999/9,999 exact on the same
supported domain**, with the same one explicit Fe(II) refusal as Chromium. This
rules out a browser-specific fingerprint discrepancy for this configuration; it
does not claim support for RDKit's coordination-bond sanitization semantics.

| Browser / metric (p50 / p95 of per-process means) | chematic current source | RDKit.js 2026.03.6 |
|---|---:|---:|
| Firefox init | 167.00 / 202.05 ms | 338.00 / 381.85 ms |
| Firefox parse | 0.019800 / 0.021635 ms/mol | 0.580700 / 0.654465 ms/mol |
| Firefox parse + canonical write | 0.586950 / 0.666580 ms/mol | 0.813550 / 0.916620 ms/mol |
| Firefox parse + radius-2 FP | 1.392189 / 1.560551 ms/mol | 0.772077 / 0.859081 ms/mol |
| WebKit init | 21.50 / 26.15 ms | 39.50 / 59.35 ms |
| WebKit parse | 0.003150 / 0.003600 ms/mol | 0.116450 / 0.129480 ms/mol |
| WebKit parse + canonical write | 0.078700 / 0.091595 ms/mol | 0.157600 / 0.179365 ms/mol |
| WebKit parse + radius-2 FP | 0.171767 / 0.214716 ms/mol | 0.127663 / 0.168977 ms/mol |

The 20-repetition raw records are
`validation/results/competitive-browser-rdkitjs-isolated-10k-firefox-20rep-v1.0.15-2026-09-16.json`,
`validation/results/competitive-browser-rdkitjs-isolated-10k-webkit-20rep-v1.0.15-2026-09-16.json`,
`validation/results/competitive-browser-rdkitjs-ecfp4-parity-10k-firefox-v1.0.15-2026-09-16.json`,
and
`validation/results/competitive-browser-rdkitjs-ecfp4-parity-10k-webkit-v1.0.15-2026-09-16.json`.

## Interpretation boundaries

- This demonstrates a smaller WASM binary, lower measured cold init, and lower
  **local no-store** download-to-ready time in this browser/environment. It
  does **not** establish internet/CDN end-to-end download cost, since both arms
  were served from the same localhost server.
- RDKit.js is faster for this parse-inclusive fingerprint operation. Do not use
  an initialization or parse result to claim universal speed superiority.
- The operation configurations are aligned and the companion bit gate agrees in
  all three browsers on 9,999 declared-supported rows; the one recorded Fe(II)
  coordination input is a typed refusal. This does not establish
  equivalent stereochemistry, process RSS, WASM linear memory, another OS, or
  a general fingerprint guarantee.
- Chrome process-tree RSS is available only as a same-host, same-launch-method
  comparison. It is not unique resident memory, and no Firefox/WebKit RSS or
  cross-host memory conclusion follows from it.
- The timing lanes have 20 repetitions per arm but remain a single-host study.
  The fixed exposed 10,000-row corpus, three local browser engines, and the
  Chromium local download-to-ready lane close only the corresponding local
  parts of the follow-up; remote-network, unique-memory, and cross-host
  evidence remain required before any broad performance claim.

## Reproduction

```sh
python3 scripts/bench_browser_wasm_vs_rdkit_isolated.py \
  --rdkit-package /path/to/node_modules/@rdkit/rdkit \
  --browser "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --rows 10000 --warmup 20 --repetitions 3 \
  --output validation/results/competitive-browser-rdkitjs-isolated-v1.0.15-YYYY-MM-DD.json
```

Replace `--browser ...` with `--engine firefox` or `--engine webkit` for the
Playwright-managed browser lanes. Keep every output as a separate browser arm.
