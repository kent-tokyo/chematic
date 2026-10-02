# Published v1.0.30 vs RDKit.js in Chromium — 20 paired blocks

This is a single-host browser follow-up to the [published npm/Node packet](2026-10-03-v1029-v1030-rdkitjs-node-isolated-paired20.md).
It uses the published `@kent-tokyo/chematic` 1.0.30 npm tarball (SHA-256
`fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201`)
and official `@rdkit/rdkit` 2026.03.6. The checked record pins both packages'
JS/WASM file hashes. The input is the first 250 rows of
`scripts/descriptor_census_corpus.smi` (whole-file SHA-256
`d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`).

Environment: Apple M4, macOS 27.0.1, Chromium 140.0.7339.186, Python
3.13.6, Playwright 1.58.0. Each observation loads one library in a fresh
Chromium process through no-store localhost routes. Twenty blocks alternate
CheMatic→RDKit.js and RDKit.js→CheMatic; each has 20 warm-up rows. To avoid
the browser timer's per-molecule granularity, each timed operation runs as a
single 250-row batch. A deterministic 10,000-resample bootstrap of paired
log-ratios gives the 95% intervals. A ratio above 1 favors CheMatic.

| Operation | CheMatic / RDKit.js median | Paired 95% CI | Claim boundary |
|---|---:|---:|---|
| Parse only | 19.48× | 19.15–19.73 | Diagnostic; browser parse output was not checked row by row |
| Parse + compatible Morgan | 7.60× | 7.53–7.65 | Fingerprint digest exact in every run |
| Prepared Morgan | 6.58× | 6.36–6.71 | Fingerprint digest exact in every run; first-use and reused molecules mixed |
| Parse + canonical write | 3.33× | 3.28–3.32 | Diagnostic only; canonical strings were not checked for equality |

The browser produces the same packed radius-2/2,048-bit fingerprint digest
(`5c006dbe`, 11,307 set bits) on both arms in every run. Separately, the
published Node preflight verifies 250/250 atom counts and Morgan bit vectors
row by row on these same artifact hashes. The browser record itself checks an
aggregate digest, not browser row-by-row equality. The prepared lane warms the
first 20 of 250 already-prepared molecules, so it is a mixed first/reuse lane,
not an isolated first-use or reuse claim. RDKit.js bit strings are packed to
the common byte representation inside its timed JS API path; this is an
end-to-end API comparison, not a claim that the underlying fingerprint cores
do equal internal work.

Whole Chromium process-tree RSS is sampled every 50 ms and retained per run.
Its median peak is 1,234,886,656 bytes for CheMatic and 1,294,548,992 bytes
for RDKit.js. Shared pages may be counted more than once; the numbers are not
library allocations or per-molecule memory. RDKit.js also does not expose
its Wasm linear-memory object to this harness. Startup and RSS are not folded
into the speed ratios.

The [raw record](2026-10-03-v1030-rdkitjs-published-chromium-paired20.json)
has SHA-256 `de0e7f1f8ba7ad208aeb906acb435b10ca210c05bf4ddc860a4217618c99275d`.
The [checked summary](2026-10-03-v1030-rdkitjs-published-chromium-paired20-summary.json)
has SHA-256 `3f27f1c17824aabfee7b52ae9c122974736290eb00da7249c29f5053823d2233`.
`scripts/check_published_browser_paired.py` recomputes the intervals and
validates all 40 raw runs, fingerprint digests, package hashes and memory
observations. The v1.0.30 artifact-packet checker pins both records.

```sh
python3 scripts/bench_browser_wasm_vs_rdkit_isolated.py \
  --rdkit-package OUT/rdkit/package \
  --schematic-dir OUT/v30/package \
  --schematic-tarball OUT/v30/kent-tokyo-chematic-1.0.30.tgz \
  --corpus scripts/descriptor_census_corpus.smi --rows 250 --warmup 20 \
  --repetitions 20 --batch-timing --engine chromium \
  --browser /path/to/chromium --measure-process-rss --output OUT/browser.json
python3 scripts/check_published_browser_paired.py \
  --input OUT/browser.json --output OUT/browser-summary.json
```

This does not transfer to other browsers, hosts, molecule sets, canonical
writing, all 63 operations, or matched library-memory accounting. The
[separate 20-pair follow-up](2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.md)
now isolates prepared first-use and reused-object calls without changing this
historical record. P0.2 stays open for the other checks.
