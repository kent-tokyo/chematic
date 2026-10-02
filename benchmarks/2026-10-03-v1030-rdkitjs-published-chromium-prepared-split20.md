# Published v1.0.30 browser prepared-fingerprint split

This follow-up separates first-use from reused-object Morgan calculations.
The earlier [Chromium record](2026-10-03-v1030-rdkitjs-published-chromium-paired20.md)
mixed those calls because its warm-up touched 20 of the 250 measured objects.
Its raw data and decision remain unchanged.

The pinned published `@kent-tokyo/chematic` 1.0.30 npm tarball and official
`@rdkit/rdkit` 2026.03.6, their JS/WASM hashes, and the first 250 rows of
`scripts/descriptor_census_corpus.smi` are checked by
`scripts/check_published_browser_paired.py`. Apple M4, macOS 27.0.1,
Chromium 140.0.7339.186, Python 3.13.6, and Playwright 1.58.0 were used.
Twenty fresh-process pairs alternate arm order. Each timed lane runs one
250-row batch; the paired log-ratio 95% intervals use 10,000 deterministic
bootstrap resamples. A ratio above 1 favors CheMatic.

For each arm, 20 *separate* prepared warm-up objects are fingerprinted before
timing all 250 measured prepared objects for the first time. The same 250
objects are then fingerprinted again for the reused-object lane. Parsing and
preparation are outside both prepared timing windows. Both arms materialize
the same 256-byte packed output; RDKit.js bit-string packing remains inside
its timed API path. This is an end-to-end API comparison, not a claim that the
internal fingerprint cores do equal work.

| Lane | Median speed ratio | Paired 95% interval | Output gate |
|---|---:|---:|---|
| Parse + compatible Morgan | 7.63× | 7.47–7.79 | Aggregate packed fingerprint digest exact |
| Prepared Morgan, first use | 6.52× | 6.46–6.76 | Aggregate packed fingerprint digest exact |
| Prepared Morgan, reused object | 6.83× | 6.57–6.95 | Aggregate packed fingerprint digest exact |
| Parse only | 19.20× | 19.14–19.65 | Diagnostic; browser output not checked row by row |
| Parse + canonical write | 3.32× | 3.29–3.37 | Diagnostic; canonical strings not equivalent |

The packed Morgan digest is `5c006dbe` with 11,307 set bits in every
fingerprint lane of every browser run. A separate published Node preflight
checks these artifact hashes at row level (250/250 atom counts and Morgan bit
vectors); this browser run itself checks only the aggregate digest. It does
not prove full browser row-level equivalence.

Whole Chromium process-tree RSS was sampled every 50 ms. Median peak was
1,235,181,568 bytes for CheMatic and 1,298,235,392 bytes for RDKit.js.
This may double-count shared pages and does not measure library allocations;
it is not incorporated into speed ratios. Init medians were 17.4 and 72.6 ms,
respectively, in this one-host local-file test.

The [raw 20-pair record](2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.json)
has SHA-256 `feaa0a80067e18e34e82c3048ef0aed626bfbe1ae4f843b2a95f6235403e2af7`.
The [checked summary](2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20-summary.json)
has SHA-256 `16cc70762572b0b953418826fadb90172529bb0fb2dee72a138df3d51cd33bd8`.
The artifact-packet checker pins both records. Reproduce with the command in
the earlier Chromium record plus `--prepared-mode split` and a new output
path, then run `scripts/check_published_browser_paired.py` on that output.

P0.2 stays open: an independent host/browser, browser row-level output gate,
equal-work perception boundaries across the remaining operations, and
matched library-memory accounting are still missing. These values are not
general browser speed or per-million-molecule memory claims.
