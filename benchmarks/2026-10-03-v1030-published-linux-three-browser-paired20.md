# Published v1.0.30 browser comparison on Linux

GitHub Actions [run 37066587611](https://github.com/kent-tokyo/chematic/actions/runs/37066587611)
on PR #710 measured the published `@kent-tokyo/chematic` 1.0.30 npm tarball
(SHA-256 `fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201`)
against official `@rdkit/rdkit` 2026.03.6. The job verifies both packages'
JS/WASM hashes and uses Ubuntu 24.04, Python 3.12, Node 22, and Playwright
1.58.0. Each browser has its own ephemeral Linux runner; runner CPU allocation
and exact browser binary versions were not frozen in the raw record.

The first 250 rows of `scripts/descriptor_census_corpus.smi` were measured in
20 alternating fresh-process pairs per browser, with 20 separate warm-up
objects and single-batch timing. A fixed-seed, 10,000-resample bootstrap of
paired log-ratios gives the 95% intervals. Ratios above 1 favor CheMatic.
Each browser separately checked every Morgan bit, molecule by molecule, for
both direct and prepared APIs: **250/250 exact for each operation in each
browser**. This is a correctness gate outside the timed path. The timed path
also checked its aggregate packed fingerprint digest on every run.

| Browser | Parse + Morgan | Prepared first use | Prepared reused object |
|---|---:|---:|---:|
| Chromium | 4.30× (3.65–4.37) | 2.49× (2.26–3.03) | 2.99× (2.74–3.82) |
| Firefox | 9.10× (9.00–9.14) | 4.78× (4.65–4.96) | 4.78× (4.61–4.83) |
| WebKit | 5.79× (4.66–6.08) | 1.65× (1.63–2.03) | 4.00× (3.82–4.51) |

The published [M4 Chromium record](2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.md)
uses different host hardware and browser build; its ratios must not be
averaged with these Linux results. These are three engine-specific results
from one CI run, not a population estimate over machines or molecules.
Parse-only and canonical-write timings remain diagnostics because the browser
runner does not gate those outputs for literal equality. RDKit.js bit-string
packing is included in its timed API path; the internal native fingerprint
cores are not proven to perform equal work.

Process-tree RSS was **not measured** in these Linux jobs. The M4 Chromium
process RSS record remains separate and does not attribute allocations to a
library. The comparison therefore cannot support a library-memory or
million-molecule-memory advantage claim.

The [workflow](../.github/workflows/published-browser-replication.yml)
reinstalls the published packages and browser runtime; it does not build the
current source tree as the benchmark candidate. The twelve checked-in raw
timing, summary, and per-row correctness files are linked from the
[benchmark index](README.md). Their SHA-256 values, package/corpus identity,
20-pair protocol, output digests, and 250/250 per-browser row checks are
verified by `python3 scripts/check_published_linux_browser_replication.py`
and included in `scripts/check_v1030_artifact_packet.py`.

P0.2 remains open for equivalent-work perception boundaries and the other
operations/bindings, representative additional corpora, and matched
library-memory accounting. This record proves a bounded Morgan result, not
that CheMatic is faster than RDKit across its complete API.
