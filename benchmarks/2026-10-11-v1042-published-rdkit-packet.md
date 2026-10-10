# Published v1.0.42 versus RDKit 2026.09.1

This packet reruns the four named Python lanes from the v1.0.42 release-period
source-candidate record using exact PyPI wheels. It also checks the equivalent
parts of the published npm surface. The result is scoped to the listed
artifacts, corpus, operations, runners, and runtimes. It is not a universal
RDKit-parity or speed claim.

## Method

- Corpus: `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`
  (`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`)
- Accuracy: all 10,000 rows before timing eligibility is decided
- Timing: first 1,000 rows, 21 alternating-order blocks
- Lanes: prepared first call and parse-inclusive pipeline, kept separate
- Interval: paired bootstrap of the per-block RDKit/CheMatic ratio, 10,000
  resamples; values above 1 favor CheMatic
- Hosts: GitHub-hosted Ubuntu 24.04 x64 and macOS 14 arm64

Each raw record includes the package filename, byte count, SHA-256, runtime,
platform, raw block samples, mismatch rows, and failure counts. A speed
interval is omitted when the output gate fails or the exposed APIs perform
different work.

## Python output gates

Both CPython 3.13 runners produced the same accuracy result against RDKit
2026.9.1:

| Operation | Result | Speed eligible |
|---|---:|---:|
| TPSA | 10,000 / 10,000 exact | yes |
| Labute ASA | 10,000 / 10,000 within `1e-9` | yes |
| Chiral Morgan R2 / 2,048 bits | 10,000 / 10,000 exact | yes |
| Ring count | 9,994 / 10,000 exact | **no** |

The six ring-count residuals are the large symmetric macrocycles fixed in the
later source candidate. Published v1.0.42 must therefore not inherit the
candidate's 10,000/10,000 ring result or any ring speed statement.

## Python timing

### Ubuntu 24.04 x64

| Operation | Prepared ratio (95% CI) | Parse-inclusive ratio (95% CI) |
|---|---:|---:|
| TPSA | 4.270x (3.933–4.518), 21/21 | 18.056x (17.612–18.243), 21/21 |
| Labute ASA | 5.717x (4.566–5.956), 16/21 | 16.834x (16.561–17.094), 21/21 |
| Chiral Morgan R2 | 1.901x (1.804–1.951), 21/21 | 5.341x (5.186–5.424), 21/21 |
| Ring count | withheld: output mismatch | withheld: output mismatch |

### macOS 14 arm64

| Operation | Prepared ratio (95% CI) | Parse-inclusive ratio (95% CI) |
|---|---:|---:|
| TPSA | 1.385x (1.129–1.499), 16/21 | 16.365x (10.250–18.211), 21/21 |
| Labute ASA | 0.547x (0.506–0.759), 2/21 | 11.627x (10.782–15.142), 21/21 |
| Chiral Morgan R2 | 1.038x (0.878–1.205), 11/21 | 4.077x (3.875–4.475), 21/21 |
| Ring count | withheld: output mismatch | withheld: output mismatch |

The parse-inclusive TPSA, Labute ASA, and chiral Morgan pipelines have a
favorable interval on both runners. The prepared-operation result is not a
cross-platform win: macOS Labute ASA loses, and macOS chiral Morgan crosses
parity. Host-specific results must remain host-specific.

## Node/WASM boundary

Published npm v1.0.42 has 10,000/10,000 Labute ASA values within `1e-9`, but
RDKit.js exposes that value only through `get_descriptors()`, which computes a
descriptor bundle. Timing it against CheMatic's single Labute getter is not
equivalent work, so no ratio is reported.

The published chiral Morgan lane agrees on 9,999/10,000 rows. Row 2960 differs,
so both prepared and pipeline speed intervals are withheld. npm v1.0.42 does
not expose the Python packet's explicit `rdkit_tpsa` or `rdkit_num_rings`
profiles; its native TPSA and `ring_count` surfaces are not substituted.

No Rust/RDKit timing ratio is added. v1.0.42 crates are source archives and
there is no published same-runtime RDKit Rust binding for this packet;
comparing a Rust binary with Python would change the runtime boundary.

## Difference from the source candidate

The source candidate's exact ring-count result does not apply to the published
v1.0.42 wheel. Its direct every-block TPSA and chiral Morgan result also does
not reproduce as a two-platform published-artifact statement: the macOS
chiral-Morgan interval crosses parity and wins only 11/21 blocks. The strictly
supported public statement is limited to the three Python parse-inclusive
pipelines on the two recorded hosted runners.

## Raw records

- [Python, Linux x64](2026-10-11-v1042-published-rdkit-python-linux-x64.json)
- [Python, macOS arm64](2026-10-11-v1042-published-rdkit-python-macos-arm64.json)
- [Node/WASM, Linux x64](2026-10-11-v1042-published-rdkit-node-linux-x64.json)
- [Node/WASM, macOS arm64](2026-10-11-v1042-published-rdkit-node-macos-arm64.json)
- [GitHub Actions run](https://github.com/kent-tokyo/chematic/actions/runs/38086615954)
