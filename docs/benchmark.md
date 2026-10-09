# Benchmarks

CheMatic benchmarks are specific to the package or source commit, comparator,
corpus, host and timed operation. A speed result counts only when the stated
outputs and work are equivalent. [The record index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md)
links the raw data and older measurements.

## Current published-package evidence

The current release is **v1.0.40**. The latest broad published-package performance
packet below is for **v1.0.30**, not a new v1.0.40 measurement. The v1.0.40
RDKit/COSMolKit and seeded ETKDG record is source-built evidence. None is a
universal speed or compatibility claim.

| Lane | Result | Boundary |
|---|---|---|
| [Python 63-operation matrix](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) | 20 operations have exact output agreement and a favorable paired confidence interval against RDKit 2026.03.6 in 20 alternating blocks. | One host and exposed corpus; other operations are not wins. |
| [Isolated Python](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1.0.30-isolated-python-time-memory.md) and [Rust](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1029-v1030-published-rust-isolated-time-memory.md) | Parse-inclusive, prepared first-use and precomputed lanes use 20 fresh-process blocks. | Process peak RSS is not per-operation allocation. Changed HBA outputs cannot support an equivalent-output v1.0.29/v1.0.30 speed claim. |
| [Published Node/WASM](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1029-v1030-rdkitjs-node-isolated-paired20.md) | Compatible Morgan is bit-identical on the 250-row speed corpus, with 20 fresh-process paired blocks against official RDKit.js. | Node is not browser evidence. |
| [Published browser replication](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md) | Chromium, Firefox and WebKit on Ubuntu 24.04 each have 250/250 exact direct/prepared Morgan rows and favorable 20-block speed intervals; the smallest lower bound is 1.63×. | Browser, host and operation scoped. A separate Chromium 10k lane has 9,999 exact outputs and one typed Fe refusal. |

The [v1.0.30 artifact audit](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md)
accounts for 10,000 chemistry rows, 310,000 SMARTS cells, 57 legacy reactions
and applicable operation outputs across Python, npm and Rust. It is an output
audit, **not** a speed or full RDKit-parity result. The
[adoption policy](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md)
retains 200 SMARTS failures, five typed CIP abstentions and four npm API gaps.

## Historical records

- The [v1.0.20 public-package fingerprint/3D record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-23-public-package-fingerprint-3d-v1.0.20.md)
  reports 9,999 supported exact Morgan rows and one typed refusal. Its
  quality-equivalent MMFF94 speed ratio was 0.944× (95% lower bound 0.861×):
  **not** an MMFF94 speed win.
- The [v1.0.15 WASM-size record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md)
  measured 4,005,280 raw / 1,460,499 gzip bytes versus official RDKit.js
  7,333,095 / 2,379,975 under one local compression method. These are
  historical assets, not current package or download sizes.
- The [v1.0.12 similarity-search record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-11-similarity-search-v1.0.12.md)
  separates native/native, RDKit-compatible/RDKit and cross-profile overlap.
  Source A/B and older operation timings remain in the
  [index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md); they do not update v1.0.40 claims.

## Reproduction and interpretation

Run the checker named in the selected record with **its** pinned packages,
corpus hash, configuration and host. For repository consistency checks:

```bash
python3 scripts/check_benchmark_index.py
python3 scripts/check_v1030_artifact_packet.py
python3 scripts/check_published_python_isolated_paired.py
python3 scripts/check_published_rust_isolated_paired.py
```

A new measurement must record artifact/source hashes; corpus and failure
accounting; hardware, OS and runtime; timed boundary and warm-up; raw paired
observations, aggregation and uncertainty; and an output-equivalence gate.
Keep source candidates, published packages, sealed accuracy data, whole-process
RSS and library allocations distinct. If outputs differ, the interval crosses
parity or work is asymmetric, report **no equivalent-output speed win**.
