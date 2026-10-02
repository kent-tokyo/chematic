# Published Rust v1.0.29 → v1.0.30 isolated timing and process RSS

This exposed 5,000-molecule ChEMBL packet compares **published crates.io
graphs**, not source builds. It answers whether v1.0.30 changed the runtime
cost of named Rust calls relative to v1.0.29. It does not compare Rust-call
time with RDKit's Python-call time, prove broad RDKit speed superiority, or
close P0.2.

## Artifact and protocol

- Apple M4, macOS 27.0.1 arm64; Rust/Cargo 1.97.0. The exact manifest and
  `Cargo.lock` under `tools/published_rust_gate*` pin every `chematic-*`
  crate to its respective published version. The runner checks lockfile and
  archived output hashes, builds each arm `--release --locked --offline`, and
  records the two executable SHA-256 values in every raw record.
- Exposed 5,000-row corpus SHA-256:
  `1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`.
  Each arm's 5,000 outputs must match its published Rust 63-operation archive
  on **every block**. HBA differs between versions on 1,359 inputs by design;
  compatible Morgan outputs are identical.
- Twenty alternating AB/BA blocks per lane, one discarded process per arm
  before measurement, and a separate warmup molecule inside each fresh
  process. `operation_ns` uses Rust `Instant` around 5,000 calls. Parse is
  included only in `parse_inclusive`; `prepared_first_use` parses beforehand,
  while `precomputed` computes once beforehand and then repeats calls on the
  same molecules. Preparation time is retained separately, not silently
  added to hot-call time.
- `/usr/bin/time -l` records each process's peak resident set size. This
  includes startup, parsing, preparation, timed calls and outputs; it is
  **not** net operation allocation. A deterministic 4,000-resample paired
  bootstrap gives 95% intervals for both time and whole-process RSS ratios.

`A/B` is v1.0.29 time divided by v1.0.30 time: values above 1 favor v1.0.30.

| Operation | Mode | Median A/B time | 95% interval | Interpretation |
|---|---|---:|---:|---|
| Compatible Morgan | Parse-inclusive | 1.006 | 0.993–1.011 | Output exact; interval crosses parity |
| Compatible Morgan | Prepared first-use | 1.002 | 0.991–1.009 | Output exact; interval crosses parity |
| Compatible Morgan | Precomputed | 0.989 | 0.973–1.009 | Output exact; interval crosses parity |
| HBA | Parse-inclusive | 1.014 | 0.997–1.028 | Output changed; no equivalent-output speed claim |
| HBA | Prepared first-use | 1.011 | 0.990–1.022 | Output changed; no equivalent-output speed claim |
| HBA | Precomputed | 0.993 | 0.944–1.006 | Output changed; no equivalent-output speed claim |

Whole-process RSS medians (v1.0.29/v1.0.30) range from 0.997 to 1.018 across
the six lanes. Raw peak bytes and paired intervals are retained in the six
JSON records; neither a per-operation memory reduction nor cross-binding
memory win follows from those ratios. This is one host and one exposed corpus.

## Reproduction and limitations

Run `scripts/bench_published_rust_isolated_paired.py` for each operation and
mode, then `python3 scripts/check_published_rust_isolated_paired.py`. The
checker revalidates pinned crate graphs and output archives, all 20 block
orders and 5,000-output digests, and recomputes both intervals and outcome
labels. [Raw HBA parse](2026-10-03-rust-v1029-v1030-hba-parse-inclusive-paired20.json),
[first-use](2026-10-03-rust-v1029-v1030-hba-prepared-first-use-paired20.json),
[precomputed](2026-10-03-rust-v1029-v1030-hba-precomputed-paired20.json),
[raw Morgan parse](2026-10-03-rust-v1029-v1030-morgan-parse-inclusive-paired20.json),
[first-use](2026-10-03-rust-v1029-v1030-morgan-prepared-first-use-paired20.json),
and [precomputed](2026-10-03-rust-v1029-v1030-morgan-precomputed-paired20.json)
retain every observation. The archived binary hashes identify the executables
used but do not substitute for an independent reproducible-build attestation.

Remaining P0.2 work includes other operations, npm/WASM, a matched
perception/precomputation boundary against RDKit, per-operation memory where
measurable, and an independent host. The older “21 operations faster” count
remains a historical diagnostic.
