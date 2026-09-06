# SDF/MOL/XYZ streaming benchmark — 2026-09-04

Status: local reproducibility evidence. This is not yet a fair cross-engine
throughput claim because the chematic runner uses file-backed `BufRead`
streaming while the available RDKit comparator uses Python block constructors.

Environment: macOS arm64, chematic workspace `1.0.3`, Rust release profile,
200 repetitions over the checked-in two-record SDF and two-frame XYZ fixtures.
The RDKit comparator was Python package `2025.09.3`.

## chematic file-backed streaming

| Format | Records | Failures | Records/s | Bytes/s |
| --- | ---: | ---: | ---: | ---: |
| SDF | 400 | 0 | 94,580 | 29,934,683 |
| MOL (SDF blocks) | 400 | 0 | 39,387 | 12,465,998 |
| XYZ | 400 | 0 | 118,231 | 9,517,591 |

Commands:

```text
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format sdf --path benchmarks/fixtures/streaming.sdf --repeats 200
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format mol --path benchmarks/fixtures/streaming.sdf --repeats 200
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format xyz --path benchmarks/fixtures/streaming.xyz --repeats 200
```

## RDKit block-parser reference

`python3 scripts/bench_streaming_formats.py --repeats 200` produced:

| Format | Records | Failures | Records/s | Bytes/s |
| --- | ---: | ---: | ---: | ---: |
| SDF | 400 | 0 | 440 | 133,641 |
| MOL | 400 | 0 | 5,775 | 1,752,616 |
| XYZ | 400 | 0 | 443,255 | 35,681,991 |

These rows are a parser reference only: they do not measure RDKit's file
supplier or a common streaming interface. A fair comparison requires a
larger identical corpus, equivalent file-backed APIs, repeated independent
runs, RSS/allocation capture, and malformed/oversized cases.

The runner now accepts explicit `--max-input-bytes`, `--max-record-bytes`,
`--max-line-bytes`, `--max-records`, and `--max-atoms` limits and records the
effective limits in its JSON output. This makes malformed and oversized SDF,
MOL, and XYZ boundaries reproducible without changing parser code. The runner
also accepts V3000, MOL2, CML, CDXML, and mmCIF as explicitly reported
`execution_mode: "materialized_one_shot"` parser rows. They are not
file-backed streaming or cross-engine throughput evidence; PDB and gzip remain
separate format-adapter work.

Example bounded runs:

```text
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- \
  --format sdf --path benchmarks/fixtures/streaming.sdf --repeats 1 \
  --max-input-bytes 1024 --max-record-bytes 4096 --max-line-bytes 1024 --max-records 2
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- \
  --format xyz --path benchmarks/fixtures/streaming.xyz --repeats 1 \
  --max-input-bytes 1024 --max-line-bytes 1024 --max-records 2 --max-atoms 100
```

Example format-parser coverage:

```text
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format v3000 --path benchmarks/fixtures/ethanol.v3000 --repeats 1
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format mol2 --path benchmarks/fixtures/ethanol.mol2 --repeats 1
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format cml --path benchmarks/fixtures/ethanol.cml --repeats 1
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format cdxml --path benchmarks/fixtures/ethanol.cdxml --repeats 1
cargo run -p chematic-mol --example streaming_benchmark --release --offline -- --format mmcif --path benchmarks/fixtures/minimal.mmcif --repeats 1
```
