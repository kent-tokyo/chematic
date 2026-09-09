# File-backed SDF comparison — 2026-09-07

Status: local, version-pinned evidence for the P1 comparison gate. This is a
same-input file-backed ingestion comparison, not a claim of chemistry-semantic
parity or a general-purpose speed ranking.

Environment: macOS arm64, workspace `1.0.8`, release-mode chematic binary,
checked-in `benchmarks/fixtures/streaming.sdf` (two records, 633 bytes), 20
repetitions. RDKit was `2025.09.3`; Open Babel was `3.2.1`.

## Results

| Engine | Operation boundary | Records | Failures | Seconds | Records/s |
| --- | --- | ---: | ---: | ---: | ---: |
| chematic | Rust `SdfFileReader`, one release process | 40 | 0 | 0.000439 | 91,150.85 |
| RDKit | `ForwardSDMolSupplier`, one Python process | 40 | 0 | 0.016774 | 2,384.68 |
| Open Babel | CLI process per repetition, SDF→SMILES | 40 | 0 | 6.251115 | 6.40 |

The chematic and RDKit rows both read the same file through file-backed APIs,
but their parser policies are not identical: chematic retains its diagnostic
path and RDKit was run with `sanitize=False` and `removeHs=False`. Therefore
the rows establish comparable ingestion evidence and record counts, not a
causal chemistry or throughput superiority claim. Open Babel includes process
startup and CLI conversion/output handling and is explicitly not a same-process
comparison.

## Reproduction

```text
cargo run --release -p chematic-mol --example streaming_benchmark --offline -- \
  --format sdf --path benchmarks/fixtures/streaming.sdf --repeats 20
python3 scripts/bench_streaming_formats.py --mode file-backed \
  --sdf benchmarks/fixtures/streaming.sdf --repeats 20 \
  --openbabel /opt/homebrew/bin/obabel
```

The existing block-constructor lane remains in
[`2026-09-04-streaming-formats.md`](2026-09-04-streaming-formats.md) and is
not combined with these file-backed rows.
