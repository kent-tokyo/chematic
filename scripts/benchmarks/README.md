# Benchmark runners

This directory contains retained measurement drivers. Run them from the
repository root. Their dated records and required inputs are documented in
`benchmarks/`; generated measurements belong there, not beside these scripts.

The Open Babel drivers cover file I/O, SDF and V3000 throughput, and MOL2
indexing. The throughput runners share the same semantic checks and keep
rejected rows visible before reporting timing. The original canonical,
SMARTS, SDF, descriptor, and streaming runners are also here, with their
existing command options. `bench_hotpath_pair.py` has a unit test in
`scripts/tests/`.

The `bench5k`, startup, SMILES parsing, browser/WASM, published Python/Rust,
and streaming-matrix drivers also run from here. Current workflow and report
commands use these paths. `scripts/bench_published_wasm_paired.mjs` remains at
its original path because the checked-in published-WASM evidence gate compares
its exact source hash. The source-pinned Python Morgan runner likewise stays
at the root; moving either would invalidate the retained packet check.

Historical benchmark records retain source hashes from the revisions that
produced them. For a new measurement, use the drivers in this directory.
