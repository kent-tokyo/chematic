# Benchmark runners

This directory contains retained measurement drivers. Run them from the
repository root. Their dated records and required inputs are documented in
`benchmarks/`; generated measurements belong there, not beside these scripts.

The Open Babel drivers cover file I/O, SDF and V3000 throughput, and MOL2
indexing. The throughput runners share the same semantic checks and keep
rejected rows visible before reporting timing. Historical benchmark records retain the source paths and hashes from the
revisions that produced them. In those records, the former
`scripts/bench_openbabel_*.py` paths refer to the historical checkout; for a
new measurement, use the drivers in this directory.
