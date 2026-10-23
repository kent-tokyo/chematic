# Benchmark runners

This directory contains measurement drivers that are still called by CI,
current benchmark documentation, or a maintained evidence check. Run them
from the repository root and write generated results under `benchmarks/`.

The retained runners cover the public Python/Rust and browser comparisons,
Open Babel file-I/O checks, chemistry dumps, conformer quality, startup,
SMILES parsing, and the current hot-path gate. Rejected rows and semantic
differences must remain visible before timing is ranked.

Release-specific experiments and superseded comparison programs are kept in
their Git tags, not in the current tree. A historical result does not by
itself keep its old runner active; use its recorded source revision when an
exact reproduction is required.

Historical benchmark records retain source hashes from the revisions that
produced them. For a new measurement, use the drivers in this directory.
