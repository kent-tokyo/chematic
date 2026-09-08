# 2026-09-08 same-input SDF streaming comparison

This is a record-accounting and boundary comparison for the checked-in
`benchmarks/fixtures/streaming.sdf` fixture. It is not a same-condition speed
claim: chematic uses its Rust file-backed `BufRead` reader, RDKit uses one
Python `ForwardSDMolSupplier`, and Open Babel is invoked as a fresh CLI
process for each repetition.

## Fixture and protocol

- Fixture: `benchmarks/fixtures/streaming.sdf`
- Fixture SHA-256: `627e7814694f7d00a9212712ec697f2174d30327b459d983d2a10f717cae0465`
- Source bytes per repetition: 633
- Repetitions: 20
- Expected records: 40
- Command: `python3 scripts/check_streaming_cross_engine.py --format sdf --repeats 20 --openbabel obabel`
- Open Babel: `Open Babel 3.2.1 -- Jul 11 2026`
- RDKit: `2025.09.3`

## Results

| Engine | Boundary | Records | Failures | Seconds | Records/s |
| --- | --- | ---: | ---: | ---: | ---: |
| chematic | Rust file-backed `BufRead` | 40 | 0 | 0.003720 | 10,753.17 |
| RDKit | Python `ForwardSDMolSupplier` | 40 | 0 | 0.353902 | 113.03 |
| Open Babel | CLI conversion per repetition, startup included | 40 | 0 | 2.929874 | 13.65 |

All three lanes agree on record accounting for this fixture. The timings must
not be used as parser-only or same-process throughput comparisons because the
execution boundaries differ.
