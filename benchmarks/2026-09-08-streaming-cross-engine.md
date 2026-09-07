# Same-input SDF cross-engine contract — 2026-09-08

Status: local v1.0.9 development-candidate evidence for P1. This report checks
record/failure agreement on one identical file-backed SDF fixture. It is not a
chemistry-semantic parity result or a general throughput ranking because the
three engines have different parser and process boundaries.

Environment: macOS arm64, checked-in `benchmarks/fixtures/streaming.sdf`
(633 bytes, SHA-256 `627e7814694f7d00a9212712ec697f2174d30327b459d983d2a10f717cae0465`),
20 repetitions; RDKit `2025.09.3`; Open Babel `3.2.1`.

| Engine | Boundary | Records | Failures | Input bytes | Records/s |
| --- | --- | ---: | ---: | ---: | ---: |
| chematic | `SdfFileReader` over file-backed `BufRead` | 40 | 0 | 12,660 | 34,355.64 |
| RDKit | `ForwardSDMolSupplier` in one Python process | 40 | 0 | 12,660 | 111.85 |
| Open Babel | CLI conversion per repetition, including process startup | 40 | 0 | 12,660 | 13.46 |

The contract passed: all engines returned the expected 40 records and zero
failures for the same input bytes. The executable check is:

```text
python3 scripts/check_streaming_cross_engine.py \
  --repeats 20 \
  --openbabel /opt/homebrew/bin/obabel \
  --output /tmp/chematic-streaming-cross-engine-20.json
```

The JSON report preserves the effective limits and each engine's explicit
comparison boundary. Process startup, conversion, parser policy, and Python
object construction are not normalized into a speed claim.
