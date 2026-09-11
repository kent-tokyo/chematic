# v1.0.12 streaming format safety gate

Current local safety evidence for the ten common-runner formats on macOS arm64.

Command:

```text
TMPDIR=/private/tmp python3 scripts/check_streaming_format_limits.py \
  --binary target/debug/examples/streaming_benchmark \
  --write-evidence --evidence-date 2026-09-11
```

| Gate | Cases | Failures |
|---|---:|---:|
| malformed negative-input attempts | 800 | 0 |
| unique malformed payloads | 800 | n/a |
| duplicate payload reuses | 0 | n/a |
| oversized inputs | 10 | 0 |
| gzip controls and post-decompression limits | 20 | 0 |
| generated parser-entry cases | 480 | 0 |

The generated parser-entry wave covers eight format-specific families with six
byte-distinct cases per family. All 80 format/category combinations produced a
typed failure kind. This is bounded parser-safety evidence; exhaustive
malformed parser-state coverage and cross-engine equivalence remain open.

The separate failure taxonomy also probes twelve Extended XYZ, mmCIF, and
MOL2 branches in addition to the 120-case base corpus, for 132 typed-error
observations.

Machine-readable evidence:
[`2026-09-11-streaming-safety-v1.0.12.json`](2026-09-11-streaming-safety-v1.0.12.json),
[`streaming-parser-entry-categories-v1.0.12.json`](../validation/results/streaming-parser-entry-categories-v1.0.12.json),
and
[`streaming-parser-entry-failure-kinds-v1.0.12.json`](../validation/results/streaming-parser-entry-failure-kinds-v1.0.12.json).
