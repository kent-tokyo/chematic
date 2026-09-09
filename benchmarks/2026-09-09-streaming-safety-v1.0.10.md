# v1.0.10 streaming format safety gate

Local safety evidence for the ten supported common-runner formats on macOS
arm64.

Command:

```text
python3 scripts/check_streaming_format_limits.py
```

The twelve checked-in base cases for each format are classified in
[`validation/streaming_format_safety_categories.json`](../validation/streaming_format_safety_categories.json).
The gate requires every format's declared categories to be represented before
the supplemental parser-path cases are executed; this guards category
coverage in addition to the aggregate case-count minimum.

| Gate | Cases | Failures |
|---|---:|---:|
| malformed negative-input attempts | 800 | 0 |
| unique malformed payloads | 800 | n/a |
| duplicate payload reuses | 0 | n/a |
| oversized inputs | 10 | 0 |
| gzip controls and post-decompression limits | 20 | 0 |

Malformed attempts by format: CDXML 74, CML 74, Extended XYZ 74, mmCIF 74,
MOL 74, MOL2 74, PDB 74, SDF 74, V3000 74, and XYZ 74.
The gate also enforces a minimum of 50 malformed attempts per format. The
supplemental waves are byte-distinct (repeated boundaries receive a carriage
return variant without adding another logical record), so all 800 attempts are
unique payloads. Expanding parser-state coverage remains part of the separate
exhaustive-corpus work.

The generated parser-entry wave adds 480 byte-distinct cases (48 per format)
and independently reports zero duplicate reuses. This is an additional
deterministic check using eight format-specific malformed shapes per format,
including a format-shaped truncated-record family; it is not a claim of
exhaustive parser-state coverage. Runner process failures now retain the
captured stderr/stdout in the diagnostic exception.

The corpus includes SDF/MOL, XYZ/Extended XYZ, V3000, MOL2, CML, CDXML,
mmCIF, and PDB. This is a bounded local safety gate; exhaustive malformed
coverage and cross-engine equivalence remain separate roadmap items.

Machine-readable evidence: [`2026-09-09-streaming-safety-v1.0.10.json`](2026-09-09-streaming-safety-v1.0.10.json).
