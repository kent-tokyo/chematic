# Open issue triage — 2026-09-08

This snapshot records the bounded completion of issue #507. It does not turn
format-specific safety cases into a throughput, compatibility, or universal
parser claim.

## Issue #507 — completed

The checked-in `validation/streaming_format_safety_cases.json` corpus now has
ten malformed-input cases for each supported runner format: SDF, V2000 MOL,
XYZ, V3000, MOL2, CML, CDXML, mmCIF, and PDB. The dependency-free gate checks
all 90 cases and requires exactly one failed record for each case.

The same gate retains nine oversized-input rejections and 18 gzip controls:
one valid decompressed control and one post-decompression input-limit rejection
for every format. CML, CDXML, and PDB remain explicitly line-limit safety
checks because their current readers are intentionally lenient about unknown,
empty, or non-record input.

## Evidence

- `python3 scripts/check_streaming_format_limits.py` — 90 negative cases,
  9 oversized cases, and 18 gzip cases passed.
- The corpus schema and exact ten-case-per-format count are validated before
  any runner invocation.

This closes the bounded malformed-corpus expansion in #507. Cross-language
streaming parity, equivalent cross-engine throughput, and broader parser
semantics remain separate open roadmap gates.
