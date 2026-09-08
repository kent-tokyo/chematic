# Open issue triage — 2026-09-08

This snapshot records the bounded completion of issue #507. It does not turn
format-specific safety cases into a throughput, compatibility, or universal
parser claim.

## Issue #303 — bounded structural slice completed

The explainable reactivity API now retains the existing epoxide, aziridine,
and Michael-acceptor findings and adds a deterministic bifunctional-electrophile
heuristic. It pairs deduplicated, non-overlapping motif sites, reports the
minimum heavy-atom bond distance, and caps output at 64 site pairs. The result
is explicitly structural triage only: it is not a genotoxicity classifier,
biological prediction, 3D geometry estimate, or cross-engine validation.

`validation/genotox_structural_fixtures.json` records three source-referenced
PubChem structures (Mitomycin C, Aflatoxin B1 8,9-epoxide, and anti-BPDE) with
retrieved isomeric SMILES and expected structural motifs. The test verifies
that the checked-in structures parse and exercise the expected motif rules;
the manifest carries provenance rather than a biological ground-truth label.

## Issue #303 evidence

- `cargo test -p chematic-chem --lib genotox --offline` — 6 tests passed,
  including all three source-referenced structures and the bifunctional case.
- The broader genotoxicity predictor, additional structural categories,
  biological fixtures, and licensing review remain outside this bounded slice.

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
