# Ordered-ring and V3000 parsing follow-ups — 2026-10-16

Two output-preserving implementation changes were accepted. A third prototype,
lazy initialization of ring-search scratch arrays, was rejected because its
paired confidence intervals included no change and is not present in the
candidate.

## Accepted changes

| Change | Fixed-corpus result | Correctness boundary |
|---|---|---|
| Ordered ring reconstruction | ChEMBL first-use ring count: 1.020x baseline/candidate, 95% CI [1.004, 1.024], 16/20 blocks faster. NCI: 1.006x [0.998, 1.015], so no NCI speed claim. | Frozen old/new searches matched 71,680 small-graph queries, 935,863 pinned-corpus queries and 116 long-cycle/chord/parallel-edge queries, including exact ring order. |
| V3000 absent/owned bond attributes | Parse improved 1.068x [1.046, 1.103] and 1.084x [1.059, 1.111] on two 1,000-record slices. Round trip improved 1.056x [1.021, 1.068] and 1.072x [1.059, 1.096]. | Full molecular reports, coordinate bits, diagnostics, errors and rewrite bytes matched the frozen implementation on 2,064 records; 23,328 endpoint-property strings also matched. |

Ratios are baseline time divided by candidate time. Each lane used fresh-process
rotating blocks on one macOS arm64 host. The ring lane used 20 blocks; V3000
used 21. Inputs were fixed slices from the repository's pinned corpora. The
measurements are source-candidate evidence, not published-package,
cross-platform or toolkit-wide superiority claims.

The V3000/Open Babel comparison was ranked only for rows where both engines
preserved the declared graph and metadata. On those selected cohorts, all six
parse/write/round-trip lanes favored the candidate in every block. The original
unfiltered diverse corpora contain Open Babel semantic differences and remain
unranked. Controlled attribute-free chains passed all three comparator gates.

## Reproduction and retention

The original measurements used source base
`235df7a14ea2d11724ccf329fef6675f67a76423`, CPython 3.13, RDKit 2026.09.1,
COSMolKit 0.5.0rc22 and Open Babel 3.2.1. The accepted code and independent
regression oracles are retained in the repository. Large transient wheels,
native binaries, duplicated source trees, copied corpora and build logs are not
versioned; they are reproducible build products rather than source artifacts.

The rejected lazy-scratch prototype measured 1.016x [0.955, 1.082] and 1.012x
[0.987, 1.042] for first-use ring count on the two slices. It was removed before
this candidate because neither interval established a gain.
