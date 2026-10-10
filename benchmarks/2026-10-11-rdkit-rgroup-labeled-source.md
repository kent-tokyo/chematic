# Label-stable R-group decomposition source check

Date: 2026-10-11  
CheMatic artifact: local source candidate based on `origin/main` `3bd9bba1`  
Oracle: RDKit Python `2025.09.3`

## Scope

This is a bounded functional check, not a claim of full RDKit
`RGroupDecomposition` parity. It covers two terminal mapped wildcard sites on a
single core and verifies that the row labels retain the same substituent
meaning as RDKit.

Core:

```text
c1cc([*:1])ccc1[*:2]
```

Inputs:

```text
Cc1ccc(CC)cc1
CCc1ccc(N)cc1
```

## Result

Both engines returned no unmatched rows and assigned:

| Input | R1 | R2 |
|---|---|---|
| `Cc1ccc(CC)cc1` | ethyl | methyl |
| `CCc1ccc(N)cc1` | ethyl | amino |

The canonical core strings place the two wildcard tokens in a different text
order on this symmetric para-disubstituted ring, but the labelled attachment
semantics and substituent graphs agree. CheMatic also reproduced the same row
after an alternative input atom ordering in its Rust and Python regression
tests.

## Deliberate boundary

The source candidate supports one core with mapped terminal wildcards or mapped
core atoms. It does not implement RDKit's multi-core MCS alignment,
tautomer/core enumeration, non-terminal or multiple R-groups per site, or GA
scoring. Ambiguous label contracts and the bounded match limit are typed
errors. Published-package validation remains pending until a release contains
this source candidate.
