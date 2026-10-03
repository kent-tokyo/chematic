# Reaction 83: post-merge source and binding gate

This is **source evidence**, not a v1.0.31 published-package rerun. On main
`4b76ff9282553549b75bb39e1b86530f450756af`, the pinned 57+26 reaction
fixtures were rerun with RDKit 2026.03.6. The legacy and added-fixture SHA-256
values remain `65b446de1ce66b7d9c355cf8b6bbbfe61d670c955e02931e8e8b8e564c558d3b`
and `4e4ce3ebcea4bcea5398f109554074d97e9f252293fa906493ae3a61d40725a1`.

| Checked Rust result | Rows |
| --- | ---: |
| Product graph, atom origin and template map all agree | 76 |
| Typed unsupported (reactant tetrahedral semantics) | 3 |
| Diagnosed product-valence refusal | 1 |
| Invalid in both engines | 3 |

The generated plain reaction rows, map rows and RDKit comparison report have
SHA-256 values `15d6a3eeaac782326b6afa1ee9d2af816a63db79d0815012b6979682f0661851`,
`2da4e45ac30c8be6c2b2a686fee7a284b6a9f92a99a704f165c2d2cb21c120dc`,
and `a310e2aa05a0b998be7e9278fe7a5f0b2d0e7fb24304f9d93e38cbd49899c633`.
The commands in the [first source audit](2026-10-03-source-reaction-83-checked-graph-origin-map.md)
reproduce them with the label changed to
`chematic 1.0.31 source 4b76ff92 checked`.

The WASM binding's `run_reactants_checked` now has a regression test over all
83 fixtures, comparing normalized product graphs with the independent pinned
RDKit oracle and checking the seven refused/unsupported/invalid statuses.
It passes 76 graph comparisons; it **does not** expose atom origin or template
map labels. The full WASM crate suite passes 359 tests on native host.

The CI source-wheel lanes now run the same 83-row checked-profile graph/status
gate on Linux and macOS after wheel installation and attach the raw report.
Those wheel outcomes are **pending CI**, not evidence in this record. The
published v1.0.31 artifacts have not been replaced or remeasured. Reaction
yield, selectivity, broad SMIRKS equivalence and cross-binding provenance
remain outside this gate.
