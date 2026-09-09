# v1.0.10 atomic-number SMARTS H1 bridge

The reaction-template normalization bridge now accepts the common
explicit-single-hydrogen forms `[#N;H1]` and `[#N;H1:map]` in addition to
`[#N]` and `[#N:map]`. It emits deterministic aliphatic/aromatic variants and
preserves the atom map. Unsupported compound SMARTS remain fail-closed;
full query-aware application semantics are not claimed by this bounded gate.

Reproduction:

```text
cargo test -p chematic-rxn --offline --test cross_binding_contract
cargo test -p chematic-rxn --offline --lib atomic_number_smirks
```

Evidence: 2 shared Rust contract tests, 3 focused Rust library tests, and the
full 209-test `chematic-rxn` package run passed with zero failures. This is a
bounded compatibility slice; query-aware SMARTS matching and broader primitive
coverage remain open.

Machine-readable evidence: [`2026-09-09-rxn-atomic-number-h1-v1.0.10.json`](2026-09-09-rxn-atomic-number-h1-v1.0.10.json).
