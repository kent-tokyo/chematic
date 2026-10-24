# MCS query SMARTS in Python and WASM

Source candidate: `codex/mcs-smarts-output-20261025`, based on `origin/main`
`a73204f91efe9df5c772a28fa95237c088fcc9e9`. Comparator: RDKit 2026.09.1
in the existing [five-case development oracle](2026-10-25-rdkit-mcs-query-alternatives.md).
This is a source/API regression check, not an independent FMCS accuracy score.

## Change

The Rust MCS query retained observed atom and bond alternatives, but Python
`find_mcs_checked()` and WASM `mcs_smiles_json_with_config()` returned a
concrete first-input projection. A user could therefore receive `CO` for a
`CO`/`NO` common query that should match both molecules.

`mcs_query_to_smarts()` serializes the MCS query subset, including branches,
rings, disconnected components, and explicit atom/bond OR expressions. Python
`find_mcs_smarts_checked()` returns `(smarts_or_none, was_timed_out)`; a resource
limit raises `ChematicInputError`. WASM `mcs_smarts_json_with_config()` returns
`smarts`, `wasTimedOut`, and `wasResourceLimited`. Existing concrete-Mol and
SMILES functions retain their behavior.

The SMARTS APIs exclude or reject charge, isotope, and chirality matching
options because those constraints are not represented in the current MCS
query. The serializer rejects unsupported predicates and malformed graphs
instead of silently weakening them.

## Verification

- Rust serializer: round-trip parse and source re-match for relaxed atom,
  relaxed bond, and ring cases; unobserved atom and bond alternatives remain
  excluded. Three targeted unit tests pass.
- WASM Rust API: two relaxed-comparator cases retain OR expressions and
  exhaustive status. One targeted unit test passes.
- `cargo check -p chematic-py -p chematic-wasm --lib` passes.
- Clippy with warnings denied, `cargo fmt --check`, and the benchmark index
  consistency check pass.
- Python and Node runtime regressions are committed for CI. The local Python
  Rust test build stopped during dependency compilation when the isolated
  worktree ran out of disk space; it did not reach the tests. Its generated
  `target` directory was removed afterward. A rebuilt extension/WASM package
  was not executed locally at this stage.

Exact RDKit SMARTS text identity, charge/isotope/chirality query encoding,
the frozen MCS holdout, and broad RDKit superiority are not established here.
