# Refactoring inventory

This inventory records the final R8 cleanup review. It distinguishes removable
debt from compatibility code that only looks obsolete. Public behavior and
published binding contracts take precedence over reducing line count.

## Removed

- Five unused direct dependencies: the unwired compressed-pattern dependency,
  three unused `chematic-ewald` dependencies, and `wasm-bindgen-test` in a
  crate with no Rust wasm-bindgen-test suite.
- The unreferenced generated compressed SMARTS table. The public
  `compress-patterns` feature name remains as a no-op so downstream manifests
  do not stop resolving.
- An uncalled private Ewald cutoff implementation hidden behind
  `#[allow(dead_code)]`.
- Duplicate Ewald cell-list traversal. Direct and damped cutoff energy now
  share one private pair iterator while retaining their existing fallback and
  accumulation rules.

The lockfile no longer carries the dependency closure used only by those
declarations.

## Retained compatibility paths

- MCP legacy protocol handling is an advertised wire contract and is exercised
  by protocol and schema tests.
- Legacy and accurate CIP modes, named standardization catalogs, and legacy 3D
  entry points remain public compatibility choices with regression coverage.
- Generated chemistry tables and ordered rule registries remain data-shaped
  exceptions. Their order is observable chemistry behavior, so mechanical
  deduplication or sorting is unsafe.
- The duplicate deterministic-math traits in `chematic-ff` and `chematic-3d`
  remain crate-private. Sharing them would require a new public cross-crate
  support API or a new workspace crate; neither is justified by the small,
  stable implementation.

## Audit method

- Scanned every crate manifest against Rust source, tests, examples, benches,
  build scripts, target-specific dependencies, and feature references.
- Searched production Rust for explicit dead-code and unused suppressions, then
  confirmed references before removal.
- Ran the refactoring budget ratchet after deleting the stale generated-file
  exclusion.
- Verified affected crates with their tests, all-feature/all-target builds,
  clippy, formatting, and workspace checks described in the main plan.

Future removals should update this inventory only when a compatibility path is
actually retired through the normal deprecation policy.
