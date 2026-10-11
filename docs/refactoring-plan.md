# Refactoring plan

This track reduces structural debt without changing chemistry, public APIs, or
published compatibility claims. Refactoring and behavioral optimization are
kept in separate changes.

## Guardrails

- `scripts/check_refactor_budget.py` scans handwritten production Rust code.
- The checked-in baseline is `validation/refactor_budget_v1.json`.
- New production files must stay at or below 1,500 lines. New functions must
  stay at or below 120 lines and 25 lexical decision points.
- Existing exceptions may shrink but may not grow. When an exception reaches
  the normal limit, its baseline entry must be removed.
- Generated tables and dedicated test files are excluded explicitly rather
  than hidden in broad path globs.
- Public paths, struct-literal compatibility, typed errors, Python names,
  TypeScript names, and JSON shapes remain stable throughout this track.

The scanner is a debt ratchet, not a correctness oracle. Chemistry tests,
cross-binding gates, output digests, and paired performance measurements remain
authoritative.

## Phases

| Phase | Scope | State |
|---|---|---|
| R0 | Refactor budget, CI ratchet, regression baselines | Complete |
| R1 | Declarative registries, generated data, and test-module separation | Complete |
| R2 | Reaction template, matching, product, sanitization, provenance, compatibility | Active |
| R3 | File parsers: lexical decode, semantic validation, builders, diagnostics | Planned |
| R4 | Python and WASM bindings by API domain | Planned |
| R5 | Descriptors and tautomer rules/search/scoring | Planned |
| R6 | Canonicalization, ring perception, and aromaticity | Planned |
| R7 | MMFF/UFF and 3D preparation, optimization, and quality policy | Planned |
| R8 | Legacy/dead-code inventory, exact helper deduplication, dependency cleanup | Planned |

R1 started by separating the MCP registry schemas from tool execution. The
former 580-line `list_tools` function is now a short ordered registry assembled
from one bounded schema function per tool. Runtime order and schema content are
covered by the existing MCP schema and protocol tests. The 3,335 lines of
reaction unit and input-boundary tests formerly embedded in `transform.rs` are
also isolated in child test modules; private-helper access and all 293 reaction
tests remain unchanged. Generated-table exclusions are now accepted only for
existing crate source files whose headers identify them as generated, preventing
the exclusion list from becoming an unchecked escape hatch.

R2 begins with the stable diagnostic boundary: transformation errors, limits,
accounting, sanitization diagnostics, compatibility refusals, and variant
diagnostics now live in `transform/diagnostics.rs` while their existing
`transform::*` public paths remain re-exported unchanged. SMIRKS normalization,
reactant-query compilation, product-template specifications, and grouped
component parsing are isolated in `transform/templates.rs`.

## Required checks per phase

1. `cargo fmt --all -- --check`, affected crate tests, and affected-crate
   clippy with warnings denied.
2. The refactor budget checker and its scanner unit tests.
3. Relevant Rust, Python, and WASM/Node contract gates.
4. Output equality for canonical SMILES, SMARTS, reactions, descriptors,
   fingerprints, formats, and 3D when those areas are touched.
5. Paired instruction/time measurements before accepting a hot-path change.

Moving code is not enough to close a phase. A phase closes only when the
remaining over-budget entries in its declared scope have either been split or
documented as generated/declarative exceptions and all applicable gates pass.
