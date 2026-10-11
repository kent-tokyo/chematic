# Refactoring plan

The structural-debt program is complete for its declared R0-R8 scope. This
page records the guardrails that remain active; implementation history belongs
in the changelog and [inventory](refactoring-inventory.md).

## Active guardrails

- `scripts/check_refactor_budget.py` scans handwritten production Rust code.
- `validation/refactor_budget_v1.json` is a debt ratchet: existing exceptions
  may shrink but may not grow.
- New production files are limited to 1,500 lines. New functions are limited
  to 120 lines and 25 lexical decision points.
- Generated and declarative tables require narrow, documented exclusions.
- Public paths, struct-literal compatibility, typed errors, Python names,
  TypeScript names, and JSON shapes stay stable during structural work.

The scanner is not a correctness oracle. Chemistry regressions, cross-binding
contracts, output digests, and paired performance gates remain authoritative.

## Completed phases

| Phase | Scope | Result |
|---|---|---|
| R0 | Budget scanner, CI ratchet, baselines | Complete |
| R1 | Registries, generated data, test separation | Complete |
| R2 | Reaction templates, matching, products, diagnostics | Complete |
| R3 | Molecular format parsers and writers | Complete |
| R4 | Python and WASM bindings | Complete |
| R5 | Descriptors and tautomer logic | Complete |
| R6 | Canonicalization, rings, aromaticity | Complete |
| R7 | MMFF/UFF and 3D pipeline | Complete |
| R8 | Dead code, dependency, and helper audit | Complete |

Notable outcomes include a bounded reaction module tree, separated
lexical/semantic format stages, API-domain binding modules, bounded descriptor
and tautomer helpers, split canonical/ring/aromaticity algorithms, staged
force-field and 3D preparation, and removal of unused dependencies and dead
private implementations.

## Required checks for future refactors

1. Run formatting, affected tests, and clippy with warnings denied.
2. Run the refactor-budget checker and its scanner tests.
3. Run relevant Rust, Python, WASM/Node, and cross-binding gates.
4. Compare outputs for every touched chemistry or format domain.
5. Measure paired time/instructions before accepting a hot-path change.

Moving code alone is not sufficient. A refactor is complete only when public
contracts remain stable, applicable gates pass, and any retained exception is
documented as generated or priority-sensitive declarative data.
