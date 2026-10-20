# Refactoring plan

The original structural-debt program is complete for its declared R0-R8
scope. A second pass, R9-R16, is now reducing the remaining checked-in debt
without changing chemistry or public contracts. Implementation history belongs
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

## Current phases

| Phase | Scope | Exit condition |
|---|---|---|
| R9 | Open-PR convergence and baseline pin | Queue empty; latest `main` is the only refactoring base |
| R10 | Test extraction and file-budget reduction | Large inline tests moved out of production files; resolved allowances removed |
| R11 | Python, WASM, CLI, and MCP boundaries | Binding declarations and adapters have bounded domain modules |
| R12 | Chemistry hot spots | Standardization, SMARTS/VF2, and reaction orchestration have staged helpers |
| R13 | I/O and depiction | Format and depiction orchestration are split from policy-free mechanics |
| R14 | Cache and compatibility algorithms | Cache ownership is explicit; generated/declarative algorithms are isolated |
| R15 | Scripts and CI | Repeated harness setup is consolidated without reducing gates |
| R16 | Final ratchet and verification | Debt baseline shrinks, full required checks pass, remaining exceptions are classified |

R9 completed on 2026-10-20 with an empty pull-request queue at `9945f1ed`.
R10 removes the CLI and SVG file allowances by moving their inline unit tests
to dedicated test modules; no production logic or public API changes.
R11 groups CLI dispatch into conversion, molecule, reaction, and batch
boundaries. The public command surface is unchanged, and the former 166-line
dispatcher allowance is removed.
R12 separates fragment scoring, SMARTS primitive families, and RDKit-compatible
reaction stereo application into named stages. The three former function
allowances are removed; chemistry ordering and public results are unchanged.

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
