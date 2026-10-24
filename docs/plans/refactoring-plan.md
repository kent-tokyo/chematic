# Refactoring closeout and maintenance rules

The R0-R28 structural-debt program is complete. Implementation history and
the reasons for retained exceptions are in the
[refactoring inventory](../reference/refactoring-inventory.md); user-visible changes belong
in the changelog.

## Current guardrails

- `scripts/checks/check_refactor_budget.py` scans handwritten production Rust code.
- `validation/refactor_budget_v1.json` is a ratchet: existing allowances may
  shrink but may not grow.
- New production files are limited to 1,500 lines. New functions are limited
  to 120 lines and 25 lexical decision points.
- Generated and declarative tables require narrow, documented exclusions.
- Public Rust paths, struct-literal compatibility, typed errors, Python names,
  TypeScript names, and JSON shapes stay stable during structural work.

The scanner is not a correctness oracle. Chemistry regressions, cross-binding
contracts, output digests, and paired performance gates remain authoritative.

## Completed program

| Phases | Scope | Result |
|---|---|---|
| R0-R8 | Scanner, baseline, core chemistry, formats, bindings, 3D, dead code | Complete |
| R9-R16 | Test extraction, boundary cleanup, hot spots, CI, first closeout | Complete |
| R17 | Format-I/O pull-request convergence and baseline freeze | Complete |
| R18-R20 | Checked conversion, standardization, depiction, and layout | Complete |
| R21-R23 | Search/reaction internals, compatibility ports, bindings, and MCP | Complete |
| R24 | Shared validation and CI harnesses | Complete |
| R25 | Final orchestration splits and release-grade verification | Complete |
| R26 | Low-risk test extraction and bounded parser/CIP ownership | Complete |
| R27 | Structural-alert catalog and test ownership | Complete |
| R28 | InChI dedup and indexed-graph ownership | Complete |

R17-R28 reduced the checked debt baseline from 14 file and 76 function
allowances to **2 file and 50 function allowances**. The program separated:

- standardization stages, reaction projection, SMARTS/VF2/MCS ownership, and
  CIP/canonical/ring orchestration;
- CJSON/MOL2 conversion policy from loss reporting;
- layout, depiction preflight, PNG/EPS rendering, and periodic-neighbor work;
- MMFF, InChI, PDB, valence, ring-stereo, and RDKit-compatibility branches;
- CLI, Python, WASM, MCP transport/schema, and validation-script boundaries.

R26 moved PDB, POSCAR, reaction-SMARTS, SMARTS-parity and SMILES-writer tests
out of production modules, isolated the reaction-SMARTS parser, and gave E/Z
completeness its own CIP module. This closed the PDB, POSCAR, reaction-SMARTS
and CIP file allowances without changing public paths or chemistry behavior.

R27 moved the immutable PAINS and Brenk SMARTS catalogs and their tests out of
the alert runtime module. Pattern text and order, cache construction, matching
budgets, public APIs, and fail-safe behavior are unchanged. This closes the
`alerts.rs` file allowance without treating ordered chemistry data as runtime
logic.

R28 moved typed indexed-graph relation comparison and the dedup test suite out
of the InChI candidate-grouping module. Existing public names are re-exported
from the same `dedup` path; graph strictness, atom-map policy, diagnostic
ordering, native-InChI evidence and feature-gated behavior are unchanged. This
closes the `dedup.rs` file allowance.

R0-R27 passed the full release-grade verification stack before merge. R28 was
validated with the affected crate's default and native-InChI test suites,
warnings-as-errors clippy, formatting, source-preservation checks, the budget
scanner, and diff checks; repository-wide release and binding gates remain the
responsibility of the PR CI before merge.

## Retained debt

The remaining allowances are ordered chemistry data, bounded graph/search
algorithms, format state machines, or compatibility ports whose branch and
accumulation order is observable. They are reviewed debt, not exemptions.

Remove one bounded stage at a time with its domain corpus. Do not split tables
or reorder chemistry merely to reduce line count. A new phase is warranted
only when it can lower the ratchet without changing public behavior.

## Required checks for future refactors

1. Run formatting, affected tests, and clippy with warnings denied.
2. Run the budget checker and scanner tests.
3. Run relevant Rust, Python, WASM/Node, and cross-binding gates.
4. Compare outputs for every touched chemistry or format domain.
5. Measure paired time or instructions before accepting a hot-path change.

Moving code alone is not completion. Public contracts must remain stable,
applicable gates must pass, and every retained exception must remain justified
in the inventory.
