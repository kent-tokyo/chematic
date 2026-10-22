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

## Completed second pass

| Phase | Scope | Result |
|---|---|---|
| R9 | Open-PR convergence and baseline pin | Complete |
| R10 | Test extraction and file-budget reduction | Complete |
| R11 | Python, WASM, CLI, and MCP boundaries | Complete |
| R12 | Chemistry hot spots | Complete |
| R13 | I/O and depiction | Complete |
| R14 | Cache and compatibility algorithms | Complete |
| R15 | Scripts and CI | Complete |
| R16 | Final ratchet and verification | Complete |

R9 completed on 2026-10-20 with an empty pull-request queue at `9945f1ed`.
R10 removes the CLI and SVG file allowances by moving their inline unit tests
to dedicated test modules; no production logic or public API changes.
R11 groups CLI dispatch into conversion, molecule, reaction, and batch
boundaries. The public command surface is unchanged, and the former 166-line
dispatcher allowance is removed.
R12 separates fragment scoring, SMARTS primitive families, and RDKit-compatible
reaction stereo application into named stages. The three former function
allowances are removed; chemistry ordering and public results are unchanged.
R13 separates CJSON conversion policy from loss reporting and EPS orchestration
from viewport, highlight, bond, and label emission. The two former function
allowances are removed; conversion reports and EPS byte ordering are unchanged.
R14 gives the RDKit ring-compatibility search one explicit state owner for its
candidate cache, budget, and fallback stages, and separates the molecule-owned
RDKit-model cache gate from per-atom compatibility rules. Two more function
allowances are removed without changing parity or budget-exhaustion results.
R15 consolidates JSON, compressed-row, digest, and failure handling shared by
four RDKit evidence gates. Every original checker still runs independently;
the refactor removes duplicated harness mechanics without weakening a gate.
R16 closes the ratchet at 14 file and 76 function allowances, classifies the
remaining debt in the inventory, and verifies the stacked result with the
workspace, binding, budget, and hosted-CI gates.

The 2026-10-21 maintenance pass retires one unreferenced exploratory probe and
places 13 retained evidence scripts under diagnostics or validation. It does
not change production chemistry or the R16 debt-ratchet baseline; see the
[script inventory](script-inventory-2026-10-20.md) for paths and retention
reasons.

## Active third pass

R17-R25 finishes the remaining structural debt without broadening chemistry
scope. Every phase must reduce or preserve the ratchet; raising an allowance is
not an acceptable completion condition.

| Phase | Scope | Exit condition | Status |
|---|---|---|---|
| R17 | Converge the open format-I/O pull-request stacks | Parent-to-child CI is green, merged state is reproducible, baseline is frozen | Complete |
| R18 | MOL2/CJSON checked-conversion orchestration | Candidate generation, loss accounting and verification are separate; two function allowances removed | Complete |
| R19 | Standardization and chemistry transforms | Long orchestration functions split with unchanged typed outcomes | Ready for review |
| R20 | Depiction and layout | Layout stages bounded; image and coordinate fixtures unchanged | Ready for review |
| R21 | Search, SMARTS and reaction internals | Search state and reaction projection stages have explicit owners | Ready for review |
| R22 | Compatibility-heavy RDKit tables and adapters | Declarative data separated from control flow; parity fixtures unchanged | Ready for review |
| R23 | Python, WASM, CLI and MCP boundaries | Binding registration and transport glue remain thin and cross-checked | Planned |
| R24 | Validation scripts, CI and documentation | Shared harnesses replace duplication without weakening gates | Planned |
| R25 | Final ratchet and release-grade verification | At most 8 file and 50 function allowances; full applicable gates pass | Planned |

R18 reduced the checked function allowances from 76 to 74. Its hosted CI,
Python/WASM bindings, security audit, CodeQL and format-contract suites passed
before merge.

The R19 candidate moves the auditable standardization pipeline and its unit
tests behind dedicated module boundaries. `standardize.rs` falls from 3,446 to
1,329 lines, removing its file allowance and reducing the checked baseline to
13 files and 74 functions. The Rust unit tests plus the standardization
idempotency, metal-disconnection, stereo-preservation and shared-contract suites
must preserve stage order, warnings and typed outcomes.

The stacked R20 candidate assigns ring placement, acyclic growth and clash
relief to separate modules, leaving the public layout facade at 397 lines.
Two one-stage placement functions are brought below the function limits, so the
ratchet reaches 12 file and 72 function allowances. The complete depiction
crate (107 unit tests, 3 public stress tests and doctests) preserves the
coordinate, crossing, SVG, PNG and EPS contracts.

The stacked R21 candidate separates reaction-query, VF2, MCS and SMARTS-parser
tests from production ownership and gives the VF2 size/element prefilter its
own module. The four production files are each below 1,500 lines, reducing the
ratchet to 8 file and 72 function allowances without changing public search or
reaction behavior. The complete `chematic-smarts` and `chematic-rxn` suites,
their doctests, and warnings-as-errors clippy checks pass.

The stacked R22 candidate separates ordered MMFF torsion cases, InChI cleanup
families, PDB residue bond rules, implicit-valence branches, ring-stereo
bookkeeping and depiction merge operations. Nine compatibility-port functions
now fit the normal limits, reducing the ratchet to 8 file and 63 function
allowances. Full `chematic-smiles` and `chematic-ff` tests preserve the pinned
PDB, InChI, stereo, depiction and force-field contracts.

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
