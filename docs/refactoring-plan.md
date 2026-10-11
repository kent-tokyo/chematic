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
| R2 | Reaction template, matching, product, sanitization, provenance, compatibility | Complete |
| R3 | File parsers: lexical decode, semantic validation, builders, diagnostics | Complete |
| R4 | Python and WASM bindings by API domain | Complete |
| R5 | Descriptors and tautomer rules/search/scoring | Complete |
| R6 | Canonicalization, ring perception, and aromaticity | Complete |
| R7 | MMFF/UFF and 3D preparation, optimization, and quality policy | Active |
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
component parsing are isolated in `transform/templates.rs`. Matching,
sanitization, stereo transfer, product construction, and prepared-reaction
orchestration now have separate modules; `transform.rs` is reduced from 8,282
lines before R1 to 528 lines of public API and shared state.
The four remaining over-budget reaction functions were then decomposed into
bounded input validation, query matching, E/Z transfer, aromaticity/Kekule/
valence diagnostics, atom assembly, and product finalization helpers. R2 now
has no refactoring-budget exceptions.

R3 starts by moving embedded format tests for CDXML, CIF, LAMMPS data, V2000,
V3000, ORCA, QCSchema, semantic documents, and XYZ into child test modules.
This removes five parser file allowances immediately and leaves the production
files organized for lexical/semantic/builder extraction.

R3 finishes with no `chematic-mol` file or function exceptions. Format codecs
are split into bounded lexical, semantic-validation, construction, and output
stages, including CDXML, KET, MOL2, LAMMPS, semantic documents, and volumetric
streams. The complete `chematic-mol` all-target test suite, clippy with warnings
denied, the refactor-budget scanner, and formatting/diff checks pass. A KET
contract regression found by the full gate was corrected so flat malformed
documents continue to report the original missing `atoms` field.

R4 finishes with no Python or WASM binding exceptions. Large WASM modules now
isolate contract tests, LAMMPS and typed-array APIs, common format conversion,
and molecule workflows. The 6,087-line Python `Mol` implementation is split
into eight API-domain modules; bulk descriptors, reports, reactions, LAMMPS,
and pipeline conversion each have separate compute and serialization stages.
Exported JavaScript/Python names and JSON shapes are unchanged. All 375 native
WASM tests, the release Node build and Node integration tests pass. The Python
wheel builds, 1,109 of 1,110 tests pass (the sole failure is the pre-existing
caffeine HBA oracle mismatch), and a representative method/descriptor/bulk
snapshot is byte-identical to the published 1.1.0 wheel. Clippy and the budget
scanner pass.

R5 finishes with no descriptor, topology-descriptor, SA-score, or tautomer
file/function exceptions. Descriptor code is separated into mass/basic, TPSA,
Crippen, ring/filter, MQN, topology/shape, composition, and information/BCUT
domains. Tautomer rules, aromatic transforms, graph transforms, search, and
audited parent selection are isolated, while the existing public module paths
remain unchanged. The TPSA nitrogen table is data-shaped instead of a large
guard chain; MQN accumulation and tautomer matching/parent search use bounded
helpers. All 938 `chematic-chem` library tests pass, including descriptor and
tautomer parity regressions. Clippy with warnings denied, formatting, the
budget scanner, and all six scanner unit tests pass.

R6 finishes with no canonicalization, Kekule assignment, SSSR, aromaticity,
or RDKit-parity file/function exceptions. Molecule storage and graph/cache
operations, Kekule matching, aromaticity tracing and ring families, RDKit
canonical ranking, parity typing, SSSR basis/candidate/symmetrized expansion,
and canonical-SMILES E/Z planning, traversal, and stereo remapping now have
bounded modules and helpers. Public paths and canonical output remain stable.
The affected crates pass 774 library tests (with two documented ignored
profiles), clippy with warnings denied, formatting, the budget scanner, and
all six scanner unit tests. R6 removes seven file and seventeen function
allowances from the debt baseline.

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
