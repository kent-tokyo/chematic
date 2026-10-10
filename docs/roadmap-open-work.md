# Open-work ledger

Updated 2026-10-10 for **v1.0.42**. This file lists unresolved work only.
Completed implementation details belong in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md);
measurements belong in the [validation report](validation.md) and
[benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks).

## Open issues

| Priority | Issue | Remaining work | Exit |
|---:|---|---|---|
| 1 | [#734 — SMIRKS differences](https://github.com/kent-tokyo/chematic/issues/734) | Broad reaction-rule validation, including the contributor corpus and classified fuzz residuals | All pinned rows classified; no wrong-confident supported result; public artifacts rerun |
| 2 | [#754 — xsmarts-autoconf](https://github.com/kent-tokyo/chematic/issues/754) | Upstream probe contribution and continued dialect/fuzz accounting | Minimized probes submitted; remaining differences documented as specification, policy, unsupported, or defect |
| 3 | [#784 — WASM binding parity](https://github.com/kent-tokyo/chematic/issues/784) | Expose the selected v1.0.40+ interoperability APIs in WASM/Node | Typed API, bounds, Node tests, and Rust/Python/WASM fixture parity |
| 4 | [#786 — reaction rejection diagnostics](https://github.com/kent-tokyo/chematic/issues/786) | Design opt-in product-rejection reasons and mapped-atom reporting | Stable typed schema, preserved indices, no fast-path or accepted-product regression |

Issues #739, #769, #779, and #785 are closed. Their implementation and
verification records remain linked from the changelog and validation pages.

## Evidence refreshes

These are recurring gates, not separate feature requests:

- **Published chemistry:** rerun the reaction, SMARTS, CIP, fingerprint, and
  3D gates from exact PyPI/npm/crates.io artifacts on their declared platforms.
- **Performance:** compare equal outputs and equal work in alternating repeated
  blocks; split parsing, perception, prepared reuse, and memory.
- **New RDKit releases:** pin distributed artifacts and hashes. Keep source,
  package, and browser lanes separate.
- **Representation limits:** retain strict/reporting paths for lossy MOL,
  V3000, CDXML, Markush/polymer, Standard InChI, and canonical identity.

## Dependency classes

- **`local-open`:** implementation, regression tests, documentation, and source gates.
- **`toolchain-open`:** cross-platform wheels, npm/WASM builds, and browser runners.
- **`external-open`:** new comparator packages, upstream review, and independent adjudication.
- **`data-sealed`:** one-time evaluation only; inspected rows never become sealed again.
- **`historical`:** retained provenance that does not establish current behavior.

## Accuracy packages

| Package | State |
|---|---|
| A0 Evaluation contract | Complete |
| A1 Perception and descriptors | Open |
| A2 Stereo and identity | Active |
| A3 Fingerprints and retrieval | Open |
| A4 Workflows and interchange | Active |
| A5 Independent adjudication | External |
| A6 3D and force fields | Active |

## Completion rule

For every gate, record the support domain, versions, corpus, options, failure
policy, and counts for success, failure, refusal, and skipped inputs. A refusal
is not a match, a local pass is not a published-package result, and registry
availability is not chemistry validation.
