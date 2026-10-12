# Open-work ledger

Updated 2026-10-13 for **v1.1.1**. This file lists unresolved work only.
Completed implementation details belong in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md);
measurements belong in the [validation report](validation.md) and
[benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks).

## Active follow-up

Issues #734 and #754 are closed. Their reaction, SMARTS, autoconf, and fuzz
corpora remain recurring regression gates. Issues #739, #769, #779, and #785
are also closed; implementation and verification records remain linked from
the changelog and validation pages.

- **v1.1.x artifact packet:** the v1.1.1 expectation now includes the 4,072-row
  HBD/QED stress lane and CI compares A6 with the current v1.1.1 wheel. Run the
  full hosted Linux/macOS/Windows published-wheel workflow before replacing
  the current broad v1.0.42 package packet.
- **Equivalent-work performance:** close prepared-operation deficits without
  moving work into parsing; add second-host and peak-RSS records.
- **File interchange:** expand the initial V3000/MOL2/CML/CDXML semantic gate
  to stratified corpora before broad Open Babel claims.
- **CDK/Indigo comparison:** the new source workflow packet covers bounded
  reaction, parent-selection, and labelled R-group lanes with full mismatch
  accounting and 21 alternating timing blocks. Expand beyond 15/39/16 rows,
  add prepared-object and memory lanes plus a second host, then rerun exact
  published artifacts. CDK still has only the earlier smoke gate.
- **R-group contract:** maintain the bounded terminal-wildcard subset across
  Rust, Python, and WASM; broader multi-core and enumeration behavior remains
  outside the contract.
- **A6 residuals:** RDKit 2026.09.1 topology bounds match within `1e-12` across
  all 768,664 frozen-corpus cells. The conformer-tail candidate fixes a real
  gradient early-exit bug, aligns relevant-ring traversal, and exposes an
  optional selective rescue profile. Once-only v3 results are 94/100 within
  1 kcal/mol at essentially tied but 0.27% slower median time; rescue is
  95/100 and 2.1% slower. A6 plus holdouts v1/v2/v3 are now development data.
  A local MMFF distance cache now funds a bounded one-candidate/100-iteration
  retry profile that was faster than RDKit in all three development blocks and
  improved the within-1-kcal/mol count to 261/265. This is not final evidence.
  Use a fresh v4 set after candidate freeze, and keep atom typing,
  convergence, stereochemistry, clashes, same-coordinate energy, basin
  diversity, minimization tails, time and memory separate.
- **Independent workflow corpora:** expand beyond the current bounded lanes
  with multi-reactant reactions, multi-core R-group cases, and real MOL2/CJSON
  files. Keep refusals and semantic-loss reports in the denominator.

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
| A1 Perception and descriptors | Active |
| A2 Stereo and identity | Active |
| A3 Fingerprints and retrieval | Maintenance |
| A4 Workflows and interchange | Active |
| A5 Independent adjudication | External |
| A6 3D and force fields | Experimental |

## Completion rule

For every gate, record the support domain, versions, corpus, options, failure
policy, and counts for success, failure, refusal, and skipped inputs. A refusal
is not a match, a local pass is not a published-package result, and registry
availability is not chemistry validation.
