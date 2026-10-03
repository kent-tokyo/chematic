# Open-work ledger

Updated 2026-10-03 for the published **v1.0.31** release. The six release
channels are verified; chemistry/benchmark results below are still from
v1.0.30. The
[roadmap](https://github.com/kent-tokyo/chematic/blob/main/ROADMAP.md) sets priority and exit criteria; this ledger names
dependencies. Completed work and raw evidence stay in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md), [validation report](validation.md), and
[benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md).

## Queue

| Order | Open work | Next verifiable exit |
|---:|---|---|
| 1 · P1/A4 | Published v1.0.30 reaction Rust audit: 73/83 graph+origin+map matches; ten other rows = one map-only, four graph, one origin, one refusal, three jointly invalid. The checked-source semantic profile is not a package result. | Adjudicate each axis, fix supported defects or return typed unsupported, add row regressions, then rebuild and rerun Rust/Python/npm artifacts. Preserve 57 legacy fixtures and verify npm E/Z JSON. Do not count refusals or jointly invalid rows as matches. |
| 2 · P1/A4 | 200/310,000 published SMARTS match-set cells differ: 194 symmetrized-ring and six organometallic semantics; 43 change the hit/no-hit Boolean. | Triage `[R2]` 64, `[R3]` 62, `[R1]` 59 cells first; classify spec/unsupported/bug, preserve native SSSR and rerun full affected published matrices with separate match-set/Boolean deltas. |
| 3 · P2/A6 | MMFF94/3D remains Experimental; the 265-row quality result is historical and 2,908 heavy-atom type differences remain in a separate source lane. | Separate typing/parameters, identical-coordinate energy/gradient, convergence/timeout, stereo/clashes and independent conformer quality. Rebuild current artifacts and require typed failures before any speed promotion. |
| 4 · P0.2 | Extend paired performance beyond the completed Python 63-operation and Node/Chromium/Firefox/WebKit Morgan lanes. | Equal-output, equal-work ≥20-block intervals for more operations/corpora; distinguish parsing, perception, prepared reuse and library memory from process RSS. The historical “21 faster” count is not a current package gate. |
| 5 · P1/A2 | CIP has 9,995 exact rows and five typed abstentions on the pinned exposed 10k lane. | Keep abstentions stable under reorder/round trips; adjudicate phosphorus and define a lone-pair convention before assigning labels. |
| 6 · P2 | CDXML, Markush/polymer, Standard InChI and canonical identity are bounded. | Per-format preservation, semantic round-trip, typed-refusal and no-false-merge gates across bindings. |
| 7 · external | A later RDKit oracle is not pinned here. | Rebaseline only against a new official artifact; preserve 2026.03.6 results separately. |

P0.1's **published-artifact accounting is complete**, not strict RDKit
parity: [the acceptance policy](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md)
retains the 200 SMARTS failures, five CIP abstentions and four npm API gaps.
The [published-artifact packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md)
contains the versioned denominators and raw-row identities.

## Accuracy package disposition

| Package | State | Queue link |
|---|---|---|
| A0 Evaluation contract | Complete | Preserve frozen/exposed-cohort accounting |
| A1 Perception and descriptors | Open | Representation and descriptor boundaries |
| A2 Stereo and identity | Active | Queue 5: CIP abstentions and adjudication |
| A3 Fingerprints and retrieval | Open | Separate compatible-Morgan/holdout gates |
| A4 Workflows and interchange | Open | Queues 1–2: reactions and SMARTS |
| A5 Independent adjudication | External | Gold data and non-maintainer review |
| A6 3D and force fields | Active | Queue 3: quality before speed |

## Cross-cutting dependencies

- **`local-open`:** parser security and batch-result accounting, browser cancellation
  and typed resource limits, stereo regression fixtures, release-document
  synchronization. These can proceed without a new oracle.
- **`toolchain-open`:** rebuilt Python/npm/WASM artifacts, browser lanes and
  cross-platform builds need their named runtimes; a source-only fix is not a
  published-package result.
- **`data-sealed`:** freeze a candidate and audit overlap before one-time
  evaluation. Exposed or inspected rows never become sealed again.
- **`external-open`:** registry publication, official future RDKit releases,
  non-maintainer review and independent security assessment require their
  respective outside artifacts or people.
- **`historical`:** source-only A/B timings, rejected candidates and older
  release-channel records remain available for provenance, not current claims.

## Completion rule

For every gate, declare the supported domain, comparator, versions, corpus,
options and failure policy. Account for success, failure, refusal and skipped
inputs; test the corrected behavior; rebuild affected bindings; and link
committed evidence. Name source, published-package and public-channel state
separately. A refusal is not an exact match, and a passing local check is not
a verified release.
