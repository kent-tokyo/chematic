# chematic 1.x Trust Release plan

Updated 2026-09-16. This document translates the priority order in
[`ROADMAP.md`](../ROADMAP.md) into verifiable delivery gates. Historical daily
notes belong in Git history and archived snapshots, not in this plan.

## Goal

Make every public compatibility claim reproducible from a pinned package,
operation, corpus, comparator, and failure policy. The Trust Release does not
promise complete RDKit replacement, full V3000 semantics, or 3D parity.

## Status

| Track | Status | Completed evidence | Remaining exit |
|---|---|---|---|
| T0 — measurement integrity | Baseline complete | Full-row SMARTS accounting and withdrawn partial-result claim | Keep baseline/candidate builds independent and reject incomplete artifacts atomically |
| T1 — compatibility contract | Partial | API profiles, V3000 baseline, frozen candidate, post-freeze attestation, 2k/8k split | Run sealed holdout; publish operation-level results and refusal coverage |
| T2 — release truth | Complete for v1.0.15 | Six-channel record and cross-platform package smoke | Repeat with new artifacts for each release |
| T3 — browser/MCP product | Partial | Published RDKit.js comparison, 10k Worker path, runtime MCP inventory | Complete cancellation, partial-failure, offline, and public-package examples |
| T4 — parser security | Partial | Fixed five-format corpus and Linux isolated CI gate | Add pinned upstream-reproducer lanes and scheduled fuzz/sanitizer evidence |
| T5 — stereo safety | Partial | 300 development structures and 5,115 spelling variants | Add unused challenge set, atom/bond permutations, full format round trips, and absolute adjudication |
| T6 — maintenance/external validation | Open | Advisory intake and dependency gates exist | Independent review, recurring release drill, competitor watch, and separate 3D program |

## T0 — Measurement integrity

Every runner must record expected rows, completed rows, failures, tool versions,
input hashes, output hashes, and whether baseline and candidate were built from
different revisions. Partial output must never be promoted to a result.

The current SMARTS experiment is a useful negative result: the candidate moved
from 12 to 21 residual cells and was not adopted. It remains evidence that the
gate catches regressions, not evidence of improved compatibility.

## T1 — Compatibility contract and sealed evaluation

The non-release annotated tag `trust-eval-candidate-20260916` freezes candidate
commit `c2682e3aa75c21566c86ce1ade9cbd052838c694`. The source used for the split was
acquired after that freeze and is bound to a maintainer unused-data attestation.
The prepared split contains 2,000 development rows and 8,000 sealed rows.

Next execution:

1. verify candidate tag, source, attestation, and split hashes;
2. run the candidate once against the sealed holdout;
3. preserve raw output before computing summaries;
4. report valid, refused, failed, and timed-out rows separately;
5. publish operation-level comparisons without tuning against the holdout.

No sealed score exists until that sequence completes.

Ordinary V3000 support is version-pinned against RDKit 2025.09.3 and Indigo
1.46.0. Coordination chemistry, haptic bonds, polymer expansion, and semantic
interpretation of ENDPTS/ATTACH remain outside the typed contract. Opaque token
retention is not equivalent to editable semantics.

## T2 — Release truth

For every release, verify GitHub tag/release, crates.io, docs.rs, PyPI, npm, and
the published site. Record digest, observed version, timestamp, URL, and runtime
smoke results. A successful build or tag is not publication proof.

Candidate tags must not publish. Immutable published artifacts are never
overwritten; failed channels are retried explicitly and remain visible until
verified.

## T3 — Browser and MCP product gate

The reference workflow is `npm install` → typed import → initialization → Worker
batch → cancellation/error handling → handle cleanup. Required browsers are
Chromium, Firefox, and WebKit; Node is a separate runtime lane.

The current public-package comparison uses a fixed exposed 10,000-row corpus.
It records package bytes, local no-store startup, parse/write, parse-inclusive
fingerprints, coverage/refusals, and bounded memory diagnostics. It does not
claim internet/CDN latency or broad fingerprint superiority.

Remaining work:

- SDF/CSV/SMI partial-failure and retry ordering;
- cancellation and UI-heartbeat budgets;
- operation parity across native scalar, native batch, and Worker paths;
- offline-after-initialization behavior;
- MCP oversized-input and cancellation contracts.

## T4 — Parser security

The fixed corpus covers SMILES, SMARTS, V2000, V3000, and SDF with valid
controls. Release evidence requires Linux process isolation, enforced wall-time
and address-space limits, complete row accounting, and raw artifacts. A macOS
functional run without enforced memory/network boundaries is diagnostic only.

Panics, crashes, signals, and limit violations must be zero. Safe refusal is
reported as refusal, not as successful parsing. Rust's memory model alone is not
security evidence because unsafe code, FFI, dependencies, and denial-of-service
risks remain in scope.

## T5 — Stereo safety

The development suite contains 300 unique structures. A separate gate checks
155 CIP corpus structures under 32 generated SMILES spellings plus the source
spelling, for 5,115 inputs. These are development regressions, not unused or
absolute-gold evaluation.

The remaining exit requires atom/bond-order permutations, MOL/SDF V2000/V3000
round trips, stereoisomer-key collision checks, and an unused challenge set.
Metrics must distinguish correct assignment, safe abstention, wrong confident
label, information loss, false merge, and unsupported representation.

## T6 — Maintenance and external validation

Track RDKit, COSMolKit, Open Babel, Indigo, CDK, and OpenChemLib by affected
operation, not by feature count. Security fixes require an advisory-to-release
drill. Independent accuracy claims require non-maintainer review and absolute
gold labels. 3D/MMFF94/UFF remains a separate experimental profile.

## Candidate decision

v1.0.16 may package completed trust infrastructure and interoperability work
without waiting for the full Trust Release. Its release notes must state which
gates are complete and which evaluation remains unrun. Completion of A0–A6 is
required only for broader compatibility/equivalence claims described in the
[accuracy plan](rdkit-accuracy-plan.md).
