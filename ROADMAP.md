# chematic roadmap

> Updated 2026-09-16. Current release: **v1.0.15**. The next candidate is
> **v1.0.16**, a Trust/Interoperability maintenance release. This document is
> the only source for priority order; detailed acceptance criteria live in the
> linked plans.

## Current position

v1.0.15 is published and verified across GitHub Releases, crates.io, docs.rs,
PyPI, npm, and GitHub Pages. The release proof also includes binary-only Python
smoke tests on Linux, macOS, and Windows. Post-release `main` adds:

- a published-package browser comparison against `@rdkit/rdkit@2026.03.6`;
- version-pinned ordinary-V3000 round trips through RDKit and Indigo;
- a 300-structure development Stereo Torture Suite and 5,115 SMILES-spelling
  invariance checks;
- an annotated candidate freeze and post-freeze unused-data attestation for a
  2,000-row development / 8,000-row sealed split;
- dependency and supply-chain updates, including the `rustls` advisory fix.

The sealed split is ready for evaluation, but **no sealed-holdout score has been
calculated**. Compatibility remains operation-, version-, and corpus-scoped.

## Priority order

| Order | Priority | Next result | Exit condition |
|---|---|---|---|
| 1 | P0 | Run the frozen candidate against the sealed 8,000-row holdout | Immutable raw output, complete denominators, hashes, and no post-hoc tuning |
| 2 | P0 | Publish the operation-level compatibility dashboard | RDKit version, profile, corpus, coverage, refusals, tolerance, and missing dimensions shown together |
| 3 | P1 | Extend stereo and parser-security challenge gates | Wrong confident labels, crashes, panics, and limit violations remain zero on the declared domain |
| 4 | P1 | Expand V3000 only where external readers prove semantics | Typed support and opaque retention are reported separately; coordination/ENDPTS boundaries stay explicit |
| 5 | P2 | Improve browser and agent adoption | Published-package Worker example, cancellation, partial-failure accounting, and reproducible 10k workflow |
| 6 | P3 | Continue 3D/MMFF94/UFF work as a separate profile | Same-coordinate numerical gates and conformer-quality evidence; no broad parity claim before A6 |

Confirmed silent corruption, panic, or resource-limit defects take precedence
over this order.

## Completed foundations

| Area | Current evidence |
|---|---|
| Release truth | `validation/results/release-channel-verification-v1.0.15.json` |
| Public RDKit.js cost comparison | `benchmarks/2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md` |
| Sealed-evaluation precondition | `validation/results/sealed-cohort-preflight-trust-eval-candidate-20260916.json` |
| V3000 ordinary interchange | `validation/results/v3000-*-v1.0.15.json` |
| Stereo development gates | `validation/results/stereo-*-v1.0.15-2026-09-16.json` |
| Parser-security gate | Linux isolated CI job and `validation/parser_security_corpus_v1.json` |

These records prove only their declared operation and environment. In
particular, the local no-store browser startup measurement is not an internet
or CDN result, and parse-inclusive fingerprint timing currently favors RDKit.js.

## Accuracy packages

All packages remain open until their full exit criteria pass. Narrower completed
gates are retained as evidence, not promoted to package completion.

- [ ] **A0 — Evaluation contract:** execute the frozen sealed split and publish
  complete provenance, denominators, failures, and paired statistics.
- [ ] **A1 — Perception and descriptors:** pass the sealed core descriptor and
  potential-stereocenter lanes without native-profile regression.
- [ ] **A2 — Stereo and identity:** resolve the remaining canonical E/Z and
  phosphorus-CIP boundaries; preserve permutation and round-trip invariance.
- [ ] **A3 — Fingerprints and retrieval:** preserve exact RDKit-compatible
  fingerprint/search results on the declared profile and sealed inputs.
- [ ] **A4 — Workflows and interchange:** finish SMARTS, standardization,
  reaction-product, and typed V3000 semantic gates.
- [ ] **A5 — Independent adjudication:** obtain non-maintainer review and
  absolute gold labels before claiming equivalence or superiority.
- [ ] **A6 — 3D and force fields:** resolve MMFF94/UFF typing, energy, gradient,
  convergence, timeout, stereo, and conformer-quality gaps separately.

The full protocol and numerical exits are in
[`docs/rdkit-accuracy-plan.md`](docs/rdkit-accuracy-plan.md).

## Product phases

Priority and phase are different: priority is execution order; P0–P6 are stable
product areas.

| Phase | Product area | Active focus |
|---|---|---|
| P0 | Trust and measurement | sealed evaluation, compatibility contract, release proof |
| P1 | Interchange and parser safety | V3000, streaming, malformed-input limits |
| P2 | Identity, descriptors, fingerprints, search | A1–A3 |
| P3 | Rust/Python/Node/WASM and agents | install, Worker, MCP, cancellation |
| P4 | Chemistry workflows | SMARTS, reactions, stereo, standardization |
| P5 | 3D and materials | A6 |
| P6 | Maintenance and external validation | advisories, independent review, competitor watch |

## v1.0.16 candidate boundary

v1.0.16 may ship as a maintenance release without claiming completion of A0–A6.
It requires:

1. synchronized versions and concise release notes;
2. clean formatting, tests, Clippy, Security Audit, and benchmark-index checks;
3. package smoke tests for affected Rust, Python, npm/WASM, and MCP surfaces;
4. explicit documentation that the sealed cohort is prepared but not scored;
5. post-publication verification of every release channel.

The release theme is strengthened trust evidence and interoperability—not full
RDKit replacement, independent chemical superiority, or complete V3000/3D
coverage.

## Reference documents

- [Trust Release execution plan](docs/trust-release-plan.md): T0–T6 status and
  next gates.
- [RDKit accuracy plan](docs/rdkit-accuracy-plan.md): comparator protocol and
  A0–A6 exits.
- [Open-work disposition](docs/roadmap-open-work.md): dependency classes and
  issue boundaries.
- [Validation report](docs/validation.md): current evidence summary.
- [Benchmark index](benchmarks/README.md): dated, version-scoped measurements.
- [Archived roadmap](docs/archive/roadmap-through-2026-09-13.md): detailed
  historical status before this consolidation.

Historical measurements retain their original source revision, package version,
corpus, hardware, and operation boundary. A newer release does not silently
upgrade an older result.
