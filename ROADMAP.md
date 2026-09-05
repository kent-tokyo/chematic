# chematic roadmap

> Revised 2026-09-05. The current published release is v1.0.6. The workspace
> version is fixed at 1.0.6 for this release.

The detailed roadmap and completed gate-by-gate evidence through 2026-09-05 is
retained in
[`docs/archive/roadmap-through-2026-09-05.md`](docs/archive/roadmap-through-2026-09-05.md).
This file is the active plan.

## North star

Make chematic a dependable, embeddable cheminformatics runtime for Rust,
Python, JavaScript/WASM, and AI agents. Compatibility claims must be named and
measured; unsupported work must fail explicitly; untrusted input must remain
bounded.

RDKit remains the broadest and most mature reference, Open Babel the format
conversion reference, CDK a reaction/SMARTS/QSAR reference, sdfrust a Rust
dataset reference, and kekule a polymer/modeling reference. COSMolKit is not
part of the active comparison program.

## Current candidate

Completed on the v1.0.6 release tree:

- [x] Close #210's five named legacy-coordinate UFF stereo-rescue residuals.
  Every returned geometry is finite, bond-sane, and independently checked
  against declared stereo.
- [x] Exceed the additional local 1.10x hot-path target: canonical SMILES
  1.115x, file-backed SDF read 1.130x, and V2000 SDF serialization 3.042x at
  the recorded medians.
- [x] Preserve output-byte, malformed-input, workspace-test, Python-binding,
  rustfmt, and warnings-as-errors clippy gates for those changes.
- [x] Reorganize public documentation: concise README files, release-focused
  CHANGELOG, active ROADMAP, corrected publication status, and archived
  detailed development evidence.

The performance figures are source-level measurements recorded before the
release artifact was built; they remain scoped to their named corpus and
configuration.

## Priority phases

| Priority | Goal | Exit evidence | Cost |
|---|---|---|---|
| P0 | Trust and measurement | Version-pinned schema, corpus hashes, typed statuses, offline validation | Critical / light |
| P1 | Interchange throughput | Same-input SDF/MOL/XYZ measurements, memory/error gates, loss-preserving streaming | Critical / heavy |
| P2 | Identity and ML primitives | Stable canonical key, held-out descriptor/fingerprint reports, explainable output | Critical / heavy |
| P3 | Portable production surface | One Rust/Python/Node/WASM fixture contract, clean-install and browser evidence | High / heavy |
| P4 | Chemistry workflows | Reaction/SMARTS/standardization provenance and ambiguity contracts | High / heavy |
| P5 | 3D and materials | Soundness, unit/frame round trips, per-class quality and failure rates | Medium / heavy |
| P6 | Ecosystem durability | Reproducible dashboards, migration paths, contributor corpus policy | Medium / external |

Complete the highest-priority reproducible local slice before starting lower-
priority breadth. A feature is not complete until implementation, tests,
documentation, and required measurement agree.

## P0 — Trust and measurement

- [x] Separate capability claims from measured results and preserve
  `unsupported`, `failure`, and `not_measured` as distinct states.
- [x] Pin benchmark protocols, corpus hashes, environment metadata, and
  historical result records.
- [x] Maintain a capability matrix for RDKit, Open Babel, CDK, sdfrust,
  kekule, and chematic.
- [x] Exclude COSMolKit from the active comparison scope.
- [x] Add a two-build Criterion null-control: independently compile `main` and
  `main-null` and compare them through the same process-level pipeline, so
  build/codegen variance is measured before a real regression can block.
- [x] Check in the #70 synthetic calibration contract for +5%, +10%,
  sub-threshold build noise, and contamination routing; hosted real-run
  evidence remains a separate gate.
- [x] Add versioned machine-readable release metadata and a tag-driven GitHub
  Release attachment. The schema, versioned raw JSON, validator, and historical
  benchmark separation are checked in under `docs/`, `release-metadata/`, and
  `scripts/`.
- [ ] Add a scorecard validator that rejects stale release versions, missing
  corpus/configuration metadata, and claims derived from unsupported rows.

## P1 — Interchange throughput and safety

- [x] Add file-backed SDF/MOL/XYZ benchmark fixtures and a resumable runner.
- [x] Keep the RDKit block-parser comparator explicitly separate from
  chematic's file-backed streaming boundary.
- [x] Record a 2,000-pass SDF/MOL/XYZ streaming lane with zero malformed
  fixture failures and explicit cross-engine boundary notes.
- [ ] Extend the common benchmark to V2000/V3000, XYZ, MOL2, CML, CDXML,
  PDB/mmCIF, and gzip, including malformed and oversized inputs.
- [ ] Add bounded streaming batch APIs with cancellation, backpressure,
  deterministic ordering, and an explicit partial-result manifest.
- [ ] Measure only equivalent operations against installed RDKit and Open
  Babel versions on identical inputs; report sdfrust separately.

## P2 — Identity and ML primitives

- [x] Make `canonical_smiles_stable_key()` the only recommended dedup/cache
  path; it fails closed when stability is not proven, including coupled E/Z
  systems using aromatic direction stashes.
- [x] Keep native and RDKit-compatible fingerprint modes separate in API names
  and documentation.
- [x] Exclude Spectrophores and Issue #464's proposed replacement pending
  independent patent/FTO review.
- [ ] Finish canonical atom-order and E/Z invariance for the supported domain.
- [x] Prevent ring-closure close-side carrier selection from erasing an E/Z
  marker; the regression fixture and focused canonical suite are green. This
  is a bounded residual fix, not completion of the full corpus gate.
- [x] Resolve coupled E/Z carrier choices jointly across shared physical
  bonds, with 19 permutation-invariance fixtures; retain fail-closed
  abstention when no conflict-free assignment is proven.
- [x] Re-audit the committed 5,000-line carrier corpus with 64 seeded
  relabelings per coupled molecule: 28 coupled components, all size 2, and 0
  correspondence failures. The audit now reproducibly exposes 3 residual
  molecules with two canonical outputs; keep #149 open until their aromatic
  carrier/stereo traversal is normalized.
- [x] Add a bounded default-Hückel fused/non-alternant fallback for
  all-carbon odd/odd envelopes, with azulene regression coverage; keep the
  broader `RdkitLike` model separately gated.
- [x] Add a held-out CIP parity gate requiring zero wrong confident labels and
  zero regressions in the non-phosphorus scope (140/140 current cases), and
  fail closed for all 15 representation-unstable phosphorus rows.
- [ ] Freeze descriptor/fingerprint shape, sparse/count, configuration,
  provenance, and explanation contracts.
- [x] Record a local #372 canonical-symmetry proxy with per-fixture nodes,
  orbit tests, leaves, pruning, timing, and old/new correctness checks; keep
  the exact downstream RENKIN witness as an external held-out gate.
- [x] Add descriptor field provenance and a shared Rust/Python/Node/WASM core
  descriptor fixture; run the 4,999-molecule MW/TPSA/HBD/HBA/heavy-atom lane.
- [ ] Add held-out parity reports for Morgan/ECFP, MACCS, topological,
  torsion, descriptors, and standardization across Rust/Python/WASM.

- [x] Make VF2 query expansion prefer mapped-neighbor and higher-degree atoms,
  preserving exhaustive semantics while resolving the known symmetric PAINS
  negative within the production visit budget; keep the typed budget outcome.

## P3 — Portable production surface

- [x] Provide typed, bounded parser and binding failures with shared
  adversarial fixtures.
- [x] Run Chromium, Firefox, and WebKit smoke/adversarial lanes for the
  published v1.0 boundary.
- [ ] Make Rust, Python, Node, and WASM consume one fixture schema and one
  versioned expected-result manifest for every shared stable operation.
- [ ] Publish current clean-install, cold-start, throughput, peak-memory, and
  WASM-size evidence with explicit platform/configuration metadata.
- [ ] Extend browser and agent adversarial cases for cancellation, malformed
  records, limits, and stable JSON errors.

## P4 — Chemistry workflow depth

- [x] Add typed reaction documents, RXN adapters, loss-preserving CDXML
  commands, and bounded semantic Markush/polymer expansion.
- [x] Preserve agents, coefficients, conditions, atom maps, source mappings,
  and unsupported-richness errors instead of flattening them.
- [x] Publish the issue #473 downstream capability matrix, including the
  explicit first-class nucleic-acid/biopolymer non-support boundary and the
  Rust/Python/WASM/Node entry points.
- [ ] Expand reaction/SMARTS/medicinal-chemistry coverage only after P0-P3
  gates have current evidence.
- [ ] Add curated reaction/query precision, recall, invalid-product, timeout,
  and ambiguity reports.

## P5 — 3D and materials

- [x] Keep 3D generation and MMFF94 Experimental, with typed failure and
  explicit force-field/fallback provenance.
- [x] Separate long-running 3D and corpus-scale canonical tests into explicit
  ignored lanes and retain their execution manifest.
- [x] Add crystal composition and materials-format foundations without
  conflating periodic structures with molecular bond graphs.
- [ ] Close MMFF94/UFF typing, charge, parameter, convergence, and stereo gaps
  with independent oracle and soundness gates. The six-molecule/32-atom
  pyridinium/macrocycle residual boundary is pinned in
  `validation/manifests/mmff94_issue337_pyridinium_sssr_residual.json`; do not
  replace this with a local atom-type heuristic.
- [x] Measure the #337 D2-root hypothesis with an all-root diagnostic; the
  six residual molecules have identical D2/all-root candidate sets, so root
  enumeration is not treated as the fix.
- [x] Probe equal-length macrocycle edge exchanges for #337; all 6 residuals
  yield same-size GF(2)-exchangeable alternatives. The candidate population
  is intentionally not promoted until a permutation-invariant representative
  tie-break is validated.
- [x] Specify the bounded #337 relevant-cycle selector boundary and fail-closed
  acceptance gates in `docs/rfcs/mmff94_relevant_cycle_selector.md`.
- [ ] Measure deterministic ensemble diversity, class-level failure rates,
  symmetry-aware RMSD/TFD, and energy sanity.

## P6 — Ecosystem durability

- [x] Publish compatibility, provenance, licensing, security, benchmark, and
  migration documents that state unsupported scope explicitly.
- [x] Record the open-issue triage and distinguish implemented bounded APIs
  from unresolved correctness, safety, performance, and CI work.
- [ ] Add stable extension points and a contributor corpus policy covering
  provenance, licensing, minimization, and oracle versioning.
- [ ] Publish a reproducible compatibility dashboard only after a clean
  checkout can regenerate it.

## Security and release gates

- [x] S0/S1: threat model, bounded public input, typed resource-limit errors,
  and parser-wide adversarial coverage.
- [x] S2: unsafe-surface inventory, optional native-InChI FFI isolation, four
  fuzz targets, focused Miri, and ASan/LSan/TSan procedures and evidence.
- [x] S3: shared binding contracts plus Python, Node/WASM, browser, MCP, CLI,
  and filesystem boundary tests.
- [x] S4: dependency/license checks, immutable Actions, SBOM/provenance,
  checksums, release-key custody, and registry verification.
- [x] S5 local gate: the repository-local review packet and executable review
  checks are complete.
- [ ] S5 independent gate: obtain a non-maintainer review or external audit of
  major parser, serialization, and binding boundaries. This is external work
  and must not be marked complete from self-review.
- [ ] S6 continuous maintenance: rehearse advisory intake, fix, backport,
  artifact publication, and supported-version synchronization.

Reproduction instructions are in
[`docs/v1.0-local-release-gate.md`](docs/v1.0-local-release-gate.md),
[`docs/security-review/`](docs/security-review/), and
[`SECURITY.md`](SECURITY.md).

## Compatibility boundaries

The v1.0 contract remains unchanged:

- `canonical_smiles()` is a representation; the stable-key API is fail-closed.
- Aromaticity and CIP expose explicit default/opt-in models and do not promise
  universal RDKit parity.
- Python `RWMol`, CDXML editing, and Markush/polymer expansion are bounded
  subsets.
- 3D/MMFF94 are Experimental.
- Pure-Rust InChI is approximate; exact standard InChI is opt-in native FFI.

See [`docs/compatibility-scope.md`](docs/compatibility-scope.md).

## Execution order

1. Finish the P0 scorecard validator.
2. Extend equivalent format-streaming benchmarks and input-safety fixtures.
3. Close canonical E/Z invariance before broadening identity-dependent APIs.
4. Unify stable cross-binding fixtures and expected results.
5. Stabilize fingerprint/descriptor contracts and held-out reports.
6. Advance reaction, 3D/materials, and ecosystem breadth only behind their
   respective correctness and measurement gates.

Do not start broad feature expansion while a shared primitive has a known
silent-corruption or unbounded-work regression.

## Definition of “ahead”

chematic is ahead of a competitor only for a named use case and a measured
dimension. Every claim must identify versions, corpus, configuration, hardware,
failure policy, and reproduction command. “Faster”, “more accurate”, “smaller”,
and “more compatible” are otherwise hypotheses, not product claims.
