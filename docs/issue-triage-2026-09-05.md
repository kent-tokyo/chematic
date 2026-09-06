# Open issue triage — 2026-09-05

This is a repository-local triage snapshot. An issue is marked **implemented**
only when the corresponding code and tests are present in the current
checkout. GitHub issue state is intentionally not changed by this document.

## Priority order

| Priority | Issue | Status in this checkout | Next action |
| --- | --- | --- | --- |
| P0 correctness | #149 shared E/Z carrier resolution | Wave 3 audit remains 3/28 divergent after 64 seeded RDKit relabelings per molecule (0 correspondence failures). Canonical-search edge coloring preserves writer-visible aromatic direction stashes during automorphism checks; the related coloring documentation is now synchronized with that contract. All 213 crate tests and 6 residual tests pass, but the divergence count is unchanged, so the residual remains open and `canonical_smiles_stable_key()` still fails closed | Continue with a component-level canonical winner proof for aromatic stash carriers; do not choose a winner by atom/bond index or claim resolution from the edge-color hardening |
| P1 search performance | #139 VF2 automorphism pruning | Closed on GitHub after constrained-query ordering proved the known symmetric negative within the 1,000,000-visit budget | Retain the typed budget contract and regression fixture |
| P1 CI reliability | #70 Criterion process-level gate | Process-level observations, ABBA/BAAB ordering, two-build null-control, strict ratio gates, and per-block artifact metadata (timestamps, execution order, load average, CPU model, and steal ticks) are implemented; local contract passes | Hosted calibration remains required: real +5%, +10%, and contamination experiments |
| P1 chemical correctness | #337 MMFF94 typing residual | Isothiocyanate/CSP sub-bug is fixed and tested; the remaining 6 pyridinium/macrocycle molecules and 32 atoms are pinned in `validation/manifests/mmff94_issue337_pyridinium_sssr_residual.json`. The all-root probe found 0 candidates missing from the existing D2-root population, so D2 root enumeration is not the cause. The symmetrized-ring path now uses a permutation-invariant candidate-root set, basis-exchange ordering no longer uses raw bond indices, and expansion now fails closed to the complete Horton basis at a 256-extra-ring cap rather than returning a partial family. Full symmetrized-ring ordering is regression-tested under relabeling, but a fresh six-fixture comparison still shows the same aromaticity/type residuals | Independently validate the relevant-cycle / minimum-cycle-basis representative family against the six-molecule oracle; do not claim #337 resolved from the safety cap or tie-break determinism alone |
| P1 safety | #185, #210 UFF/3D residuals | #210 closed for the five named legacy-coordinate witnesses. #185's bounded slice is implemented: UFF rejects unsound line-search proposals and exposes additive `sound` state through the Rust result, Python binding, and WASM JSON; full UFF torsion/OOP terms remain open | Keep unsupported experimental cases typed and fail-closed; do not equate `converged` with geometrical soundness |
| P2 performance | #372 canonical Boc/tBu symmetry | Local minimized-equivalent lane and stage counters are now recorded in `validation/results/canonical_issue372_local_2026-09-05.md`; correctness remains green | Obtain the exact RENKIN held-out witness, then compare it against the already-safe exact twin/orbit path; do not claim the preferred 2x target from the local proxy |
| P2 layout | #246, #255, #256 | #246 resolved in this checkout: bridged-ring anchoring now scores both regular-polygon sides against all already-placed ring atoms, and uses a deterministic closest-pair fallback when no shared edge exists. #255/#256 remain resolved via the connectivity-ordered 3D engine; fresh 33-molecule evaluation remains raw sound 33/33, deterministic 33/33, and UFF-only success 33/33 | Keep the bridged-ring bond-length regression and the 3D differential harness before future placement changes |
| P3 scope | #460–#463, #473 | Bounded typed APIs and downstream boundary documentation exist; full rich semantics remain intentionally unsupported | Keep the bounded contract; split any future full RXN/CDXML/polymer/biopolymer work into separate schemas and fixtures |
| P1 release hygiene | #474 | Versioned release metadata schema, v1.0.7 document, no-dependency validator, and tag-driven GitHub Release attachment are implemented and included in the v1.0.7 release path | Verify the v1.0.7 attached asset and keep registry/artifact measurements explicitly version-pinned |
| P3 crystal identity | #477 | Resolved in this checkout: added versioned `PeriodicStructure::identity_bytes()` plus a pure-Rust SHA-256 `identity_digest()` for deterministic exact-identity cache/provenance keys; crystal tests pass | Keep the version byte in the hashed identity bytes; the digest is not symmetry canonicalization or a material-similarity score |
| P1 ingestion | #478 | Resolved in this checkout: `chematic-smiles::SmilesBatchCanonicalizer` provides lazy iterator and newline-delimited `BufRead` results with reusable parser limits and per-record accepted/rejected diagnostics; `build_identity_index()` uses only `canonical_smiles_stable_key()`, preserves duplicate positions, and fails closed for unstable identities. Shared Rust/Python/Node/WASM fixtures and versioned partial-result envelopes now cover the JSON wrappers | Add optional parallel execution only with equivalent ordering/error fixtures |

## Already present and verified locally

- #299: batch fingerprint CLI and partial per-record error manifest exist in
  `chematic-cli`, with CLI tests and documentation; the common versioned
  envelope now covers every batch operation while preserving input order and
  effective limits. GitHub issue closed after review of the current checkout.
- The Python SDF batch iterator now rejects zero/oversized batch sizes,
  supports explicit cancellation, and exposes deterministic progress JSON;
  WASM and cross-language streaming parity remain open rather than inferred.
- The SDF/MOL/XYZ benchmark runner now accepts explicit resource limits and
  records them in JSON; a bounded input-limit run produced 0 records and 1
  failure as expected. V3000/MOL2/CML/CDXML/mmCIF parser rows are now covered
  as `materialized_one_shot`; PDB/gzip and true streaming parity remain open.
- #463: occupancy-aware `PeriodicStructure::composition()` and its disorder,
  zero-occupancy, deterministic-order, and explicit-supercell tests exist;
  GitHub issue closed after review of the current checkout.
- #460–#462: typed reaction documents, loss-aware RXN adapters,
  loss-preserving CDXML commands, and bounded semantic expansion exist.
- #473: the downstream capability matrix is documented in
  `docs/issue-473-rich-document-boundary.md`.

These entries are not treated as new implementation work merely because the
repository retains historical issue references.

## Local evidence run for this triage

- `bash scripts/test_criterion_gate.sh` — passed.
- The Criterion block artifact now persists execution order, UTC start/end timestamps,
  load average, CPU model, and `/proc/stat` steal ticks when available; the local
  contract test asserts these fields without changing baseline/candidate values.
  Each record is schema version 2 and explicitly identifies one Criterion process
  point estimate as its measurement unit.
- `cargo test -p chematic-rxn --offline` — 193 passed.
- `cargo test -p chematic-mol --lib --offline` — 538 passed.
- `cargo test -p chematic-smiles --lib ez_shared_carrier --offline` — passed.
- `cargo test -p chematic-smarts --lib --offline` — 184 passed.
- `cargo test -p chematic-chem --lib --offline` — 840 passed, 1 ignored.
- `cargo test -p chematic-chem --lib test_pains_di_tert_butylphenol_resolves_negative_within_budget --offline -- --nocapture` — passed in 1.66s; `NotFound` within the production budget.
- `TMPDIR=/private/tmp bash scripts/test_criterion_gate.sh` — passed, including the two-build null-control workflow contract checks.
- `validation/criterion-gate-calibration.json` — checked-in synthetic +5%/+10%/noise/contamination expectations; this is a contract fixture, not hosted execution evidence.
- `PYTHONDONTWRITEBYTECODE=1 python3 scripts/ez_shared_carrier_coupling_mechanism_diagnosis.py --relabelings 64 scripts/descriptor_census_corpus.smi` — current Wave 3 audit: 522 component rows, 28 coupled components, all size 2, 550 ends, 0 correspondence failures, and 3 coupled molecules with two canonical outputs.
- `cargo test -p chematic-smiles --test canonical_ez_residual --offline` — 5 passed; the three-family aromatic-stash residual is now rejected by `canonical_smiles_stable_key()` rather than accepted as a dedup/cache key.
- `cargo test -p chematic-smiles --test canonical_ez_residual --offline` — now includes the three-family exact-output/fail-closed residual contract.
- `validation/manifests/mmff94_issue337_pyridinium_sssr_residual.json` — checked-in six-molecule/32-atom residual contract; this is a local fixture, not an RDKit re-run or a claim that #337 is resolved.
- `validation/results/mmff94_issue337_all_root_cycle_probe_2026-09-05.md` — D2-root versus all-root diagnostic: all six residual molecules had identical candidate sets and zero missing same-size candidates; a fresh 2026-09-06 run reproduced the same boundary (`symm_count` 8/10/10/10/10/9 and no missing same-size candidates).
- `validation/results/mmff94_issue337_edge_exchange_probe_2026-09-05.md` — bounded relevant-cycle probe: all 6 residuals expose same-size GF(2)-exchangeable alternatives; exact-cycle populations (4 or 16) exceed the RDKit representative populations (2–4), so a permutation-invariant relevant-cycle selector is still required.
- `cargo test -p chematic-perception --lib --offline` and `cargo test -p chematic-ff --lib --offline` — passed after adding the bounded symmetrized-ring fallback; these are local regression results, not #337 resolution evidence.
- `docs/rfcs/mmff94_relevant_cycle_selector.md` — selector boundary and acceptance gates fixed before any production candidate-family change.
- `validation/results/canonical_issue372_local_2026-09-05.md` — local multi-Boc/multi-pivaloyl proxy: 432 exhaustive leaves to 1 orbit-pruned leaf, 8 nodes, and zero old/new correctness mismatches; the exact RENKIN witness remains external.
