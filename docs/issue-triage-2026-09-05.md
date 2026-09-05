# Open issue triage — 2026-09-05

This is a repository-local triage snapshot. An issue is marked **implemented**
only when the corresponding code and tests are present in the current
checkout. GitHub issue state is intentionally not changed by this document.

## Priority order

| Priority | Issue | Status in this checkout | Next action |
| --- | --- | --- | --- |
| P0 correctness | #149 shared E/Z carrier resolution | Wave 3 audit strengthened to 64 seeded RDKit relabelings per molecule and found 3/28 coupled components with two canonical outputs (0 correspondence failures); exact variants are pinned in `validation/manifests/canonical_issue149_aromatic_stash_residual.json`, and `canonical_smiles_stable_key()` rejects them | Normalize aromatic direction-stash interpretation/carrier election independently of input spelling; do not choose a winner by atom/bond index |
| P1 search performance | #139 VF2 automorphism pruning | Closed on GitHub after constrained-query ordering proved the known symmetric negative within the 1,000,000-visit budget | Retain the typed budget contract and regression fixture |
| P1 CI reliability | #70 Criterion process-level gate | Two-build null-control and checked-in synthetic calibration manifest are implemented; local contract passes | Hosted calibration remains required: real +5%, +10%, and contamination experiments |
| P1 chemical correctness | #337 MMFF94 typing residual | Isothiocyanate/CSP sub-bug is fixed and tested; the remaining 6 pyridinium/macrocycle molecules and 32 atoms are pinned in `validation/manifests/mmff94_issue337_pyridinium_sssr_residual.json`. The all-root probe found 0 candidates missing from the existing D2-root population, so D2 root enumeration is not the cause | Keep out of atom-typing fixes and root-set expansion; implement and independently validate a relevant-cycle / minimum-cycle-basis tie-break |
| P1 safety | #185, #210 UFF/3D residuals | #210 closed for the five named legacy-coordinate witnesses; broader experimental 3D/MMFF94 scope remains separate | Keep unsupported experimental cases typed and fail-closed |
| P2 performance | #372 canonical Boc/tBu symmetry | Local minimized-equivalent lane and stage counters are now recorded in `validation/results/canonical_issue372_local_2026-09-05.md`; correctness remains green | Obtain the exact RENKIN held-out witness, then compare it against the already-safe exact twin/orbit path; do not claim the preferred 2x target from the local proxy |
| P2 layout | #255, #256 | New connectivity-ordered path contains regression coverage; legacy `generate_coords` remains compatibility behavior | Do not claim legacy defects closed until the legacy API itself is changed or explicitly retained as a documented residual |
| P3 scope | #460–#463, #473 | Bounded typed APIs and downstream boundary documentation exist; full rich semantics remain intentionally unsupported | Keep the bounded contract; split any future full RXN/CDXML/polymer/biopolymer work into separate schemas and fixtures |
| P1 release hygiene | #474 | Versioned release metadata schema, v1.0.7 document, no-dependency validator, and tag-driven GitHub Release attachment are implemented and included in the v1.0.7 release path | Verify the v1.0.7 attached asset and keep registry/artifact measurements explicitly version-pinned |

## Already present and verified locally

- #299: batch fingerprint CLI and partial per-record error manifest exist in
  `chematic-cli`, with CLI tests and documentation; GitHub issue closed after
  review of the current checkout.
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
- `validation/results/mmff94_issue337_all_root_cycle_probe_2026-09-05.md` — D2-root versus all-root diagnostic: all six residual molecules had identical candidate sets and zero missing same-size candidates.
- `validation/results/mmff94_issue337_edge_exchange_probe_2026-09-05.md` — bounded relevant-cycle probe: all 6 residuals expose same-size GF(2)-exchangeable alternatives; exact-cycle populations (4 or 16) exceed the RDKit representative populations (2–4), so a permutation-invariant relevant-cycle selector is still required.
- `docs/rfcs/mmff94_relevant_cycle_selector.md` — selector boundary and acceptance gates fixed before any production candidate-family change.
- `validation/results/canonical_issue372_local_2026-09-05.md` — local multi-Boc/multi-pivaloyl proxy: 432 exhaustive leaves to 1 orbit-pruned leaf, 8 nodes, and zero old/new correctness mismatches; the exact RENKIN witness remains external.
