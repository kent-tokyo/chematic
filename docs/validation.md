# Validation report

Updated 2026-09-25. The current release is **v1.0.25**. Each result keeps its
recorded version, comparator, corpus, and operation; no result is silently
upgraded to the current source revision.

## Current trust evidence

| Area | Evidence | Boundary |
|---|---|---|
| v1.0.25 release-source identity | `benchmarks/2026-09-25-perf-digest-diff-v1.0.24-issue-650.json` | Six pre-existing-output operations across 46,736 molecules; zero differences. The HBA profile and plain writer corrections are intentional behavior changes and are outside this unchanged-output claim. |
| Source performance identity | `benchmarks/2026-09-24-perf-speed3-output-identical.md` | Differential against v1.0.23 over 1,402,080 output rows; zero differences. Shared-VM source timing only, not a package or cross-platform result. |
| Release channels | `validation/results/release-channel-verification-v1.0.24.json` | Historical v1.0.24 channel observation. v1.0.25 requires a separate post-publication record; it is not implied by this source release document. |
| RDKit output agreement | `benchmarks/2026-09-24-rdkit-agreement-accuracy-branch.md` | RDKit 2026.03.6; three source lanes of up to 5,000 rows. Operation-specific agreement and residuals only; shared 2-vCPU VM, not a package/WASM/cross-platform result. |
| Source operation matrix | `benchmarks/2026-09-24-python-op-matrix-vs-rdkit-perf-branch.md` | RDKit 2026.03.6, 5,000 molecules, 43 comparable operations; 21 exact-output faster medians and 18 all-repeat wins. Shared 2-vCPU source run only; not a package, WASM, or cross-platform claim. |
| RDKit.js browser comparison | `benchmarks/2026-09-23-public-package-fingerprint-3d-v1.0.20.md` | Registry-installed v1.0.20 on the fixed exposed 10k corpus; both compatible-Morgan lanes pass and all 9,999 supported rows are bit-exact; not internet/CDN latency |
| MMFF94 same-coordinate energy | `benchmarks/2026-09-23-mmff94-current-source-energy.md` | Current-source 265-row packet against pinned RDKit 2026.03.6: 262 comparable, p90 absolute delta 1.144679 kcal/mol, two residuals above 5 kcal/mol; not conformer quality, convergence, stereo, speed, or a published-package claim |
| A6 source speed/convergence candidate | `benchmarks/2026-09-23-a6-analytic-mmff94-source-candidate.md` | Historical source packet; use the v1.0.20 public-package record for current package claims |
| A6 MMFF94 stereo-safe quality | `benchmarks/2026-09-23-public-package-fingerprint-3d-v1.0.20.md` | Public v1.0.20 closes the prior 12 typed failures and four gross-clash rows at 265/265 independently sound, stereo-clean, clash-free successes. Paired speed is 0.944x with a 0.861x lower bound, so this is not an MMFF94 speed win. |
| RDKit 2026.03.6 rebaseline | `validation/results/rdkit-rebaseline-*-v1.0.19-vs-2026.03.6-2026-09-22.*` plus the 2026-09-23 native-availability record | Exact Python/npm artifacts and 10,000 complete exposed rows; correspondence-correct CIP is 9,770/10,000 (the older 9,880 engine-local-index result is withdrawn), with 18 classified SMILES stereo-writer regressions, one typed Morgan contract difference, and unresolved CIP/SMARTS residuals. The independent native/C++ lane remains explicitly unavailable; historical diagnostics, not oracle adoption |
| Issue #632 SMILES release-source diagnostic | `validation/results/smiles-ez-semantic-issue632-v1.0.20-candidate-vs-rdkit-2026.03.6-2026-09-23.json` | v1.0.21 contains source `9808f54f`, which reduces the same exposed 10k lane from 18 semantic differences to zero with zero graph differences; CIP, Morgan, and SMARTS counts are unchanged. The long 28 x 1,024 relabel audit was interrupted, so this is not a completed long-audit or published-package remeasurement claim. |
| Parser security | `validation/parser_security_corpus_v1.json` plus hosted Linux gate | Fixed five-format corpus with process/time/memory boundaries |
| V3000 interchange | `validation/results/v3000-*-v1.0.15.json` | Ordinary structures and declared SGROUP/stereo contracts; no coordination/haptic/polymer semantic claim |
| Stereo development | `validation/results/stereo-*-v1.0.15-2026-09-16.json` | 300 development structures and 5,115 spelling variants; not independent gold |
| Sealed evaluation | `validation/results/a0-core-eight-sealed-acceptance-20260922.json` plus historical rejection/TPSA summaries | Candidate `5e9211a6` passed all eight declared fields on a separately sourced one-time 8k holdout after post-freeze acquisition and overlap audit. All raw sealed rows remain local-only and exposed/ineligible for later candidates. |

## Accuracy snapshots

The exposed 4,999/5,000-molecule ChEMBL-derived lanes report exact or
tolerance-matched results for molecular weight, HBA/HBD, TPSA, LogP, molar
refractivity, Fsp3, ring families, rotatable bonds, and related descriptors.
These are regression and compatibility evidence for their recorded versions,
not results from the sealed 8,000-row holdout.

For v1.0.19, the replacement candidate `bac7ae44` was frozen as
`trust-eval-candidate-20260921` before a new 14,764-row ChEMBL source was
acquired. Canonical/parent/scaffold overlap was audited against 20,000 exposed
rows, leaving 12,345 eligible rows. TPSA then passed 2,000/2,000 development
rows and the one-time sealed holdout at 8,000/8,000 with `1e-6` tolerance and
maximum absolute error `0.0`. Only TPSA was measured in that sealed run; the
run remains a historical TPSA-only result.

The post-v1.0.19 candidate `5e9211a6` was frozen as
`trust-eval-candidate-a0-multifield-r3-20260922` before a third independent
source acquisition. After canonical/parent/scaffold overlap exclusion against
all exposed sources, the eight-field profile passed 2,000/2,000 development
rows and 8,000/8,000 sealed rows for every field, with parse failures,
unsupported values, and mismatches all zero. LogP and molar-refractivity
floating-point deltas remained below `5e-14` and `2.1e-12`, respectively, under
the declared `1e-6` tolerance. Run `python3 scripts/check_a0_core_eight_sealed.py`
to validate the commit-safe aggregate. This completes A0 core-eight adoption;
additional descriptor families, potential centers, stereo, retrieval,
interchange, independent gold, and 3D retain their separate A1–A6 gates.

The modern CIP snapshot reports 4,171/4,186 resolved labels agreeing with the
pinned RDKit labeler; 15 phosphorus rows fail closed as representation-unstable.
Safe abstention is counted separately from a correct assignment.

The RDKit-compatible Morgan/search profile has exact exposed-corpus lanes for
configured folded bits, sparse counts, bit information, and top-k retrieval.
Native ECFP4 uses a different definition; cross-profile recall is diagnostic and
is not a compatibility percentage.

## Browser comparison summary

The 2026-09-23 registry-installed v1.0.20 Chromium record passes both compatible-
Morgan speed gates: 1.398x parse-inclusive and 3.511x prepared, with 95% lower
bounds 1.363x and 3.407x. Both paths preserve 9,999/9,999 configured-bit
agreement; one Fe(II) coordination input remains a typed refusal.

These numbers do not establish internet download latency, unique process memory,
all-browser performance superiority, or unmeasured fingerprint configurations.

## Reproduction entry points

```bash
# Core workspace checks
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings

# Documentation and evidence consistency
python3 scripts/check_release_docs_consistency.py
python3 scripts/check_benchmark_index.py
python3 scripts/check_compatibility_profiles.py
python3 scripts/check_rdkit_rebaseline_execution.py
python3 scripts/check_rdkit_rebaseline_evidence.py

# Development accuracy snapshot (requires the pinned RDKit environment)
python3 scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi
```

Use the exact command and environment recorded by an artifact for formal
reproduction. The commands above are entry points, not substitutes for its
pinned metadata.

## Known limits

- Canonical SMILES is not always a safe identity key; use
  `canonical_smiles_stable_key()` where fail-closed behavior is required.
- Coupled aromatic E/Z and phosphorus-CIP cases retain explicit residuals.
- Coordination/haptic V3000 semantics, broad polymer expansion, and full CDXML
  editing are outside the stable contract.
- Pure-Rust InChI is approximate; standard InChI requires the optional native
  feature.
- 3D/MMFF94/UFF remains experimental until the separate A6 gates pass.

See [compatibility scope](compatibility-scope.md),
[RDKit migration](rdkit-migration.md), [accuracy plan](rdkit-accuracy-plan.md),
and the [benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)
for exact boundaries.
