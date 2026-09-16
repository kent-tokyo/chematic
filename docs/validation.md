# Validation report

Updated 2026-09-16. Current public release: **v1.0.16**. Results below keep the
version, comparator, corpus, and operation written in each artifact; they are not
silently upgraded to the current source revision.

## Current trust evidence

| Area | Evidence | Boundary |
|---|---|---|
| Release channels | `validation/results/release-channel-verification-v1.0.15.json` | GitHub, crates.io, docs.rs, PyPI, npm, Pages; platform smoke uses the actually published wheels |
| RDKit.js browser comparison | `benchmarks/2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md` | Fixed exposed 10k corpus and local no-store assets; not internet/CDN latency |
| Parser security | `validation/parser_security_corpus_v1.json` plus hosted Linux gate | Fixed five-format corpus with process/time/memory boundaries |
| V3000 interchange | `validation/results/v3000-*-v1.0.15.json` | Ordinary structures and declared SGROUP/stereo contracts; no coordination/haptic/polymer semantic claim |
| Stereo development | `validation/results/stereo-*-v1.0.15-2026-09-16.json` | 300 development structures and 5,115 spelling variants; not independent gold |
| Sealed evaluation | `validation/results/sealed-cohort-preflight-trust-eval-candidate-20260916.json` | Candidate, attestation, and 2k/8k split are sealed; no score has been calculated |

## Accuracy snapshots

The exposed 4,999/5,000-molecule ChEMBL-derived lanes report exact or
tolerance-matched results for molecular weight, HBA/HBD, TPSA, LogP, molar
refractivity, Fsp3, ring families, rotatable bonds, and related descriptors.
These are regression and compatibility evidence for their recorded versions,
not results from the sealed 8,000-row holdout.

The modern CIP snapshot reports 4,171/4,186 resolved labels agreeing with the
pinned RDKit labeler; 15 phosphorus rows fail closed as representation-unstable.
Safe abstention is counted separately from a correct assignment.

The RDKit-compatible Morgan/search profile has exact exposed-corpus lanes for
configured folded bits, sparse counts, bit information, and top-k retrieval.
Native ECFP4 uses a different definition; cross-profile recall is diagnostic and
is not a compatibility percentage.

## Browser comparison summary

The published-package Chromium local no-store download-to-ready p95 is 51.62 ms
for chematic and 98.495 ms for RDKit.js in the recorded environment. chematic's
raw/gzip WASM assets are smaller and its parse/write lane is faster there;
RDKit.js is faster for the parse-inclusive fingerprint operation. One Fe(II)
coordination input is a typed fingerprint refusal, leaving 9,999 supported rows.

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
and the [benchmark index](../benchmarks/README.md) for exact boundaries.
