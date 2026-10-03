# Validation report

Updated 2026-10-03 for **v1.0.31**. All six release channels are verified.
The broad chemistry and speed packet remains pinned to published v1.0.30;
the separately named A6 quality rows use published v1.0.31 artifacts.
The v1.0.31 WASM formula and reaction-JSON fixes have source
regression tests; published npm output parity has not been rerun. Every
comparison is limited to its recorded artifact, comparator, corpus and operation.

## Current evidence

| Gate | Result and limit | Record |
|---|---|---|
| Published-artifact output audit (P0.1) | Python/npm/Rust reran 10,000 chemistry rows, 310,000 SMARTS cells and 57 legacy reactions. Published Python/Rust cover 63 operation outputs; only HBA and its bundle changed from v1.0.29. The audit is complete, **not** strict RDKit parity: 200 SMARTS cells fail, five CIP rows abstain and four npm operations are unexposed. | [Artifact packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md) · [adoption policy](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md) |
| RDKit-compatible HBA | The hash-verified published v1.0.30 macOS arm64 wheel matches RDKit 2026.03.6 on 5,000/5,000 exposed ChEMBL rows; published v1.0.29 matched 3,641/5,000. HBA only. | [Release-channel evidence](https://github.com/kent-tokyo/chematic/blob/main/validation/results/release-channel-verification-v1.0.30.json) |
| Paired speed (P0.2) | Published Python: 20 operations pass exact-output and paired-interval gates in 20 alternating blocks on one host. Isolated Python and Rust lanes split parse-inclusive, first-use and precomputed calls. Whole-process RSS is not library allocation; changed HBA output is not an equivalent-output speed win. | [Python matrix](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) · [isolated Python](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1.0.30-isolated-python-time-memory.md) · [Rust](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1029-v1030-published-rust-isolated-time-memory.md) |
| Published browser Morgan | On Ubuntu 24.04, Chromium, Firefox and WebKit each have 250/250 direct and prepared bit parity and favorable 20-block speed intervals versus official RDKit.js (smallest lower bound 1.63×). A separate Chromium 10k lane has 9,999 exact and one typed Fe refusal. No universal browser or memory claim. | [Three-browser record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md) · [M4 Chromium](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.md) |
| Reactions (P1) | Legacy 57/57 fixtures match. Published v1.0.30 Rust graph/origin/template-map identity matches RDKit on 73/83 stratified rows. The unpublished checked-source profile has 76 exact on all three axes, three typed unsupported, one diagnosed refusal and three jointly invalid. A local Python source extension reproduces that classification; release-profile CI and published-package provenance are pending. No general SMIRKS parity claim. | [Published packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md) · [Python source gate](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-reaction-83-python-provenance-source.md) |
| Release channels | GitHub Release, npm, PyPI, the crates.io publish workflow, docs.rs and Pages are verified for v1.0.31. Direct crates.io API access returned 403 on this host; docs.rs rendered the version-pinned crate. Channel verification does not establish chemical parity. | [v1.0.31 channel record](https://github.com/kent-tokyo/chematic/blob/main/validation/results/release-channel-verification-v1.0.31.json) |

## Other bounded results

- **A0 core-eight:** frozen candidate `5e9211a6` passed 2,000/2,000
  development and a one-time 8,000/8,000 sealed holdout on eight declared
  descriptors after source and overlap checks. The raw sealed rows remain
  local-only and cannot be reused for tuning. [Commit-safe summary](https://github.com/kent-tokyo/chematic/blob/main/validation/results/a0-core-eight-sealed-acceptance-20260922.json).
- **CIP:** pinned exposed 10k comparison has 9,995 exact and five typed
  abstentions (four `oracle_unstable`, one `lone_pair_center`). This is
  neither a guessed label nor full CIP parity. [Source adjudication](https://github.com/kent-tokyo/chematic/blob/main/validation/results/rdkit-rebaseline-issue634-v1.0.27-candidate-vs-rdkit-2026.03.6-2026-09-28.json).
- **Morgan:** the published M4 Chromium lane matches on 9,999/9,999
  supported 10k inputs; one Fe(II) coordination input is a typed refusal.
  Native ECFP4 is a different profile. [Row-level record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-chromium-morgan-row-parity.md).
- **SMARTS source candidate:** An opt-in Python source-wheel gate matches
  309,982/310,000 RDKit 2026.03.6 match sets; the other 18 are typed
  unsupported. This does not change the published v1.0.30 baseline of 200
  differing sets. [Source packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-python-source-smarts-optin-310k.md).
- **3D/MMFF94:** Experimental. Published v1.0.31 macOS quality checks pass
  geometry/stereo/clash on 265/265, but only 100/265 converge at the declared
  limit; 262 comparable same-coordinate total energies are within 1 kcal/mol.
  Linux/Python 3.9 published and source wheels both have typed stereo failures
  on rows 53/246, so cross-platform equivalence and independent conformer
  quality remain open. [Quality record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-a6-published-v1031-mmff94-quality.md)
  · [platform diagnosis](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-a6-source-mmff94-termination.md).

Earlier source diagnostics and per-release channel checks remain in the
[benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md) and
[versioned validation results](https://github.com/kent-tokyo/chematic/tree/main/validation/results/). They are not promoted
to v1.0.31 measurements.

## Reproduce and interpret

```bash
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
python3 scripts/check_release_docs_consistency.py
python3 scripts/check_benchmark_index.py
python3 scripts/check_v1030_artifact_packet.py
```

Formal reruns must use each record's exact versions, hashes, corpus and
runtime. Canonical SMILES is not a universal identity key; use the bounded
`canonical_smiles_stable_key()` when its documented domain applies. Rich
CDXML, coordination/haptic V3000, broad Markush/polymer expansion and
Standard InChI outside the optional native feature remain bounded. See the
[compatibility scope](compatibility-scope.md) and
[accuracy plan](rdkit-accuracy-plan.md).
