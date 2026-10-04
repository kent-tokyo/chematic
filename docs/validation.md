# Validation report

Updated 2026-10-04 for **v1.0.34**. The broad chemistry
and speed packet remains pinned to published v1.0.30; the separately named A6
quality rows use published v1.0.31 artifacts. v1.0.34 is verified on six
publication channels, while its package-output reruns remain open. Every
comparison is limited to its recorded artifact, comparator, corpus and operation.

## Current evidence

| Gate | Result and limit | Record |
|---|---|---|
| RDKit 2026.09.1 npm rebaseline | Published RDKit.js 2026.03.6 → 2026.09.1 on exposed 10k: both 9,999/9,999 exact CheMatic 1.0.33 Morgan bits with one typed Fe refusal. Graph-checked CIP: 9,988 exact, five typed abstentions, six pre-existing imine E/Z mismatches and one unproven row in both versions. Published CheMatic WASM SMARTS atom sets: 194 → 195 mismatches among 310,000 cells, with 31 index-unproven cells in each lane; the two RDKit versions differ on 12 `[R2]`/`[R3]` cells. No old/new parse/canonical/Morgan changes. Isolated 20-run browser timing and hashes recorded. Python/nanobind/native 2026.09.1 lanes are unmeasured pending a distributed artifact. | [npm rebaseline](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-rdkit-2026-09-1-npm-rebaseline.md) |
| RDKit Python 2026.03.6 baseline for rebaseline | Published macOS arm64 CPython 3.13 wheels, same exposed 10k: 9,989 CIP exact, five typed abstentions and six imine E/Z mismatches after correcting the cis/trans comparator; 9,999 Morgan exact plus one typed Fe refusal; 200/310,000 SMARTS atom-set differences; all semantic SMILES round-trips pass. Runtime is Boost.Python; seven old-version-only boundary timings recorded. This is **not** a 2026.09.1 Python result. | [Python baseline](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-rdkit-2026-09-1-python-baseline.md) |
| Published-artifact output audit (P0.1) | Python/npm/Rust reran 10,000 chemistry rows, 310,000 SMARTS cells and 57 legacy reactions. Published Python/Rust cover 63 operation outputs; only HBA and its bundle changed from v1.0.29. The audit is complete, **not** strict RDKit parity: 200 SMARTS cells fail, five CIP rows abstain and four npm operations are unexposed. | [Artifact packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md) · [adoption policy](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md) |
| RDKit-compatible HBA | The hash-verified published v1.0.30 macOS arm64 wheel matches RDKit 2026.03.6 on 5,000/5,000 exposed ChEMBL rows; published v1.0.29 matched 3,641/5,000. HBA only. | [Release-channel evidence](https://github.com/kent-tokyo/chematic/blob/main/validation/results/release-channel-verification-v1.0.30.json) |
| Paired speed (P0.2) | Published Python: 20 operations pass exact-output and paired-interval gates in 20 alternating blocks on one host. Isolated Python and Rust lanes split parse-inclusive, first-use and precomputed calls. Whole-process RSS is not library allocation; changed HBA output is not an equivalent-output speed win. | [Python matrix](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) · [isolated Python](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1.0.30-isolated-python-time-memory.md) · [Rust](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1029-v1030-published-rust-isolated-time-memory.md) |
| Published browser Morgan | On Ubuntu 24.04, Chromium, Firefox and WebKit each have 250/250 direct and prepared bit parity and favorable 20-block speed intervals versus official RDKit.js (smallest lower bound 1.63×). A separate Chromium 10k lane has 9,999 exact and one typed Fe refusal. No universal browser or memory claim. | [Three-browser record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md) · [M4 Chromium](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.md) |
| Reactions (P1) | Legacy 57/57 fixtures match. Published v1.0.30 Rust: 73/83 graph/origin/map exact. PR #755 release-profile **source** wheels on Linux/macOS: 80/83 exact and three jointly invalid; the WASM Node 83-row test passed. Published v1.0.34 packages need their own rerun. No general SMIRKS parity claim. | [Published packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md) · [PR #755 and CI](https://github.com/kent-tokyo/chematic/pull/755) |
| Release channels | v1.0.34 GitHub Release, npm, PyPI, 20 crates.io archives, docs.rs and Pages were directly verified. This confirms availability, not chemical accuracy. | [v1.0.34 record](https://github.com/kent-tokyo/chematic/blob/main/validation/results/release-channel-verification-v1.0.34.json) |

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
  Linux/Python 3.9 published wheels have typed stereo failures on rows
  53/246; unreleased source, which takes its math from the `libm` crate,
  passes 265/265 on Linux ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-stereo-integrity-smarts-a6-followups.md)).
  A macOS rerun of that build and independent conformer quality remain open. [Quality record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-a6-published-v1031-mmff94-quality.md)
  · [platform diagnosis](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-a6-source-mmff94-termination.md).

Earlier source diagnostics and per-release channel checks remain in the
[benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md) and
[versioned validation results](https://github.com/kent-tokyo/chematic/tree/main/validation/results/). They are not promoted
to v1.0.34 measurements.

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
