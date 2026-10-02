# Validation report

Updated 2026-10-03. The current public release is **v1.0.30**. Every claim
below is limited to its pinned artifact, comparator, corpus and operation.
Source-only changes in [Unreleased](../CHANGELOG.md) are not package results.

## Current evidence

| Gate | Result and limit | Record |
|---|---|---|
| Published-artifact output audit (P0.1) | Python/npm/Rust reran 10,000 chemistry rows, 310,000 SMARTS cells and 57 legacy reactions. Published Python/Rust cover 63 operation outputs; only HBA and its bundle changed from v1.0.29. The audit is complete, **not** strict RDKit parity: 200 SMARTS cells fail, five CIP rows abstain and four npm operations are unexposed. | [Artifact packet](../benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md) · [adoption policy](../benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md) |
| RDKit-compatible HBA | The hash-verified published v1.0.30 macOS arm64 wheel matches RDKit 2026.03.6 on 5,000/5,000 exposed ChEMBL rows; published v1.0.29 matched 3,641/5,000. HBA only. | [Release-channel evidence](../validation/results/release-channel-verification-v1.0.30.json) |
| Paired speed (P0.2) | Published Python: 20 operations pass exact-output and paired-interval gates in 20 alternating blocks on one host. Isolated Python and Rust lanes split parse-inclusive, first-use and precomputed calls. Whole-process RSS is not library allocation; changed HBA output is not an equivalent-output speed win. | [Python matrix](../benchmarks/2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) · [isolated Python](../benchmarks/2026-10-03-v1.0.30-isolated-python-time-memory.md) · [Rust](../benchmarks/2026-10-03-v1029-v1030-published-rust-isolated-time-memory.md) |
| Published browser Morgan | On Ubuntu 24.04, Chromium, Firefox and WebKit each have 250/250 direct and prepared bit parity and favorable 20-block speed intervals versus official RDKit.js (smallest lower bound 1.63×). A separate Chromium 10k lane has 9,999 exact and one typed Fe refusal. No universal browser or memory claim. | [Three-browser record](../benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md) · [M4 Chromium](../benchmarks/2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.md) |
| Reactions (P1) | Legacy 57/57 fixtures match. Published Rust graph/origin/template-map identity matches RDKit on 73/83 stratified rows, plus one map-only difference. An unpublished checked-source profile has 76 matches, three typed unsupported, one refusal and three jointly invalid rows; it does not establish published or general SMIRKS parity. | [Artifact packet](../benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md) · [checked-source profile](../benchmarks/2026-10-02-reaction-checked-source-profile.md) |
| Release channels | GitHub Release, npm, PyPI, crates.io, docs.rs and Pages were independently verified for v1.0.30. A new release needs a new check. | [Channel record](../validation/results/release-channel-verification-v1.0.30.json) |

## Other bounded results

- **A0 core-eight:** frozen candidate `5e9211a6` passed 2,000/2,000
  development and a one-time 8,000/8,000 sealed holdout on eight declared
  descriptors after source and overlap checks. The raw sealed rows remain
  local-only and cannot be reused for tuning. [Commit-safe summary](../validation/results/a0-core-eight-sealed-acceptance-20260922.json).
- **CIP:** pinned exposed 10k comparison has 9,995 exact and five typed
  abstentions (four `oracle_unstable`, one `lone_pair_center`). This is
  neither a guessed label nor full CIP parity. [Source adjudication](../validation/results/rdkit-rebaseline-issue634-v1.0.27-candidate-vs-rdkit-2026.03.6-2026-09-28.json).
- **Morgan:** the published M4 Chromium lane matches on 9,999/9,999
  supported 10k inputs; one Fe(II) coordination input is a typed refusal.
  Native ECFP4 is a different profile. [Row-level record](../benchmarks/2026-10-03-v1030-published-chromium-morgan-row-parity.md).
- **3D/MMFF94:** Experimental. The historical v1.0.26 wheel quality packet
  retained 265 rows; neither same-coordinate energy nor a successful
  minimization proves broad conformer quality or speed parity.
  [Quality record](../benchmarks/2026-09-26-mmff94-public-v1.0.26.md).

Earlier source diagnostics and per-release channel checks remain in the
[benchmark index](../benchmarks/README.md) and
[versioned validation results](../validation/results/). They are not promoted
to v1.0.30 measurements.

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
