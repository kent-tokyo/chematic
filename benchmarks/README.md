# chematic benchmark records

This directory stores dated benchmark evidence. Every result is scoped to the
named source or package, comparator, corpus, host, runtime and timed boundary.
The methodology and interpretation rules are in
[`docs/benchmark.md`](../docs/benchmark.md); the complete record inventory
is in [`archive-index.md`](archive-index.md).

## Current evidence

| Question | Record | Boundary |
|---|---|---|
| Are the former Labute ASA and first-use ring-count deficits closed? | [v1.1.0 source follow-up](2026-10-12-v1.1.0-source-labute-ring-speed.md) · [raw JSON](2026-10-12-v1.1.0-source-rdkit-labute-ring-speed.json) | Yes on the recorded macOS arm64 source-built wheel: both are 10,000/10,000 equivalent, with prepared medians of 1.71x and 1.22x versus RDKit 2026.09.1 and favorable paired intervals. Not a published-package or cross-platform claim. |
| How does current `main` compare with RDKit? | [v1.1.0-source / RDKit 2026.09.1 Python packet](2026-10-11-v1.1.0-main-vs-rdkit-2026.09.1.md) · [focused raw JSON](2026-10-11-v1.1.0-main-rdkit-python-macos-arm64.json) · [63-operation raw JSON](2026-10-11-v1.1.0-main-vs-rdkit-63op-paired20.json) | Source-built wheel at `44338e16`; four 10,000-row output gates plus 21 fresh-process blocks, and a 63-operation/20-block supplementary matrix on macOS arm64. Not a published-package result. |
| What is the latest published-package comparison? | [v1.0.42 / RDKit 2026.09.1](2026-10-11-v1042-published-rdkit-packet.md) | Exact PyPI/npm artifacts on hosted Linux x64 and macOS arm64. |
| What is the broad Python operation matrix? | [v1.0.30 63-operation matrix](2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) | Historical published-wheel result against RDKit 2026.03.6; only output-equivalent lanes with favorable confidence intervals count as wins. |
| How does browser/WASM compare? | [v1.0.30 three-browser packet](2026-10-03-v1030-published-linux-three-browser-paired20.md) | Published npm artifacts, 250 exact Morgan rows and 20 blocks in Chromium, Firefox and WebKit. |
| How does similarity search compare? | [v1.0.13 similarity-search packet](2026-09-12-similarity-search-a3-v1.0.13.md) | Separates native/native, RDKit-compatible/RDKit and cross-profile overlap. |
| How does file interchange compare? | [Open Babel file-I/O program](../docs/openbabel-file-io.md) · [same-process candidate](2026-10-10-openbabel-file-io-same-process-v1.0.42.md) | Small pinned fixtures and source candidate; not broad-corpus or published-package evidence. |
| How do CDK and Indigo compare? | [bounded CDK/Indigo baseline](2026-10-10-cdk-indigo-bounded-baseline.md) | 24-row common matrix and 21 timing blocks; source-candidate smoke evidence. |
| How do reaction, parent selection, and R-group workflows compare with Indigo? | [Indigo workflow packet](2026-10-12-indigo-workflows-source.md) · [raw JSON](2026-10-12-indigo-workflows-source.json) | Source-built, one-host bounded lanes. CheMatic has higher exact counts for reaction and parent selection; R-group is tied. Common-correct medians favor CheMatic in all three lanes. This is not a toolkit-wide or published-package claim. |
| How does COSMolKit compare? | [RDKit/COSMolKit 0.5 comparison](2026-10-10-cosmolkit-050.md) | Source-built 161-operation and unusual-SMILES evidence; not a published-package claim. |

## What counts as a speed win

A lane is called an equivalent-output speed win only when:

1. both engines perform the same declared work;
2. every compared output satisfies the lane's exact or numeric contract;
3. failures and unsupported rows are fully accounted for;
4. the lower bound of the paired 95% confidence interval is above `1.0`;
5. source, package, corpus, runtime and host provenance are recorded.

Parse-inclusive and prepared-input timings answer different questions. A
cached RDKit property lookup is not compared as though it included perception,
and a source-built result is not presented as a registry-package result.

## Reproduce the current Python packet

Build a wheel from the recorded revision and install it with the pinned RDKit
wheel in an isolated Python 3.13 environment. Then run:

```bash
python scripts/run_published_v1042_rdkit_python_packet.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --artifact-dir /path/to/exact-wheels \
  --output /tmp/chematic-rdkit.json \
  --expected-chematic 1.1.0 \
  --expected-rdkit 2026.9.1 \
  --schema rdkit-labute-ring-optimization/v1 \
  --blocks 41 \
  --chematic-revision f2d030fc6a2760fc9ad3f168f98713584c20828c
```

The runner checks all 10,000 outputs before emitting speed intervals and uses
41 alternating fresh-process blocks over the first 1,000 rows.

Validate the generated record with:

```bash
python scripts/check_published_v1042_rdkit_packet.py /tmp/chematic-rdkit.json \
  --expected-schema rdkit-labute-ring-optimization/v1 \
  --expected-chematic 1.1.0 \
  --expected-rdkit 2026.9.1 \
  --expected-revision f2d030fc6a2760fc9ad3f168f98713584c20828c \
  --expected-blocks 41 \
  --minimum-artifacts 1
```

## Record maintenance

- Add new claims as a Markdown summary plus machine-readable raw data.
- Keep superseded records; label them historical instead of rewriting them.
- Link every checked-in root-level `.md` and `.json` record from this file or
  [`archive-index.md`](archive-index.md).
- Validate the indexes with `python3 scripts/check_benchmark_index.py`.
- Do not infer universal accuracy or speed from one corpus, host or API lane.

The [complete inventory](archive-index.md) contains the full dated record
list, streaming/format matrices, 3D evidence, artifact audits and earlier
candidate measurements.
