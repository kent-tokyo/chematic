# chematic benchmark records

This directory stores versioned evidence, not timeless product claims. Every
result is limited to its named artifact or source revision, comparator, corpus,
host, runtime, options, and timed boundary.

Start with the [benchmark method](../docs/reference/benchmark.md). Use the
[complete inventory](archive-index.md) only when you need a historical record.

## Current entry points

| Question | Best entry | Boundary |
|---|---|---|
| Did the current public Python artifacts pass chemistry gates? | [v1.2.0 hosted chemistry run](https://github.com/kent-tokyo/chematic/actions/runs/39511354905) | Eleven hosted jobs; pass/fail evidence, not a speed result or retained row packet |
| What is the latest durable cross-platform wheel packet? | [v1.1.3 artifact packet](2026-10-14-v113-published-artifact-packet.md) | Four exact PyPI wheels; later hosted runs and npm are separate |
| What is the broadest published RDKit speed packet? | [v1.0.42 comparison](2026-10-11-v1042-published-rdkit-packet.md) | Historical exact artifacts; do not transfer ratios to current releases |
| What public-wheel descriptor result is reproducible? | [v1.1.8 Labute ASA](2026-10-19-v118-published-labute-asa-cosmolkit.md) | One macOS arm64 wheel, host, descriptor, and corpus |
| What are the latest standardization counts? | [Cleanup](2026-10-22-rdkit-cleanup-metal-aromaticity.md) · [Reionize](2026-10-22-rdkit-reionize-in-memory-stereo.md) | Source evidence against RDKit 2026.09.1; not v1.2.0 package evidence |
| Is MCS input-order independent and are incomplete searches explicit? | [stereo parity](2026-10-22-mcs-stereo-neighbor-parity.md) · [resource status](2026-10-21-mcs-non-induced-order-status.md) | Focused source cases; broad FMCS parity remains open |
| Does a relaxed-comparator MCS query re-match every input? | [observed query alternatives](2026-10-25-rdkit-mcs-query-alternatives.md) · [Python/WASM SMARTS output](2026-10-25-rdkit-mcs-smarts-output.md) | Source API and focused round trips; charge/isotope/chirality constraints and broad FMCS parity remain open |
| Did the v1.2.1 priority fixes survive unused data? | [priority follow-ups and holdout v2](2026-10-24-v121-priority-followups.md) · [machine summary](2026-10-24-v121-priority-followups.json) | MCS speed and Cleanup parity improved; frozen MCS/3D residuals remain explicit |
| What is the current reaction/R-group boundary? | [workflow packet](2026-10-18-rdkit-workflows-followup.md) · [automatic columns](2026-10-22-rdkit-rgroup-auto-columns.md) | Bounded fixtures; symmetric-ring residuals remain |
| What is the current 3D/MMFF boundary? | [adaptive stereo follow-up](2026-10-25-rdkit-search-adaptive-stereo.md) · [holdout v2](2026-10-24-v121-priority-followups.md) · [MMFF94 atom typing](2026-10-24-mmff94-fullerene-typing-rdkit.md) · [torsion-search follow-up](2026-10-24-etkdg-torsion-search-rdkit.md) · [alignment and ETKDG follow-up](2026-10-24-3d-alignment-etkdg-followup.md) | v2 separates embedding success from coordinate parity; source MMFF94 types match RDKit on fixed corpora, while broad package and equal-compute claims remain open |
| Are the current fingerprint and USR residuals resolved? | [search and fingerprint regression](2026-10-25-rdkit-search-adaptive-stereo.md) · [aromatic ring triple and same-coordinate USR](2026-10-25-rdkit-core-usr-accuracy.md) · [machine fingerprint packet](2026-10-25-rdkit-core-usr-accuracy.json) | Five 2048-bit fingerprints exact on three fixed corpora; USR matches on 24 fixed coordinate sets, not conformer-quality evidence |
| How does file interchange compare with Open Babel? | [comparison scope](../docs/reference/openbabel-file-io.md) · [SDF follow-up](2026-10-14-openbabel-sdf-read-followup.md) | Semantic equality is required before speed ranking |
| How does browser/WASM Morgan compare? | [three-browser packet](2026-10-03-v1030-published-linux-three-browser-paired20.md) | Published v1.0.30 npm artifacts on one hosted environment |
| Where are COSMolKit comparisons? | [0.5 comparison](2026-10-10-cosmolkit-050.md) · [seeded ETKDG](2026-10-09-cosmolkit-parity-2.md) | Source-built operation evidence, not universal superiority |

The [validation report](../docs/reference/validation.md) summarizes current
correctness boundaries. The inventory keeps older experiments discoverable
without making this page a second changelog.

## What counts as a speed win

A lane is ranked only when:

1. both engines perform the same declared work;
2. outputs satisfy the stated exact, numeric, or semantic contract;
3. failures, refusals, unsupported rows, and skipped work remain visible;
4. repeated paired measurements support the result; and
5. artifact, corpus, runtime, host, and raw samples are recorded.

Parse-inclusive, first-use, and prepared/reuse timings answer different
questions. Source results are never relabelled as registry-package results.

## Retention policy

Keep current release evidence, frozen inputs, compact machine-readable results,
and records referenced by CI, tests, release metadata, or active docs. Remove
copied wheels, source archives, build trees, and large intermediate output once
their summarized experiment is retired. On 2026-10-23, 88 unreferenced
intermediate records through 2026-10-09 were removed after confirming that all
of them remain available in the immutable `v1.0.41` tag. Git history and release
tags, rather than the active directory, are the archive for retired material.

A retained historical record is not evidence about current source.

## Adding or changing a record

- Put the conclusion, limits, exact command, revisions, artifact/corpus hashes,
  host, and validation rule in the dated Markdown record.
- Keep compact raw data beside the summary when it is needed to reproduce the
  calculation.
- Link every retained root-level record from this file or `archive-index.md`.
- Run `python3 scripts/checks/check_benchmark_index.py` after changes.
- Do not infer universal accuracy or speed from one corpus, host, or API lane.
