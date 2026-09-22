# Public-package fingerprint and 3D comparison — 2026-09-22

This record separates released artifacts from local source candidates. It does
not establish that chematic is universally faster than RDKit.

## Environment and fixed inputs

- Host: macOS 26.5.2, arm64; browser lane uses Google Chrome/Chromium.
- Published artifacts: `@kent-tokyo/chematic@1.0.19`,
  `chematic==1.0.19`, `@rdkit/rdkit@2026.3.6`, and `rdkit==2026.3.6`.
- Fingerprints: fixed exposed 10,000-SMILES corpus, radius 2, 2,048 bits,
  packed 256-byte LSB-first output, 20 fresh browser processes per arm.
- 3D: fixed 265-molecule Tier A+B corpus, seed `20260801`, at most 8 embedding
  attempts, and the same external geometry/stereo scorer.

## Fingerprint result

The published 1.0.19 npm package already wins the parse-inclusive lane, but its
old prepared lane repeats compatibility preprocessing and loses to RDKit.js.
The source candidate adds an immutable prepared handle and a conservative
linear aromaticity fast path.

| Artifact / operation | chematic ms/mol (p50 process mean) | RDKit.js ms/mol | Geometric speedup | 95% lower bound | Result |
|---|---:|---:|---:|---:|---|
| Published 1.0.19, parse + FP | 0.036654 | 0.098465 | 2.658x | 2.617x | pass |
| Published 1.0.19, prepared FP | 0.034148 | 0.023857 | 0.687x | 0.671x | fail |
| Source candidate, parse + prepare + FP | 0.075238 | 0.105301 | 1.397x | 1.333x | pass |
| Source candidate, prepared FP | 0.007106 | 0.024427 | 3.551x | 3.384x | pass |

For the source candidate, all 20 paired repetitions favored chematic: the
minimum observed speedups were 1.114x parse-inclusive and 2.771x prepared.
The prepared API matched RDKit.js bit-for-bit on all 9,999 supported rows. One
high-coordinate Fe(II) row remains the same explicit typed refusal. The source
candidate is not a registry release and must be rerun from the eventual
published package before becoming a release claim.

Evidence:

- `validation/results/competitive-browser-rdkitjs-published-npm-10k-chromium-20rep-prepared-v1.0.19-2026-09-22.json`
- `validation/results/competitive-browser-rdkitjs-source-candidate-10k-chromium-20rep-prepared-2026-09-22.json`
- `validation/results/competitive-browser-rdkitjs-ecfp4-prepared-parity-10k-source-candidate-2026-09-22.json`

## 3D result

Published 1.0.19 has scoped speed wins, but not a general 3D win:

| Lane | chematic usable | RDKit usable | Geometric speedup | 95% lower bound | Interpretation |
|---|---:|---:|---:|---:|---|
| Raw embedding | 212/265 | 264/265 | 3.287x | 2.971x | faster, inadequate stereo coverage |
| UFF | 261/265 | 264/265 | 2.733x | 2.496x | faster, three fewer usable rows |
| UFF best-of-10 | 262/265 | 264/265 | 1.878x | 1.745x | faster, two fewer usable rows |
| MMFF94 stereo-safe | 264/265 | 264/265 | 0.0689x | 0.0610x | equal usable count, about 14.5x slower |

The UFF source candidate fixes the extreme finite-gradient step problem. Its
best-of-10 lane is usable on 265/265 inputs versus RDKit's 264/265 and is
2.327x faster geometrically (95% lower bound 2.168x). This is candidate-only
evidence. Single-conformer UFF still has two typed failures, and MMFF94 speed
remains open. The fingerprint and 3D candidate artifacts have separate recorded
source-diff and binary/wheel hashes; this report does not pretend they were one
published package.

Evidence:

- `validation/results/public-package-3d-summary-comparable-v1.0.19-vs-rdkit-2026.3.6-2026-09-22.json`
- `validation/results/public-package-3d-summary-stereo-safe-v1.0.19-vs-rdkit-2026.3.6-2026-09-22.json`
- `validation/results/source-candidate-3d-summary-uff-vs-rdkit-2026.3.6-2026-09-22.json`

## Decision

The source candidate demonstrates a scoped RDKit.js speed win for both
parse-inclusive and correctly prepared compatible Morgan, and a scoped 3D UFF
best-of-10 speed/coverage win. General 3D superiority is not established.
Before publication: run full workspace/clippy gates, merge the candidate, then
install the new npm/PyPI artifacts from registries and repeat these exact
measurements.
