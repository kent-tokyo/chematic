# v1.1.3 evaluation follow-up

This record measures the source candidate at commit `90542885` after the
v1.1.3 evaluation follow-up. It addresses the reported Rust source break,
descriptor/fingerprint residuals, and ordinary-input speed regressions. It is
source-candidate evidence, not a new published-package result.

## Correctness

The frozen 4,072-row unusual-SMILES corpus was compared with RDKit 2026.03.6.

| Operation | Result |
|---|---:|
| RDKit molecular weight | 4,072 / 4,072 exact |
| RDKit-compatible Morgan R2/2048 | 4,072 / 4,072 exact |
| RDKit-compatible ring count | 4,072 / 4,072 exact |
| Labute ASA | 3,776 exact; all 4,069 finite comparable rows within `1e-9`; 3 rows are NaN in both engines |
| CIP | 4,062 / 4,072 exact |

The ten CIP differences are classified rather than silently counted as
parity: seven are square-planar tags for which RDKit emits tetrahedral R/S,
two contain wildcard ligands and now return typed `oracle_unstable`, and one
is an isotope case where direct Sequence Rule 2 priority is applied before
remote Rule 5 pseudoasymmetry. No guessed wildcard label is emitted.

The machine-readable summary is
[`validation/results/v1.1.3-evaluation-fixes-stress-summary.json`](../validation/results/v1.1.3-evaluation-fixes-stress-summary.json).

## Source compatibility

`Chirality::Allene1` and `Chirality::Allene2` no longer add enum variants.
They are source-compatible associated constants backed by the existing two
parity variants; molecular context identifies allene geometry. An integration
test compiles the exact exhaustive v1.1.2-style match, preventing another
patch-release `E0004` break.

## Performance

The timings use macOS 27.0.1 arm64, Python 3.13.6, the first 1,000 rows of
`rdkit-js-browser-10k-v1.smi` (SHA-256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`),
20 fresh-process AB/BA blocks, one untimed warm-up, and output-digest equality.

| Lane | A / B median | Paired bootstrap 95% interval | Interpretation |
|---|---:|---:|---|
| public v1.1.2 / candidate parse | 0.9728 | 0.9254–1.0162 | statistically tied |
| public v1.1.2 / candidate parse + TPSA | 0.9719 | 0.9084–1.0209 | reported 21% regression reduced to a statistical tie |
| public v1.1.2 / candidate parse + Morgan | 0.9826 | 0.9304–1.0061 | reported 10% regression reduced to a statistical tie |
| public v1.1.3 / candidate prepared QED | **1.1609** | **1.1319–1.1779** | candidate is 16.1% faster |

The public v1.1.2 wheel SHA-256 is
`89b7269d2cc7051ad43ba9915ef7ce00813923c94724c7c7cfcabb8e7f354cda`;
the public v1.1.3 wheel is
`5569b170b037108d13e9af51386173050ac3b4ac4c798bf706d2a522ddd0e17b`;
the candidate wheel is
`01cca77a8a647c83a4d8c28fc6a250d4c28479902a7669a7547c6e8017622a94`.
Raw paired observations are stored beside this record.

The QED gain comes from reusing an already-compatible graph. The parse lanes
avoid a whole-molecule allene scan when no allene geometry is possible and
replace a complete SSSR construction in the aromatic-macrocycle shortcut with
a local shortest-cycle search.

## Verification

- all tests for `chematic-core`, `chematic-smiles`, `chematic-perception`,
  `chematic-cip`, `chematic-fp`, and `chematic-chem` passed;
- Python wheel: 1,036 passed, 77 skipped;
- Clippy passed for all changed Rust targets with warnings denied;
- formatting and `git diff --check` passed.

Hosted Linux, Windows, WASM/Node, and published-artifact gates remain CI work.
