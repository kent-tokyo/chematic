# v1.0.42 regression follow-up — 2026-10-10

This record separates the published v1.0.42 baseline from commit `549338b2`,
which is an unpublished source candidate. It addresses the false ring counts,
explicit-hydrogen reaction residuals, Rust API compatibility, and stale
published-artifact gates identified after v1.0.42.

## Fixed inputs

- Published artifact: PyPI macOS arm64 CPython 3.13 wheel
  `chematic-1.0.42-cp313-cp313-macosx_11_0_arm64.whl`, SHA-256
  `1f9a0e7878e2b3957e3d0dae7c70eaa5ea5218efca6141e2285cba63301a71f0`.
- Source candidate: wheel built from `549338b2`, SHA-256
  `f6ad4dfdd724b053e3e39fbc0edbf379147be6b1424f179e4f66e05e408402f3`.
- Comparator: RDKit 2026.03.6.
- Host: macOS arm64, CPython 3.13.6.
- Published expected values:
  [`published-wheel-chemistry-gates-expected-v1.0.42.json`](../validation/published-wheel-chemistry-gates-expected-v1.0.42.json).

The source wheel is not a published-package result. Linux, Windows, npm/WASM,
and crates.io publication lanes remain separate checks.

## Published v1.0.42 baseline

The published wheel retains 80 exact reaction rows plus three inputs rejected
by both engines, 310,000/310,000 SMARTS cells, the recorded exposed and ChEMBL
CIP/hybridization/MMFF counts, and lossless MOL stereo round trips. Its full
BioTransformer explicit-hydrogen lane contains 6,486 exact rows, 153 differing
rows, and 165 rows whose raw RDKit products cannot be sanitized.

The same wheel passes all 265 A6 rows: all are independently sound and
stereo-clean, none has a gross clash, and the coordinate digest is
`b1ce930f542ee3ac7a674acadc4caa43bc2d63df570381557db76fb5fd9de804`.

## Ring-count regression

v1.0.42 treated one SMILES closure pair as one ring even when the closure
created only a coordination or disconnected-component connection. The source
candidate verifies that removing the closure bond leaves a path between its
endpoints before using the one-closure shortcut.

| SMILES | v1.0.42 | source candidate | RDKit |
|---|---:|---:|---:|
| `[Fe]1(Cl)(Cl)Cl.Cl1` | 1 | 0 | 0 |
| `C=C(C)N1C=NC<-C1` | 1 | 0 | 0 |
| `N1CCN->[Cu]1` | 1 | 0 | 0 |
| `C1CCN->[Cu]1` | 1 | 0 | 0 |

Rust, Python, and WASM fixtures cover these cases. Mutated molecules do not use
the parser shortcut, and explicit dative/query bonds retain the general graph
path. A diagnostic 1,000-molecule run with 20 calls per molecule showed no
measurable normal-corpus regression: the two public-wheel medians were 22.588
and 22.094 ms; the two source-candidate medians were 22.242 and 22.767 ms. This
short run is not a performance claim.

## Explicit-hydrogen reaction residuals

The RDKit compatibility profile had suppressed caller-supplied hydrogen graph
atoms when any explicit hydrogen was present. A deletion/edit case could then
emit a carbon radical. The source candidate preserves those graph atoms and
adds overlapping residual facets for deletion, bond breaking, radicals, and
re-sanitization or valence refusal.

| Fixed corpus | Metric | v1.0.42 | source candidate |
|---|---|---:|---:|
| 500 templates / 4,000 applications, seed 1039 | differing explicit-H rows | 52 | 19 |
| 500 templates / 4,000 applications, seed 1039 | exact explicit-H rows | 239 | 290 |
| 2,000 templates / 16,000 applications, seed 754 | differing explicit-H rows | 141 | 70 |
| 2,000 templates / 16,000 applications, seed 754 | exact explicit-H rows | 1,051 | 1,184 |
| full BioTransformer lane | differing explicit-H rows | 153 | 35 |
| full BioTransformer lane | exact explicit-H rows | 6,486 | 6,604 |

The implicit-hydrogen lane and the other chemistry gates are unchanged. The
remaining 35 BioTransformer differences and 70 generated-property differences
are classified residuals, not parity claims. Expected generated-corpus counts
and facets are fixed in
[`smirks_property_fuzz_expected.json`](../validation/smirks_property_fuzz_expected.json).

## Public API and failure contracts

- `TracedReactionTransformReport` again has its original two public fields, so
  downstream Rust struct literals compile. Rejected RDKit-profile products are
  available through a separate opt-in detailed report.
- Balaban J, IPC, and legacy MMFF94 BCI charges have checked entry points with
  typed unsupported/size errors. Legacy calls retain their previous sentinel
  behavior for source compatibility.
- Published-wheel workflow defaults and A6 pins now point to v1.0.42. The old
  v1.0.38 expected file remains historical evidence.

## Verification and boundary

- Rust library suites for `chematic-chem`, `chematic-rxn`, and
  `chematic-wasm` pass; the downstream struct-literal integration test passes.
- The source wheel passes 1,105 Python tests.
- Relevant clippy and formatting checks pass.
- The source chemistry packet retains 80+3 reaction rows, 310,000 SMARTS
  cells, and the published CIP/hybridization/MMFF/MOL counts.

This record does not replace the published v1.0.42 evidence with source output,
does not claim complete SMIRKS parity, and does not establish a cross-platform
performance result. Publication and other-host reruns must name their artifact
hashes separately.
