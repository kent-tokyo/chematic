# RDKit-compatible multi-conformer ETKDG/MMFF94 on A6

This source-candidate experiment evaluates a deterministic multiple-conformer
path that shares ETKDG setup and one immutable MMFF94 parameterization across
candidates. Native builds process candidates in parallel; WASM uses the same
candidate order sequentially.

The implementation targets RDKit 2026.03.1 ETKDGv3 coordinates and seed
progression. It is compared here with RDKit 2026.09.1 as the current quality
reference. A partial attempt to expose 2026.09 bounds was withdrawn because
only its 1-3 changes had been ported; no 2026.09 coordinate-parity claim is
made.

## Provenance

- Base revision: `035f479ecfd5aef10ebba89b862ef73471662813`
- Candidate: uncommitted source diff on branch `codex/a6-multiconformer-etkdg`
- Final validation wheel SHA-256:
  `afd40db7c3573bc6787c890e67119c7e53a76e75fb292214bffc7f78177a19a1`
- Wheel: `chematic-1.1.1-cp313-cp313-macosx_11_0_arm64.whl`
- Comparator: RDKit `2026.09.1`
- Host: macOS 27.0.1 arm64, CPython 3.13.6, 10 logical CPUs
- Corpus: frozen A6 tier A+B, 265 molecules
- Tier A SHA-256: `6a478ea0f5d4ef067a4d1739e77a7209e8f76ecaa837e3f487c723dd6f465d6b`
- Tier B SHA-256: `ec93ac160aa91c3c4416d57d3000d6e344e0b75cdf5c122211364c7495dce239`
- Candidate count: 20 per molecule
- Initial MMFF budget: 200 iterations per candidate

The concise machine-readable result is
[`2026-10-12-a6-multiconformer-etkdg-mmff-summary.json`](2026-10-12-a6-multiconformer-etkdg-mmff-summary.json).

## Quality

| Lane | Median CheMatic - RDKit | p90 | Maximum | At/below RDKit | Within 1 kcal/mol | Within 2 kcal/mol |
|---|---:|---:|---:|---:|---:|---:|
| Baseline, 20 x 200 iterations | -0.00000000017 | 0.859 | 5.769 | 175/265 | 243/265 | 253/265 |
| Adaptive, lowest five non-converged + 2,000 iterations | -0.0000000133 | 0.138 | 5.637 | 211/265 | 256/265 | 258/265 |

Adaptive retry changed 112 rows in the measured result; all changed rows had
lower or equal returned energy. It raises converged CheMatic candidates from
2,835 to 3,602. The remaining maximum residual shows that neither lane beats
RDKit on every molecule.

## Throughput

The two nearby one-repeat baseline runs produced opposite winners:

| Run | CheMatic | RDKit | RDKit / CheMatic |
|---|---:|---:|---:|
| A | 80.233 s | 86.089 s | 1.073x |
| B | 79.761 s | 72.250 s | 0.906x |

The adaptive run took 96.269 s versus RDKit's 75.325 s (0.782x). These runs
do not satisfy the repository's 21 alternating-block, confidence-interval,
two-host, memory, or exact-public-artifact gate. Therefore this record makes
no accepted speed-win claim.

## API and regression checks

- Rust exposes deterministic multi-conformer ETKDG, combined ETKDG/MMFF94,
  and opt-in bounded adaptive refinement.
- Python exposes the same paths through `Mol.rdkit_embed_multiple_confs()` and
  `Mol.rdkit_embed_multiple_confs_mmff()`.
- Native execution uses deterministic candidate-order parallel collection;
  browser WASM uses the sequential fallback.
- Frozen 2026.03 seed progression, repeated-call determinism, shared-field
  equivalence, adaptive monotonicity, and Python binding behavior have direct
  regression tests.

## Reproduce

Build and install the source wheel in an isolated environment containing the
pinned RDKit comparator, then run:

```bash
python scripts/a6_rdkit_compatible_etkdg_mmff.py \
  --output /tmp/a6-etkdg-mmff.json \
  --conformers 20 \
  --repeats 21

python scripts/a6_rdkit_compatible_etkdg_mmff.py \
  --output /tmp/a6-etkdg-mmff-adaptive.json \
  --conformers 20 \
  --retry-candidates 5 \
  --retry-max-iters 2000 \
  --repeats 21
```

Before promotion, repeat from exact public artifacts on two hosts, record peak
memory, evaluate an unused holdout, and either close the 5.637 kcal/mol tail or
report it as a typed quality boundary.
