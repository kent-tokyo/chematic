# A6 RDKit 2026.09 bounds and conformer follow-up

This source-candidate record asks two bounded questions:

1. Does CheMatic construct the same ETKDG topology bounds as RDKit 2026.09.1?
2. With 20 candidates and the same 200-iteration MMFF94 limit, how close is
   CheMatic's best conformer to RDKit's?

It does not establish universal conformer superiority or a published-package
speed win.

## Provenance

- CheMatic commit: `7e6a524054cc5950c858bd7acb11405da03988a5`
- Source-built wheel: `chematic-1.1.1-cp313-cp313-macosx_11_0_arm64.whl`
- Wheel SHA-256:
  `0ad3af2781f92a9de20d90ce9b2728250c3261a1194a2b4bcee787078d7c3508`
- Comparator: RDKit `2026.09.1`
- Host: macOS 27.0.1 arm64, Python 3.13.6
- Corpus: frozen A6 tiers A+B, 65+200 molecules
- Tier-A manifest SHA-256:
  `6a478ea0f5d4ef067a4d1739e77a7209e8f76ecaa837e3f487c723dd6f465d6b`
- Tier-B manifest SHA-256:
  `ec93ac160aa91c3c4416d57d3000d6e344e0b75cdf5c122211364c7495dce239`

The implementation adds RDKit 2026.09.1 angle-derived 1-3 bounds, 1-4
interval merging, relevant-ring handling, macrocycle amide behavior, forced
amide/ester torsions, and the ordinary-bond covalent-radius fallback. The
boundary fixture contains 776 focused cases.

## Results

| Gate | CheMatic versus RDKit 2026.09.1 |
|---|---:|
| Raw bounds, molecules within `1e-12` | 265/265 |
| Raw bounds, cells within `1e-12` | 768,664/768,664 |
| Raw bounds, maximum absolute delta | `1.78e-15` |
| Triangle-smoothed bounds, molecules within `1e-12` | 265/265 |
| Triangle-smoothed bounds, cells within `1e-12` | 768,664/768,664 |
| Triangle-smoothed bounds, maximum absolute delta | `1.42e-14` |
| Conformer rows compared | 265 |
| Failures | 0 |
| Best-energy median delta | `-1.50e-15` kcal/mol |
| Best-energy p90 delta | `+0.128` kcal/mol |
| Best-energy maximum delta | `+3.111` kcal/mol |
| At or below RDKit | 168/265 |
| Within 1 kcal/mol | 258/265 |
| Within 2 kcal/mol | 263/265 |

Seven rows remain more than 1 kcal/mol above RDKit. Six have zero converged
candidates in both engines under the common 200-iteration budget. The seventh
has 9/20 converged candidates in both engines. These rows are a basin-diversity
or minimization-tail problem, not evidence of a same-coordinate MMFF energy
formula mismatch.

One timing pass took 51.88 s for CheMatic and 56.39 s for RDKit, a nominal
1.087x ratio. The host was running unrelated coverage jobs and the benchmark
has only one alternating block, so this is diagnostic data, not an accepted
speed claim.

## Reproduce

After installing the recorded wheel and RDKit 2026.09.1 in an isolated Python
environment:

```bash
python scripts/a6_rdkit_bounds_matrix_parity.py \
  --tiers AB --no-triangle-smoothing --output /tmp/a6-bounds-raw.json
python scripts/a6_rdkit_bounds_matrix_parity.py \
  --tiers AB --output /tmp/a6-bounds-smoothed.json
python scripts/a6_rdkit_compatible_etkdg_mmff.py \
  --tiers AB --conformers 20 --max-mmff-iters 200 --repeats 1 \
  --output /tmp/a6-quality.json
```

The compact machine-readable result is
[`2026-10-12-a6-rdkit-2026-09-bounds-followup-summary.json`](2026-10-12-a6-rdkit-2026-09-bounds-followup-summary.json).

## Remaining acceptance work

- close the seven greater-than-1-kcal/mol rows without asymmetric work;
- evaluate a previously unused holdout;
- run at least 21 alternating blocks on an idle host and a second host;
- record peak memory and exact public PyPI/npm/crates.io artifacts.
