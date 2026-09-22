# A6 MMFF94 stereo-safe quality candidate (2026-09-23)

This packet evaluates source candidate `dd7fe3e9` against the pinned published
`rdkit==2026.3.6` wheel. It closes the previously reported 12 typed stereo
failures and four gross-clash rows on the fixed 265-molecule A/B corpus. It is
source-candidate evidence, not a result from a published CheMatic package.

> **Published confirmation:** the registry-installed v1.0.20 rerun preserves
> the 265/265 stereo-safe, independently sound and clash-free result. See
> [`2026-09-23-public-package-fingerprint-3d-v1.0.20.md`](2026-09-23-public-package-fingerprint-3d-v1.0.20.md).

## Change under test

The benchmark now exposes the production `PipelineV2Config::stereo_safe`
contract instead of treating `StereoPolicy::RepairAndVerify` alone as the
production quality lane. The pipeline retains expanded-hydrogen stereo during
MMFF94 line search when heavy-only and expanded stereo evaluation disagree.
Rejected trial steps are never committed; failure to find an admissible step is
reported as ordinary non-convergence.

The historical `Ignore` and `RepairAndVerify` arms remain unchanged for
comparison.

## Protocol

- Corpus: fixed tier A/B manifests, 65 + 200 molecules, seed `20260801`.
- CheMatic: source revision
  `dd7fe3e9b04fb143f82c8306f5079f92607a262d`, clean tree, release build,
  strict MMFF94, maximum 300 iterations, 25-second per-row hard timeout.
- Comparator: published `rdkit==2026.3.6`, runtime `2026.03.6`, ETKDGv3 plus
  MMFF94, using the already recorded wheel and result hashes.
- Judge: the same external geometry/stereo scorer for both engines. Every input
  remains in the denominator.
- Speed: paired common successes only; geometric mean of
  `RDKit elapsed / CheMatic elapsed`, with a required lower 95% bound above 1.0.
- Timing limitation: the two engines were measured on the same Apple-arm64 host
  family but not interleaved in one multi-repetition run.

## Results

| Measure | CheMatic stereo-safe | RDKit |
|---|---:|---:|
| successful rows | **265/265** | 264/265 |
| independently sound | **265/265** | 264/265 |
| stereo-clean | **265/265** | 264/265 |
| gross-clash success rows | **0** | 0 |
| p50 elapsed | 41.639 ms | 28.733 ms on successful rows |
| p95 elapsed | 227.495 ms | 162.375 ms on successful rows |

The paired speed ratio is **0.798x** with a **0.730x** 95% lower bound. The
stereo-safe quality lane therefore does not beat RDKit on speed in this packet.
The lighter `StereoPolicy::Ignore` source lane remains the separately measured
1.315x throughput result, but it is not stereo-quality equivalent.

This result closes the specific 12-failure/four-clash cohort. It does not close
A6: two same-coordinate energy residuals above 5 kcal/mol, term-level oracle
adjudication, timeout behavior, conformer-quality acceptance, and a published-
package rerun remain open.

## Evidence

- CheMatic rows:
  `validation/results/a6-mmff94-stereo-safe-dd7fe3e9-2026-09-23.jsonl`
- CheMatic provenance:
  `validation/results/a6-mmff94-stereo-safe-dd7fe3e9-2026-09-23.meta.json`
- Common-scored rows:
  `validation/results/a6-mmff94-stereo-safe-dd7fe3e9-vs-rdkit-2026.3.6-common-scored-2026-09-23.jsonl`
- Summary:
  `validation/results/a6-mmff94-stereo-safe-dd7fe3e9-vs-rdkit-2026.3.6-summary-2026-09-23.json`
- Comparator rows and provenance:
  `validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.{jsonl,meta.json}`

The CheMatic binary SHA-256 is
`97a29f3339914aa9bd96c7457c7890eafc4c15c060b49167a55ef168670f2633`.
The tier A/B manifest SHA-256 values are
`6a478ea0f5d4ef067a4d1739e77a7209e8f76ecaa837e3f487c723dd6f465d6b`
and `ec93ac160aa91c3c4416d57d3000d6e344e0b75cdf5c122211364c7495dce239`.
