# MMFF94 stereo-safe source-candidate performance (2026-09-23)

Source candidate `e9f178fa` removes the measured MMFF94 stereo-safe throughput
deficit without relaxing the 300-iteration budget or the external geometry and
stereo judge. This is source-candidate evidence, not a published-package result
and not a claim of general 3D superiority.

## Before and after

| Artifact / corpus | CheMatic usable | RDKit usable | Geometric speedup | 95% lower bound |
|---|---:|---:|---:|---:|
| Published v1.0.19, fixed 265 | 264/265 | 264/265 | 0.0689x | 0.0610x |
| Published v1.0.20, fixed 265 | 265/265 | 264/265 | 0.944x | 0.861x |
| Source `e9f178fa`, fixed 265, three runs | 265/265 | 264/265 | **3.633–3.707x** | **3.336–3.402x** |
| Source `e9f178fa`, post-freeze holdout 100, three runs | 100/100 | 100/100 | **2.195–2.209x** | **1.998–2.014x** |

The historical v1.0.19 row is the approximately 14.5x deficit that motivated
this work. The v1.0.20 row is the current registry boundary. Only a future
registry rerun may promote the source-candidate rows into a public-package
claim.

## Frozen protocol

- Host: macOS 26.5.2 arm64; Rust 1.97.0; CPython 3.13.6.
- Comparator: `rdkit==2026.3.6` / runtime `2026.03.6`.
- CheMatic arm: `chematic_pipeline_v2_mmff94_strict_stereo_safe`, seed
  `20260801`, at most eight embedding attempts, MMFF94 maximum 300 iterations.
- RDKit arm: ETKDGv3 + MMFF94, the same seed and attempt limit, maximum 200
  force-field iterations.
- Quality: the same independent geometry/stereo scorer for both engines.
- Speed: paired common successes, geometric mean of
  `RDKit elapsed / CheMatic elapsed`; the lower 95% bound must exceed 1.0.
- The 265-row comparator is the already-pinned public RDKit result. The new
  RDKit holdout run was measured once and the frozen CheMatic binary three
  times sequentially. Confidence bounds are across molecule pairs, not
  process-level repetition or cross-host variance.

Candidate `e9f178fa7e9d969f594a6ace0a8712d3b0f472e2` was committed before the
final holdout was selected. The 100-row manifest contains 50 declared-stereo
and 50 general molecules selected from
`scripts/chembl_accuracy_corpus_4999.smi` (source SHA-256
`1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`).
Its corpus SHA-256 is
`503bffd4714817e1259256b14dff71c8059a71664f1fccfabeeb2b280b45fc48`.
Tier A, Tier B, and the first development holdout were excluded. Selection
used RDKit parse metadata and deterministic hashes only; CheMatic outcomes
were not available until after the manifest was frozen.

## What changed

- Final residual-force reporting now uses the already-validated prepared
  analytic MMFF94 gradient instead of an extra central-difference sweep that
  performed 6N full energy evaluations.
- Declared tetrahedral centers materialize only the implicit H required for
  stereo geometry. Multi-ring junction centers retain full-H expansion where
  the narrower representation did not preserve the existing quality gate.
- Molecules with at most 64 atoms scan the prepared van der Waals pair list
  directly. This avoids rebuilding spatial hash bins for every line-search
  energy and gradient probe; the cutoff and prepared pair set are unchanged.
- The L-BFGS direction reuses its coordinate buffer and fixed-size alpha
  storage. This is allocation cleanup, not the primary measured speedup.

## Result details

On the independent 100-row holdout, both engines were 100/100 successful,
independently sound, and stereo-clean. CheMatic p50 was 10.823–11.222 ms versus
RDKit's 23.5 ms; p95 was 43.876–45.926 ms versus 116.65 ms. The three lower
confidence bounds were 1.998x, 2.008x, and 2.014x.

On the fixed 265 rows, CheMatic remained 265/265 usable versus RDKit's 264/265.
CheMatic p50 was 10.557–10.818 ms versus RDKit's 28.732 ms. The three lower
confidence bounds were 3.402x, 3.336x, and 3.348x.

Not every molecule is faster: the minimum per-molecule RDKit/CheMatic ratio on
the holdout was 0.731–0.744. The acceptance claim is therefore aggregate for
the declared corpus and operation, not universal per-molecule dominance.

## Evidence

- Final manifest:
  `validation/manifests/mmff94_stereo_safe_post_freeze_holdout_v2.json`
- Holdout raw rows, common-scored rows, and summaries:
  `validation/results/a6-mmff94-stereo-safe-e9f178fa-holdout-v2-*`
- Fixed-265 candidate rows, common-scored rows, and summaries:
  `validation/results/a6-mmff94-stereo-safe-e9f178fa-fixed-265-*`
- Fixed-265 RDKit comparator:
  `validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.jsonl`

The post-freeze lane can be rerun with:

```bash
python scripts/gen_mmff94_stereo_safe_holdout_manifest.py \
  --candidate-commit e9f178fa7e9d969f594a6ace0a8712d3b0f472e2 \
  --source scripts/chembl_accuracy_corpus_4999.smi \
  --output validation/manifests/mmff94_stereo_safe_post_freeze_holdout_v2.json \
  --salt chematic-mmff94-stereo-safe-holdout-v2-final \
  --exclude-manifest validation/manifests/mmff94_stereo_safe_post_freeze_holdout_v1.json

python scripts/pipeline_v2_vs_rdkit_oracle.py \
  --manifest validation/manifests/mmff94_stereo_safe_post_freeze_holdout_v2.json \
  --only-arm rdkit_etkdgv3_mmff94

SCHEMATIC_MMFF94_MAX_ITERATIONS=300 \
  target/release/examples/pipeline_v2_vs_rdkit_dump \
  --manifest validation/manifests/mmff94_stereo_safe_post_freeze_holdout_v2.json \
  --only-arm chematic_pipeline_v2_mmff94_strict_stereo_safe
```

Use `pipeline_v2_vs_rdkit_common_scorer` for the shared judge and
`scripts/summarize_public_package_3d.py --pair mmff94-stereo-safe
chematic_pipeline_v2_mmff94_strict_stereo_safe rdkit_etkdgv3_mmff94` for the
paired summary. Redirect each raw run to a separate JSONL file; do not average
or discard failed rows before scoring.

## Decision and remaining boundary

The source candidate passes the scoped MMFF94 stereo-safe speed and quality
gate on both the fixed corpus and a post-freeze independent holdout. The A6
speed deficit can be closed for current source. The two same-coordinate energy
residuals, term-level adjudication, timeout behavior, broader conformer
quality, and public-package rerun remain open. No result here changes the
Experimental 3D support label.
