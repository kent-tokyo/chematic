# A6 analytic MMFF94 source candidate (2026-09-23)

This packet evaluates source candidate `af7c0c44` against the already-pinned
published RDKit `2026.3.6` wheel. It is current-source evidence, not a claim
about published CheMatic `1.0.19` and not yet a release result.

> **Follow-up:** candidate `dd7fe3e9` closes this packet's 12 typed MMFF94
> failures and four gross-clash rows in the production stereo-safe lane. See
> [`2026-09-23-a6-mmff94-stereo-safe-quality.md`](2026-09-23-a6-mmff94-stereo-safe-quality.md).
> The lighter speed lane below remains historical and is not quality-equivalent.

## Protocol

- Corpus: the fixed A/B 3D manifests, 265 molecules, seed `20260801`.
- Candidate: analytic prepared MMFF94 gradient, L-BFGS history 80, maximum 300
  iterations; UFF uses the same 300-iteration pipeline budget.
- Timeout: one subprocess per molecule, 25 seconds hard limit. All 265 inputs
  retain a terminal row and their original per-tier index.
- Comparator: `rdkit==2026.3.6`, ETKDGv3 plus MMFF94 or UFF, 200 force-field
  iterations. The wheel digest and runtime are recorded in the existing
  public-package metadata.
- Judge: the same independent geometry/stereo scorer for both engines.
- Speed gate: paired common successes; geometric mean of
  `RDKit elapsed / CheMatic elapsed`; the lower 95% bound must exceed 1.0.
- Timing: floating-point milliseconds. An earlier diagnostic used integer
  milliseconds and omitted sub-millisecond rows from the paired ratio; that
  packet was rejected and is not committed.

## Results

| Lane | CheMatic / RDKit successes | CheMatic / RDKit converged | paired geometric speedup | 95% lower bound | quality boundary |
|---|---:|---:|---:|---:|---|
| MMFF94, `StereoPolicy::Ignore` | 265 / 264 | 148 / 129 | **1.315x** | **1.227x** | 265/265 independently sound and zero gross-clash rows, but only 211/265 molecules are stereo-clean versus 264/265 for RDKit |
| MMFF94, `RepairAndVerify` | 253 / 264 | 142 / 129 | **1.318x** | **1.226x** | all 253 successes are stereo-clean; 12 typed failures and four successful rows with gross clashes remain, so quality non-inferiority is not established |
| UFF | 264 / 264 | 31 / 137 | **3.002x** | **2.771x** | both sides have 264/265 independently sound, stereo-clean successes; CheMatic has zero gross-clash success rows |

The MMFF94 candidate exceeds RDKit's measured throughput and convergence count,
but it does **not** yet establish a complete 3D win: the default ignore-stereo
lane is not quality-equivalent, while the repair lane loses 11 additional
inputs and retains four clash rows. UFF is the current lane where the measured
quality denominator is equal and the speed gate passes.

The convergence flags are engine-specific diagnostics, not a shared numerical
tolerance proof. Geometry/stereo quality therefore remains the acceptance
judge. Timings were collected on the same Apple-arm64 host family but on
different days and are not an interleaved multi-repetition benchmark; a
published-package rerun is still required.

## Evidence

- Candidate rows and provenance:
  `validation/results/a6-mmff94-analytic-m80-300-af7c0c44-2026-09-23.{jsonl,meta.json}`
- Repair rows and provenance:
  `validation/results/a6-mmff94-analytic-m80-300-repair-af7c0c44-2026-09-23.{jsonl,meta.json}`
- UFF rows and provenance:
  `validation/results/a6-uff-300-af7c0c44-2026-09-23.{jsonl,meta.json}`
- Machine summaries and common-scored rows:
  `validation/results/a6-*-vs-rdkit-2026.3.6-{summary,common-scored}-2026-09-23.*`
- Comparator:
  `validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.{jsonl,meta.json}`

Rebuild and rerun the candidate before interpreting the timing rows:

```bash
CARGO_BUILD_JOBS=1 cargo build --release -p chematic-3d \
  --example pipeline_v2_vs_rdkit_dump \
  --example pipeline_v2_vs_rdkit_common_scorer
SCHEMATIC_MMFF94_MAX_ITERATIONS=300 \
  python3 scripts/run_pipeline_hard_timeout.py \
  --binary target/release/examples/pipeline_v2_vs_rdkit_dump \
  --arm chematic_pipeline_v2_mmff94_strict \
  --timeout-s 25 \
  --output /tmp/chematic-a6-mmff94.jsonl \
  --metadata-output /tmp/chematic-a6-mmff94.meta.json
```
