# A6 conformer-tail development record

Date: 2026-10-13

This record covers source-candidate development against RDKit 2026.09.1. It
does not establish universal or published-package superiority.

## Fixed comparison

- Oracle: RDKit 2026.09.1
- Host: macOS arm64
- Work: 20 ETKDGv3 candidates per molecule, seed 42, MMFF94, 200 initial
  minimizer iterations
- Quality: best retained MMFF94 energy, with convergence, stereo and failure
  accounting kept separate
- Timing: alternating engine order; the formal gate requires repeated paired
  measurements rather than a single run

The original 265-row A6 development set and the former 100-row v2 holdout have
both been inspected during diagnosis. They are development evidence from this
point onward and are not reused as untouched final evaluation data.

## Candidate changes

1. Degenerate torsion or inversion terms no longer stop evaluation of later
   valid gradient terms.
2. ETKDG relevant rings use RingDecomposerLib traversal order.
3. MMFF optimization edits the conformer coordinate storage in place.
4. A separately named rescue profile retries low-energy non-converged
   candidates only when no initial candidate converged. Existing behavior is
   unchanged unless callers opt in.

## Development results

| Lane | Within 1 kcal/mol | At or below RDKit | Maximum delta | Timing |
|---|---:|---:|---:|---:|
| 265-row development, normal profile | 257/265 | 167/265 | +3.134 | one run previously favored CheMatic; not an accepted speed gate |
| used 100-row v2 set, normal profile | 94/100 | 59/100 | +6.253 | development-only |
| used 100-row v2 set, selective rescue (1 × 200) | 97/100 | 75/100 | +5.482 | CheMatic median 25.307 s; RDKit 24.571 s, so CheMatic was about 3% slower |

The rescue profile improves the used v2 quality tail but does not meet the
simultaneous speed target. It remains optional and is not evidence that the
default beats RDKit.

## Frozen v3 evaluation

The candidate was frozen at `0b7a197db7595e079911b147c96f3eec53a85480`.
The outcome-blind generator then selected v3 while excluding Tier A, Tier B,
holdout v1 and holdout v2. The manifest contains 50 declared-stereo and 50
general molecules, has no overlap with v1/v2, and has SHA-256
`a2962ce0f33fc10523cf33e642714772a708d19468814a7bd6bb551026050608`.

The exact candidate wheel had SHA-256
`ade612cda0ef976b95b8847caa1488fd83c9338fa517a85fa149718c0de3e04d`.
The two profiles below were declared before either result was inspected:

- **performance profile:** 20 candidates, 200 MMFF iterations, no retry;
- **quality-rescue profile:** the same work, plus one 200-iteration retry only
  when no initial candidate converged.

| v3 profile | Within 1 kcal/mol | Within 2 kcal/mol | At or below RDKit | Maximum delta | Converged | Median seconds, CheMatic / RDKit |
|---|---:|---:|---:|---:|---:|---:|
| performance | 94/100 | 97/100 | 61/100 | +3.521 | 887 / 891 | 29.845 / 29.764 |
| quality rescue | 95/100 | 99/100 | 73/100 | +3.365 | 901 / 891 | 30.860 / 30.214 |

The performance profile is effectively tied on this three-repeat host result
but is 0.27% slower by the median. The rescue profile improves convergence and
the quality tail but is 2.1% slower. Neither profile establishes simultaneous
quality-and-speed superiority. Do not tune from v3 or use it as a sealed set
again. A broad claim additionally requires a fresh evaluation set, a second
host, peak memory, 21 alternating blocks, and exact public artifacts.
