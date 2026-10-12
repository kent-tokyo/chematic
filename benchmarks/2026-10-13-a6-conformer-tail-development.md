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

## Frozen next evaluation

Before seeing any CheMatic result, generate a new v3 manifest from the source
corpus while excluding Tier A, Tier B, holdout v1 and holdout v2. Freeze the
candidate commit and then evaluate these two predeclared profiles exactly once:

- **performance profile:** 20 candidates, 200 MMFF iterations, no retry;
- **quality-rescue profile:** the same work, plus one 200-iteration retry only
  when no initial candidate converged.

Report both profiles independently. Do not choose a default or tune thresholds
from v3. A broad speed-and-quality claim additionally requires a second host,
peak memory, 21 alternating blocks, and exact public artifacts.
