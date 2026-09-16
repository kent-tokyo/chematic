# Roadmap open-work disposition

Updated 2026-09-16. This file records dependency classes and issue boundaries;
priority order belongs only in [`ROADMAP.md`](../ROADMAP.md), and numerical exits
belong in the [accuracy plan](rdkit-accuracy-plan.md).

## Dependency classes

| Class | Meaning | Current examples |
|---|---|---|
| `local-open` | Can be implemented and tested from the repository | canonical E/Z residuals, binding regressions, documentation synchronization |
| `toolchain-open` | Needs a pinned external runtime or hosted environment | RDKit/Indigo readers, Linux namespace/memory enforcement, browser engines |
| `data-sealed` | Input exists but may not be inspected before the frozen run | 8,000-row sealed holdout |
| `external-open` | Requires another person or external decision | independent chemical adjudication, external security review |
| `historical` | Retained evidence; not active work | abandoned additional 1.10× speed stretch and older version snapshots |

## Open issue boundaries

- #149 and #503: coupled canonical E/Z carrier systems remain fail-closed; the
  accepted exit is invariant semantic identity without index-based heuristics.
- #372: optimize Boc/tBu symmetry only after exact-output and budget-exhaustion
  gates remain unchanged.
- #337 and #227: MMFF94 typing/coverage need same-input RDKit availability and
  numerical validation; partial type coverage is not parity.
- #185: UFF convergence requires geometry and energy soundness, not only a
  finite return value.
- #70: process-level benchmark calibration remains necessary before performance
  verdicts are trusted.

## Current external dependencies

- A5 needs non-maintainer review and independent gold labels.
- Parser-security comparison lanes need pinned upstream releases and legally
  reproducible minimized cases.
- Browser memory claims need a defined measurement API and equivalent process
  boundary; summed process-tree RSS is diagnostic only.
- Published release completion requires registry and site observation after
  upload; local packaging cannot substitute for it.

## Evidence policy

Local, historical, and externally verified records stay distinct. A narrower
passing fixture does not close a broader package. Every open item must point to
an implementation/test change or an immutable result; prose alone does not
advance status.

Detailed pre-consolidation status is preserved in
[`docs/archive/roadmap-through-2026-09-13.md`](archive/roadmap-through-2026-09-13.md).
