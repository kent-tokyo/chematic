# MCS complete-ring candidate selection — 2026-10-22

This source change starts from `origin/main` at
`671d99cb9d2f61c0f3e4fe638ed89052beecad0c`. It corrects one bounded
MCS search result. It is not a published-package or broad RDKit comparison.

## Reproducer and oracle

Inputs are `C1CCCCC1.CCCC` and `C1CCCC1.CCCC`, in both orders. The first
components are different-sized saturated rings; both molecules also contain
the same four-carbon chain. With `complete_rings_only=true`, a partial ring
cannot be the answer, but the common chain is valid.

RDKit 2026.09.1 `rdFMCS.FindMCS` with both
`AtomCompareParameters.CompleteRingsOnly` and
`BondCompareParameters.CompleteRingsOnly` set to `True` returns an exhaustive
4-atom, 3-bond chain in both orders. The local macOS arm64 CPython 3.13 RDKit
wheel has SHA-256
`29b6176ace255989abde23d40d2187b5ed5effdb2f00dadaced24f181ab27046`.

| Input order | Previous CheMatic source | Corrected source | RDKit 2026.09.1 |
| --- | ---: | ---: | ---: |
| six-ring first | 0 atoms / 0 bonds | 4 atoms / 3 bonds | 4 atoms / 3 bonds |
| five-ring first | 0 atoms / 0 bonds | 4 atoms / 3 bonds | 4 atoms / 3 bonds |

The previous search selected its largest unconstrained mapping, then removed
partial-ring atoms at the end. That could erase the winning mapping without
recovering a smaller valid candidate already visited. The corrected search
keeps the best mapping that satisfies the complete-ring constraint while
exploring candidates. The pre-fix focused test failed at 0 versus 4 atoms;
the same test passes after the change.

## Verification and limits

- `cargo test --locked -p chematic-smarts`: all 232 unit tests and integration
  tests passed, including the 1,980-row pinned default MCS comparator profile.
- `cargo clippy --locked -p chematic-smarts --all-targets -- -D warnings`,
  `cargo fmt --all -- --check`, `git diff --check`, and the refactor-budget gate
  passed.
- The existing complete-ring and resource-limit tests passed. Both reproduced
  searches report `McsOutcome::Exhaustive`.

The correction can explore more candidates because a large invalid mapping no
longer sets the branch-and-bound threshold. No broad timing or memory result
was measured. The existing ring definition is checked on the first input's
ring set; full RDKit FMCS parity, stereochemistry, fingerprints, and independent
holdout accuracy remain open.
