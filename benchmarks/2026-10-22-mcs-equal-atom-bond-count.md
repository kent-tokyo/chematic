# MCS equal-atom bond-count correction

This source change starts from `origin/main` at
`9b085d59965a5bb918a0b38d9b66bfd59e3fb7cb`. It fixes a search
completeness error in `chematic-smarts`; it is not a published-package or
toolkit-wide accuracy result.

## Reproducer

For `CCC.C1CC1` and `C1CC1`, the first chain candidate contains three atoms
and two common bonds. The later cyclopropane candidate has the same atom count
and three bonds. The search's atom-count upper bound previously discarded the
later candidate when its bound equaled the current best atom count.

| Input order | Previous source | Corrected source | RDKit 2026.09.1 `FindMCS` |
| --- | ---: | ---: | ---: |
| `CCC.C1CC1`, `C1CC1` | 3 atoms / 2 bonds | **3 atoms / 3 bonds** | 3 atoms / 3 bonds |
| `C1CC1`, `CCC.C1CC1` | not separately measured | **3 atoms / 3 bonds** | 3 atoms / 3 bonds |

The corrected search only prunes a branch when its atom-count upper bound is
strictly below the best atom count. Both the seed and recursive bound use this
condition. The same-size, higher-bond candidate can therefore reach the
existing `maximize_bonds` tie breaker.

## Verification and limits

- The chain-first reproducer failed before the correction. Both input orders
  pass after the correction.
- All 1,980 atom/bond-count rows in the pinned RDKit 2026.09.1 MCS comparator
  profile still match its stated native-policy expectations; all searches
  finish exhaustively. Its four explicit-hydrogen policy differences from raw
  RDKit remain separate.
- The 47 focused MCS unit tests pass. The complete `chematic-smarts` test
  suite, Clippy, format and refactor-budget checks pass.
- The expanded equal-size search space can cost more work on larger symmetric
  molecules. No broad performance comparison was made.

This tests a bounded MCS search correction. It does not establish parity for
generalized SMARTS, all FMCS options, stereochemistry, or fingerprints.
