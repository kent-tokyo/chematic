# MCS tetrahedral neighbor parity, 2026-10-22

## Reproduction

- Base: `978669e4264683bb73f5a7a3389a945b5fdf2b87`.
- Oracle: RDKit `2026.09.1`, `rdFMCS.FindMCS`, default `MCSParameters`
  with `AtomCompareParameters.MatchChiralTag=True`.
- CheMatic: `find_mcs_with_config_checked`, default `McsConfig` with
  `match_chiral_tag=true`; both molecule orders were checked. All reported
  searches were exhaustive.
- Measurement: query atom and bond counts. This is a correctness comparison,
  not a runtime benchmark.

| Inputs | RDKit atoms/bonds | CheMatic before | CheMatic after |
| --- | ---: | ---: | ---: |
| `F[C@H](Cl)Br` / `Cl[C@@H](F)Br` (same configuration) | 4/3 | 1/0 | 4/3 |
| `F[C@H](Cl)Br` / `Cl[C@H](F)Br` (opposite configuration) | 3/2 | not recorded | 3/2 |
| `F[C@H](Cl)Br` / `Br[C@@H](Cl)F` (same configuration) | 4/3 | not recorded | 4/3 |
| `F[C@](Cl)(Br)I` / `Cl[C@@](F)(Br)I` (same configuration) | 5/4 | not recorded | 5/4 |
| `F[C@](Cl)(Br)I` / `Cl[C@](F)(Br)I` (opposite configuration) | 4/3 | not recorded | 4/3 |
| `N[C@@H](C)C(=O)O` / `N[C@H](C)C(=O)O` (opposite configuration) | 5/4 | not recorded | 5/4 |
| `N[C@@H](C)C(=O)O` / `C[C@H](N)C(=O)O` (same configuration) | 6/5 | not recorded | 6/5 |

The first six pairs above are checked in both input orders: 12/12 count
comparisons match RDKit. The seventh alanine pair is also checked in both
orders: 14/14 total. Three equivalent spellings of the first molecule as a
three-input MCS also return 4/3 with both implementations.

## Scope and remaining difference

The fix compares the ligand permutation under the actual MCS atom mapping.
Partially mapped centers remain eligible; allene and non-tetrahedral stereo
retain their prior policy. This result does not establish FMCS-wide parity,
fingerprint accuracy, or overall superiority to RDKit.

For the three-input set `F[C@H](Cl)Br`, `Cl[C@@H](F)Br`, and
`Cl[C@H](F)Br`, RDKit returns 2/1 while CheMatic returns 3/2. This difference
is unresolved and is not counted as a pass. Future work should determine
whether RDKit's multi-input stereo rule intentionally excludes a partial
center before changing either implementation.
