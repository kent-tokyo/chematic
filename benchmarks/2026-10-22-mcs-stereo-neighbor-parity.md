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
- [Machine-readable counts](2026-10-22-mcs-stereo-neighbor-parity.json) are
  consumed by the Rust regression test. To rerun the oracle, parse each JSON
  pair with `Chem.MolFromSmiles`, set `MCSParameters().AtomCompareParameters.MatchChiralTag=True`,
  and call `rdFMCS.FindMCS` with the two molecules. The Rust check is
  `cargo test --locked -p chematic-smarts chiral_mcs`.

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

## Scope and adjudicated oracle difference

The fix compares the ligand permutation under the actual MCS atom mapping.
Partially mapped centers remain eligible; allene and non-tetrahedral stereo
retain their prior policy. This result does not establish FMCS-wide parity,
fingerprint accuracy, or overall superiority to RDKit.

For the mixed three-input set `F[C@H](Cl)Br`, `Cl[C@@H](F)Br`, and
`Cl[C@H](F)Br`, RDKit 2026.09.1 is input-order dependent: the six
permutations return either 2/1 or 3/2 depending on which molecule is first.
CheMatic returns 3/2 for all six permutations. This row is therefore classified
as `oracle_order_unstable`, not tuned to one RDKit ordering. A Rust regression
test fixes CheMatic's permutation-invariant result.
