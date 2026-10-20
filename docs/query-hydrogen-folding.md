# Query hydrogen folding with tetrahedral parity

[RDKit #9675](https://github.com/rdkit/rdkit/issues/9675) is open as checked on
2026-10-14. RDKit 2026.09.1 reproduces its inversion: the query
`[C@@]([H])(C)(N)C(=O)O` matches `C[C@H](N)C(=O)O` before H folding when the
target has explicit H nodes, but its `mergeHs=True` form matches the opposite
configuration instead.

This local candidate starts from main commit `cac480fb`. It adds a conservative
native query transformation and mapped-neighbour tetrahedral matching. It is
not a package publication or a claim of general RDKit `MergeQueryHs` parity.

## Python

```python
import chematic

original = chematic.SmartsQuery("[C@@]([H])(C)(N)C(=O)O")
query = original.merge_hydrogens()
assert query.atom_count == original.atom_count - 1
assert query.matches(chematic.from_smiles("C[C@H](N)C(=O)O"), use_chirality=True)
assert not query.matches(chematic.from_smiles("C[C@@H](N)C(=O)O"), use_chirality=True)
```

The operation returns a new immutable `SmartsQuery`. Its `source` remains the
original input text, not serialized folded SMARTS; `repr` includes the folding
operation. `matches`, `matches_many`, `matches_batch`, and `find_matches` accept
`use_chirality=False` by default to preserve existing behavior. Set it to `True`
when tetrahedral configuration must constrain the match. These are compiled
query APIs; no new option is added to molecule string-query convenience APIs.

## Rust and scope

`chematic_smarts::merge_query_hydrogens(&QueryMolecule)` returns
`Result<QueryMolecule, QueryHydrogenError>`. It folds only plain, unmapped,
terminal `[H]` or `[#1]` nodes attached by a single or implicit bond.
For one removed H it adds `!H0` to the neighbour; for two it also adds `!H1`,
and so on. This is a minimum **total** H count, preserving existing predicates.
Targets are ordinary chemical graphs with monovalent H; the operation is not
an equivalence for arbitrary malformed multivalent-H graphs.

Isotopic, mapped, charged or additionally constrained H nodes survive, as do
H-H bonds and conditional/non-single bonds. Surviving atom-map values and
non-H bond predicates are preserved. Returned query indices are compacted;
match mappings refer to those new query indices, while target indices still
refer to the supplied molecule. The input query is unchanged and successful
folding is idempotent.

An affected chiral atom must have exactly one folded H and a complete
four-neighbour tetrahedral frame. Its `@`/`@@` sense is adjusted for the
virtual-H position. Partial/ambiguous frames and non-tetrahedral chirality
codes return an error. Recursive SMARTS are rejected by this transformation.
No H folding is performed inside recursive predicates.

Public graph consistency is checked before cloning. Limits are 1024 atoms,
4096 bonds, 16384 total atom/bond-expression nodes, and expression depth 32;
the folded output is checked against the same limits. Python maps folding
errors to `ValueError`.

The native VF2 path defers exact chirality until a complete mapping is known.
It compares neighbour permutations rather than raw tags, including the missing
implicit/explicit H slot. Three-valued candidate checks avoid pruning logical
NOT/OR expressions prematurely. Query ring openings reserve lexical slots,
including multiple open rings whose closures arrive in either order. Partial
unfolded chiral queries and targets without recorded stereo-neighbour order
retain their previous raw-tag fallback. The parity-preservation checks here use
SMILES-parsed targets with that order present. This does not
certify every SMARTS stereo class or the separate RDKit-parity matcher.

## Validation

- All `chematic-smarts` tests passed: 234 passed, 2 existing ignored tests.
  Ten new folding tests cover the issue, H positions, explicit/isotopic target H,
  minimum counts, retained constraints/maps, idempotence, invalid inputs/limits,
  and ring-opening/closure order.
- Python: 22 new cases and 68 existing module-function cases passed using the
  locally compiled source extension.
- `scripts/validation/validate_query_hydrogens.py` generates 48 acyclic patterns from all
  four-neighbour permutations and both tags, tested on four target spellings,
  plus six ring cases. All 198 folded chematic results match unmerged RDKit
  queries against explicit-H targets. RDKit 2026.09.1's folded queries differ
  from that reference in 102 of these cases. This is a fixed diagnostic set,
  not an estimate of real-world accuracy.
- Affected-crate all-target Clippy passed with warnings denied. Refactor budget
  passed without increasing existing allowances. Dated source/module hashes and
  oracle examples are in `validation/results/rdkit-query-hydrogens-2026-10-14.json`.

The implementation is in a separate working checkout so the original dirty
workspace is preserved. The changes are uncommitted and unreleased.
