# Bounded RDKit path enumeration

The isolated candidate adds ordered fixed-length path enumeration to Rust and
Python, corresponding to
[`Chem.FindAllPathsOfLengthN`](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.FindAllPathsOfLengthN).
It returns either bond-index paths or atom-index paths, with root, hydrogen
and shortest-path options. This complements the existing adjacency matrix,
shortest-path and radius-environment APIs. It does not add branched-subgraph
enumeration or chemical equivalence filtering.

## Interfaces

```python
import chematic

mol = chematic.from_smiles("C1CC1")
assert mol.rdkit_paths_of_length(2) == [[0, 1], [2, 1], [0, 2]]
assert mol.rdkit_paths_of_length(4, use_bonds=False) == [[0, 1, 2, 0]]
assert mol.rdkit_paths_of_length(3, only_shortest_paths=True) == []

paths = mol.rdkit_paths_of_length(
    2,
    use_bonds=True,
    use_hs=False,
    rooted_at_atom=0,
    only_shortest_paths=False,
    max_visits=1_000_000,
    max_frontier=100_000,
    max_results=100_000,
    max_path_indices=1_000_000,
)
```

Rust exports `RdkitPathParams`, `RdkitPathLimits`, `rdkit_paths_of_length`
and `rdkit_paths_of_length_with_limits` from `chematic_smiles`:

```rust
use chematic_smiles::{parse, rdkit_paths_of_length, RdkitPathParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mol = parse("C1CC1")?;
    let paths = rdkit_paths_of_length(&mol, 2, &RdkitPathParams::default())?;
    assert_eq!(paths, vec![vec![0, 1], vec![2, 1], vec![0, 2]]);
    Ok(())
}
```

`length` is an unsigned 32-bit integer. Python `rooted_at_atom=None` and Rust
`None` mean unrooted; a nonnegative root starts every returned path at that
atom. An out-of-range root returns an empty list, matching the pinned oracle.
Python negative roots are outside this API; use `None` instead of RDKit's
negative sentinel. Results are newly owned vectors/lists and do not mutate the
source molecule or its properties. WASM and the `rdkit_compat` wrapper have no
new adapter in this increment.

## Exact path semantics

| Option or boundary | Contract |
|---|---|
| `use_bonds=True` | Length counts bonds; output entries are bond indices |
| `use_bonds=False` | Length counts atoms; output entries are atom indices |
| Default options | Bond output, graph H excluded during expansion, unrooted, shortest filter disabled |
| Root | Constrains the first atom, rather than requiring the root anywhere in a path |
| Zero length | Empty result for either output mode |
| Atom length one | One singleton for each allowed start, including isolated atoms and graph Hs |
| Neighbors | Ascending sanitized atom index, independent of bond insertion order |
| Duplicate removal | Same set of bond indices means duplicate, even for atom output; retain the first ordered path |
| Last step | May revisit an earlier atom except the immediately preceding one, closing a cycle |
| Shortest option | Every retained prefix is a shortest unweighted path from its start |

Atom and bond indices refer to the existing sanitized RDKit model used by
`rdkit_adjacency_matrix`, `rdkit_shortest_path`, `rdkit_fragment_smiles` and
`rdkit_atom_environment`. They can differ from the raw chematic graph when
ordinary explicit hydrogens are removed. `rdkit_num_atoms()` gives the model's
atom-index range. Explicitly added H atoms and retained isotope H atoms follow
the established model contract.

Two upstream details are deliberate:

- With `only_shortest_paths=True`, RDKit constructs adjacency from its
  unfiltered distance matrix. Graph H atoms therefore participate even when
  `use_hs=False`. For `[2H]C`, ordinary bond-path enumeration excludes the
  H–C edge by default; shortest-mode enumeration includes it.
- A final cycle closure can revisit an interior atom, not just the starting
  atom. Starting from the tail atom 4 of `C1CCC1C`, atom length six returns
  `[4, 3, 0, 1, 2, 3]`. Bond-set deduplication retains one representation of
  that tail-plus-cycle path.

These are index-based paths, not unique chemical fragments. Symmetry-equivalent
paths with different bond sets remain distinct. Bond order, aromaticity and
dative direction do not supply weights to this operation. The shortest option
retains all enumerated geodesics after bond-set deduplication, not a single
representative shortest path.

## Bounded computation

The existing preparation guard checks at most 1024 input atoms and 16,384
input bonds before hydrogen removal. Unsupported inputs from the existing
sanitized model remain explicit errors. The search is iterative and uses no
depth-dependent recursion.

| Limit | Default | Counted work/storage |
|---|---|---|
| `max_visits` | 1,000,000 | Neighbor extensions examined, including rejected repeats and shortcuts |
| `max_frontier` | 100,000 | Prefix count in either breadth-first frontier, before deduplication |
| `max_results` | 100,000 | Unique returned paths |
| `max_path_indices` | 1,000,000 | Total indices in either frontier or in returned paths, checked separately |

Initial singleton prefixes count toward frontier limits. For example, an
unrooted `CCO` search needs room for three initial prefixes even when its final
result contains two bonds. The two frontiers, result storage and deduplication
keys have separate bounded storage; these settings are counts, not an exact
whole-process byte limit. Vector/hash-table growth uses fallible reservation.

Shortest-mode distances use unweighted all-pairs BFS over the bounded model,
equivalent to the default RDKit topological distance matrix. This preparation
is bounded by the input caps and is separate from `max_visits`, which counts
path expansion work.

Exceeding any limit returns `RdkitSmilesError::Unsupported` in Rust or
`ValueError` in Python, without a partial result. Zero length, an outside root
and a provably impossible length return empty without expanding a frontier.
Impossible lengths, including `UINT32_MAX`, terminate after bounded model
preparation. The oracle's loop does not stop when its frontier becomes empty,
so huge lengths are tested only in chematic; no upstream timing or overflow
compatibility is claimed.

## Pinned validation, 2026-10-15

The oracle is **RDKit 2026.09.1**, release `Release_2026_09_1`, resolved commit
`fece8caa860bdf6c9c82bdb9253f22a0cf87150c`. The port follows
[`Subgraphs.cpp`](https://github.com/rdkit/rdkit/blob/fece8caa860bdf6c9c82bdb9253f22a0cf87150c/Code/GraphMol/Subgraphs/Subgraphs.cpp),
specifically `extendPaths`, `pathFinderHelper`, `findAllPathsOfLengthsMtoN`
and `findAllPathsOfLengthN`. Retrieved source SHA-256:
`aa8295ce04042f809d93a29f8f6ca0beca2c656b4f3865270237131119a55597`.
Copyright (C) 2003–2026 Greg Landrum and other RDKit contributors;
the BSD-3-Clause attribution is recorded in `THIRD_PARTY_NOTICES.md`.

All checked-in graph recipes are original synthetic MIT OR Apache-2.0 data.
The manifest records provenance, redistribution, minimization, oracle settings,
expected behavior and normalized UTF-8 SMILES hashes. Eight explicit golden
queries additionally store their pinned expected vectors.

| Check | Result |
|---|---|
| Original graph recipes | 24 |
| Explicit golden queries | 8 |
| Graph/H/order variants | 134 |
| Exact ordered-vector comparisons | 107,904 |
| Paths / indices compared | 176,992 / 672,982 |
| Shortest-mode comparisons | 53,952 |
| Permuted-atom comparisons | 71,488 |
| Returned atom paths with final repeated atom | 1,758 |
| Unexpected mismatches | 0 |

Comparisons cover all bond/atom, H and shortest flags; every root, unrooted
search and two outside roots; H removal/addition; identity, reversed and seeded
atom permutations where distinct; lengths zero through `min(n+2, 7)` plus
`n+1` and `n+2`. Adjacency matrices and atom counts are checked before each
variant's path comparisons. Native-only tests cover `UINT32_MAX-1` and
`UINT32_MAX` without invoking upstream's large-length loop.

Local validation passed 386 Rust unit tests and 23 integration tests, with
one existing manual test ignored. Python 130 tests passed, including path,
atom-environment, graph-query, renumbering and MOPAC regressions. Formatting,
scoped Clippy with warnings denied, the unchanged refactor allowance budget,
binding inventory, product naming and strict MkDocs are recorded in the evidence.
No hosted CI or publication is claimed by these local results.

## Reproduction and delivery

```sh
cargo test -p chematic-smiles --lib --test rdkit_path_enumeration \
  --test rdkit_graph_queries --test rdkit_atom_environment
python -m pytest crates/chematic-py/tests/test_rdkit_path_enumeration.py
python scripts/validate_rdkit_paths.py
```

Build/install the candidate Python extension first. The oracle script requires
RDKit 2026.09.1 and checks that exact version. Fixtures live in
`validation/rdkit_issues/path-enumeration/fixtures.json`; the report is
`validation/results/rdkit-path-enumeration-2026-10-15.json`.

The source candidate is `/private/tmp/chematic-paths`, based on
`/private/tmp/chematic-mopac`. The incremental patch is
`validation/rdkit_issues/path-enumeration-2026-10-15.patch`.
Only this document, patch and evidence are copied to the primary checkout.
Unrelated primary-source work remains intact. The feature is uncommitted,
unpushed and unpublished.
