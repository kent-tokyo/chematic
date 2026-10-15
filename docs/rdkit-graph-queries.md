# RDKit-compatible graph queries

This source candidate adds two previously absent public operations:
`rdkit_adjacency_matrix` and `rdkit_shortest_path` in Rust, and methods of the
same names on Python `Mol`. The earlier candidate provided distance matrices;
its private layout/path helpers did not expose these operations.

The reference is [RDKit 2026.09.1's graph API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html)
and the adjacency/shortest-path routines in
[`Release_2026_09_1` Matrices.cpp](https://github.com/rdkit/rdkit/blob/fece8caa860bdf6c9c82bdb9253f22a0cf87150c/Code/GraphMol/Matrices.cpp).
The compatibility implementation retains the existing BSD-3-Clause notice in
`THIRD_PARTY_NOTICES.md`. No RDKit binary dependency is added.

## Python

```python
import chematic

mol = chematic.from_smiles("N->[Cu]")
assert mol.rdkit_adjacency_matrix() == [[0.0, 1.0], [1.0, 0.0]]
assert mol.rdkit_adjacency_matrix(use_bo=True) == [[0.0, 0.0], [1.0, 0.0]]
assert mol.rdkit_shortest_path(0, 1) == [0, 1]
assert mol.rdkit_shortest_path(1, 0) == [1, 0]
```

Adjacency returns independent `list[list[float]]` rows, including for the
binary profile. RDKit returns a NumPy array, with integer dtype for its binary
profile. An empty molecule gives `[]`; callers converting this to NumPy should
reshape it to `(0, 0)` when that shape is required.

Shortest paths are independent lists including both endpoints. A path counts
edges rather than bond order and treats dative bonds as undirected. Endpoints
in different fragments give `[]`. Equal or out-of-range endpoints raise an
error, matching RDKit's refusal of these calls; error classes/messages are not
claimed identical. Python extraction also rejects negative/oversized unsigned
indices. Supported computation/model errors are surfaced as `ValueError`.

## Rust

```rust
use chematic_smiles::{RdkitSmilesError, parse, rdkit_adjacency_matrix, rdkit_shortest_path};

fn graph_queries() -> Result<(), RdkitSmilesError> {
    let molecule = parse("CC.O").map_err(|e| RdkitSmilesError::Unsupported(e.to_string()))?;
    let adjacency = rdkit_adjacency_matrix(&molecule, false, 0)?;
    assert_eq!(adjacency[0][1], 1.0);
    assert!(rdkit_shortest_path(&molecule, 0, 2)?.is_empty());
    Ok(())
}
```

Both operations return `Result<_, RdkitSmilesError>`. Size, allocation,
unsupported-model and invalid-endpoint errors are explicit. The shared dense
matrix helper uses checked dimensions and fallible storage allocation. The
shortest-path search uses fallibly reserved linear buffers.

## Atom order and chemistry contract

These methods reuse `rdkit_mol_for_writing`, the sanitized model used by
`rdkit_distance_matrix`. Thus graph hydrogens removed by RDKit's SMILES
preparation disappear before indexing results. For example,
`[H]C([H])([H])[H]` has a one-row result. `rdkit_add_hydrogens()` explicitly
retains the generated hydrogens, with the same order used by other named
RDKit profiles. Indices are model indices and may differ from raw chematic
atom indices. The reader does not silently return a raw-graph matrix instead.

The unweighted adjacency is 1 for any supported bond and 0 otherwise.
The weighted profile uses each endpoint's valence contribution, including:

| Bond | Forward entry | Reverse entry |
|---|---|---|
| Single / double / triple / quadruple | 1 / 2 / 3 / 4 | same |
| Aromatic-typed bond | 1.5 | 1.5 |
| Dative donor → acceptor | 0 | 1 |

Shortest paths search neighbors in the model's RDKit bond-index order.
Different equal-length paths can be chosen when the call direction is
reversed. For a four-membered ring `C1CCC1`, path `(1, 3)` is `[1, 0, 3]`,
while path `(3, 1)` is `[3, 2, 1]`. This ordering is checked against the oracle
after atom permutations as well as for original inputs.

The operations do not mutate the molecule, its properties or a caller-owned
result from an earlier call. Binary and weighted calls allocate fresh results;
there is no shared matrix cache and no `force`/`prefix` parameter.

## Bounds and unsupported branches

- At most 1024 input atoms and 16,384 input bonds, checked before model
  preparation or hydrogen removal. The prepared matrix is also checked
  against the existing 1,048,576-entry dense cap.
- Self bonds, zero/query bonds and inputs unsupported by the existing
  sanitized model fail explicitly. The existing model's restrictions on
  non-SMILES stereo annotations also apply. These APIs do not add generic
  query-graph, unsanitized-molecule or arbitrary CTAB support.
- Adjacency supports `empty_value=0` only. Nonzero values return an explicit
  error, including on an empty graph. The pinned upstream implementation
  initializes its double buffer with `memset(emptyVal)`. A measured weighted
  `emptyVal=1` produces `7.748604185489348e-304` on the diagonal rather than
  numeric 1; other byte values can yield NaNs. This behavior is not reproduced.
- RDKit's cache can also retain the previous fill when `emptyVal` changes
  without `force=True`; the oracle gate forces recalculation for comparison.
- This delivery exposes native Rust and Python methods. It does not add
  WASM exports, a NumPy convenience method, or `rdkit_compat` facade functions.

## Validation and regeneration

`scripts/validate_rdkit_graph_queries.py` pins RDKit **2026.09.1** and verifies
the checked-in manifest hash. It uses 27 original minimal structures under
`MIT OR Apache-2.0`; no external molecular dataset is redistributed. Each
fixture records provenance, license, minimization, oracle/settings, expected
results and a normalized UTF-8 payload hash in
`validation/rdkit_issues/graph-queries/fixtures.json`.

The gate compares both binary and weighted matrices, repeats binary calls to
check profile independence, and compares every ordered pair of distinct atoms.
Identity, reversed and seeded shuffled atom orders are used where distinct.
It covers explicit-H removal and addition, dative direction, aromatic/Kekulé
input, cleanup, isotope and map metadata, fragments, salts, bridges, cubane,
and quadruple bonds. Path lengths are also checked against RDKit distances.
Molecule serialization and a sentinel property must remain unchanged.

| Local oracle check, 2026-10-15 | Result |
|---|---|
| Original fixtures / graph variants | 27 / 139 |
| Matrix calls / exact entries | 417 / 30,033 |
| Exact ordered shortest paths | 9,072 |
| Disconnected paths among those calls | 1,154 |
| Direction-dependent tied paths | 798 |
| Unexpected mismatches | 0 |

Rust validation additionally checks the exact 1024-atom boundary, a graph
exceeding 65,535 atoms, an edge-limit violation before invalid-valence
sanitization, pre-removal H limits and unsupported branches. The relevant
library/integration suites pass 398 tests, with one pre-existing manual test
ignored. Python graph-query tests pass 19 cases; with existing matrix,
renumbering and EXYZ regressions the scoped Python run passes 94 cases.

```sh
python scripts/validate_rdkit_graph_queries.py
cargo test -p chematic-smiles --lib --test rdkit_graph_queries --test rdkit_matrix_limits
```

The Python command requires the built candidate extension on `PYTHONPATH` and
the pinned RDKit oracle. Measurements describe this local source candidate;
they are not published-package evidence, a performance result, or general
RDKit parity. The candidate inherits the previous EXYZ and Layered feature
work. Its source is isolated from the primary dirty checkout; a portable patch
and evidence report are delivered without committing or pushing.
