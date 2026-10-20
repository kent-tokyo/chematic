# Bounded native heavy-atom distance matrices

This isolated, unreleased candidate extends the allocation boundary motivated by
[RDKit #9647](https://github.com/rdkit/rdkit/issues/9647) to the native graph API.
The issue was OPEN when checked on 2026-10-15; its latest update was
2026-10-10T16:48:18Z. The upstream report describes 32-bit multiplication
overflow and heap corruption at 65,536 atoms. That C++ crash reproducer was
not executed here. chematic's safe Rust had a different remaining problem:
the Python native matrix method called an unbounded, infallibly allocating
native Rust function.

The earlier [RDKit-compatible matrix guards](rdkit-distance-matrix-limits.md)
already protected `rdkit_distance_matrix` and `rdkit_distance_matrix_3d`;
the [adjacency API](rdkit-graph-queries.md) was also bounded. This candidate
adds a checked native entry point and routes the existing Python
`Mol.topological_distance_matrix()` through it.

## API

```python
import chematic

mol = chematic.from_smiles("CCC")
assert mol.topological_distance_matrix() == [
    [0, 1, 2], [1, 0, 1], [2, 1, 0]
]
assert mol.topological_distance_matrix(
    max_atoms=3, max_bonds=2, max_entries=9
) == [[0, 1, 2], [1, 0, 1], [2, 1, 0]]

# Keyword budgets can tighten the hard limits.
# mol.topological_distance_matrix(max_entries=8) raises ValueError.
```

Rust exports from `chematic_chem`:

```rust
use chematic_chem::{
    NativeDistanceMatrixLimits,
    try_topological_distance_matrix,
    try_topological_distance_matrix_with_limits,
};

let mol = chematic_smiles::parse("CCC")?;
let distances = try_topological_distance_matrix(&mol)?;
assert_eq!(distances[0][2], 2);
let limits = NativeDistanceMatrixLimits {
    max_atoms: 3,
    max_bonds: 2,
    max_entries: 9,
};
assert_eq!(try_topological_distance_matrix_with_limits(&mol, &limits)?, distances);
# Ok::<(), Box<dyn std::error::Error>>(())
```

Errors use `NativeDistanceMatrixError`, with `ResourceLimit`,
`DimensionOverflow` and `Allocation` variants. Python maps these to
`ValueError`; negative or out-of-range Python integers fail unsigned conversion
with `OverflowError`. Budget arguments are keyword-only. Existing successful
no-argument calls retain their nested integer-list results.

## Limits and allocation order

| Resource | Hard cap | Counted before exclusion |
| --- | ---: | --- |
| Input atoms | 1024 | All atoms, including H and isotope H |
| Input bonds | 16384 | All bonds, including H bonds |
| Output entries | 1048576 | Square of non-H atom count |

Caller budgets are intersected with these hard caps; oversized budgets cannot
raise them. Zero is a valid tighter budget. An empty native graph satisfies
all zero budgets, while any input atom violates `max_atoms=0`.

Atom and bond counts are checked before graph traversal or buffer creation.
The heavy count then supplies a checked dimension product, checked against
the entry budget before allocating. Atom-position mapping, heavy indices,
outer rows, each distance row and the reusable BFS queue all use
`try_reserve_exact`. The reserved queue holds at most one entry per heavy atom.
Every failure returns an error and no partial matrix.

The dense u32 payload is at most **4 MiB**; this excludes row headers, mapping,
queue, source graph, allocator overhead and Python objects. Runtime is
O(H × (V + E)), with bounded input atoms and bonds. No sanitization, prepared
RDKit model, SMILES provenance or derived-cache population is required.

## Native semantics

- Rows/columns follow original atom order with every atomic-number-1 atom
  omitted, including isotopic H. Wildcard atoms (atomic number zero) remain.
- Distances count bonds in the induced non-H subgraph. Hydrogen cannot connect
  two retained atoms. Bond order, aromaticity, dative direction and coordinates
  do not weight these distances.
- The diagonal is zero. Disconnected entries are **4294967295 (`u32::MAX`)**.
  The Python stub previously printed the incorrect decimal 2147483647; its
  documentation is corrected without changing the value.
- Empty and H-only native graphs return an empty matrix within input budgets.
  The SMILES parser itself still rejects an empty string.

This is the existing native heavy-subgraph contract. RDKit's ordinary matrix
uses a different disconnected sentinel and can retain explicit H. The parity
oracle therefore explicitly deletes all H vertices without sanitizing and
normalizes RDKit's disconnected `1e8` to `u32::MAX`.

The legacy infallible Rust `topological_distance_matrix` remains available for
source compatibility. Native descriptors that call it, other descriptor
allocations, WASM bindings and parser input budgets are outside this change.
This candidate does not claim repository-wide dense-allocation protection.

## Validation

- Rust: **957 passed**, comprising 948 existing library tests and nine new
  integration tests; one pre-existing library test remains ignored.
- Python: **39 passed**, comprising 16 new tests, four existing RDKit matrix
  guard tests and 19 graph-query tests.
- RDKit **2026.09.1**: **1104 exact matrices / 48702 elements**, zero mismatches.
  The packet exhausts all 1024 simple labeled five-carbon graphs plus 18
  small boundary recipes and chains/cycles through 32 atoms. It includes
  disconnected graphs, rings, multiple/aromatic/dative bonds, explicit and
  isotope H, H bridges, wildcard atoms and stereochemical annotations.
- Native and Python reject 1025/46341/65535/65536/65537 input atoms before
  dense allocation, including attempts to raise the cap. Rust also tests exact
  1024-atom connected/disconnected matrices, exact 16384-bond acceptance,
  16385-bond rejection and tighter entry limits. Graphs built without SMILES
  provenance prove the native path needs no prepared model.

Native allocation-failure branches are implemented with fallible reservation;
allocator failure was not injected. The upstream oversized C++ crash was not
executed. Small-graph agreement does not establish general descriptor parity.

Synthetic recipes, original licensing, normalization and oracle settings are in
`validation/rdkit_issues/native-distance/fixtures.json`. Each generated record's
input and expected-matrix SHA-256 are saved in
`validation/results/native-distance-matrix-2026-10-15.json`.
The validator is `scripts/validation/validate_native_distance_matrix.py`.

```sh
cargo test -p chematic-chem --lib --test native_distance_matrix --offline --quiet
cargo clippy -p chematic-chem -p chematic-py --all-targets --offline --no-deps -- -D warnings
python3 -m pytest crates/chematic-py/tests/test_native_distance_matrix.py crates/chematic-py/tests/test_rdkit_matrix_limits.py crates/chematic-py/tests/test_rdkit_graph_queries.py -q
python3 scripts/validation/validate_native_distance_matrix.py
```

The candidate is `/private/tmp/chematic-native-matrix`, inherited from
`/private/tmp/chematic-pdb-selection`. Its incremental patch is
`validation/rdkit_issues/native-distance-matrix-2026-10-15.patch`, applicable
to that inherited source. Only documentation, patch and evidence are copied
back to the primary checkout; source changes remain isolated, uncommitted
and unpushed. The sparse Git baseline commit is a local comparison anchor,
not an integration commit.
