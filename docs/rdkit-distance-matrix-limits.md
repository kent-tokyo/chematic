# RDKit-compatible distance matrix limits

The unreleased source adds an allocation boundary to
`rdkit_distance_matrix` and `rdkit_distance_matrix_3d`, motivated by
[RDKit #9647](https://github.com/rdkit/rdkit/issues/9647).
The issue reports 32-bit dimension overflow and heap corruption in upstream
matrix APIs at 65,536 atoms. chematic's safe Rust implementation does not
have that C++ out-of-bounds write, but previously attempted unchecked dense
allocations and used infallible buffer copies.

## Contract

Both existing APIs now accept at most **1024 input atoms** and at most
1024 atoms in their prepared RDKit-compatible model. A matrix can therefore
contain at most **1,048,576 entries**.

- The input check runs before model preparation, sanitization, or H removal.
- The prepared-model check runs before allocating the flat matrix.
- Dimension multiplication uses `checked_mul`.
- Flat matrices, both Floyd-Warshall work buffers, and nested output rows
  use fallible capacity reservation. Reservation errors are returned.
- Existing return types are unchanged: Rust returns
  `RdkitSmilesError::Unsupported` with a dimension/limit/allocation message;
  Python raises `ValueError`.

The topological path has three dense f64 buffers during Floyd-Warshall,
whose combined payload is at most 24 MiB. This is a bound on those buffers,
not total process memory, molecule metadata, or Python object overhead.

```python
import chematic

mol = chematic.from_smiles("C=CO")
assert mol.rdkit_distance_matrix(use_bo=True)[0] == [0.0, 0.5, 1.5]
```

## Compatibility and scope

Values and atom ordering within the supported size remain unchanged,
including bond-order weights, atomic-weight diagonals, `1e8` between
fragments, and the port's explicit-H removal behavior.

The new input limit is conservative: explicit H nodes count towards it
even when model preparation would remove them. For example, 513 separate
`[H]C` fragments contain 1026 native atoms and are rejected before preparation,
although their H-removed matrix would be smaller. Larger inputs that were
previously attempted now receive a documented error.

This change covers the two named RDKit-compatible distance APIs. It does
not add an adjacency-matrix API or change native topological descriptors,
distance-geometry bounds, coordinate validation, or numeric edge-case policy.

## Source validation

Validation passed:

- 379 library tests and four public-API integration tests for
  `chematic-smiles` (383 passed, zero failed, one pre-existing ignored test).
- 11 Python tests: four new size/error checks and seven existing RDKit writer
  and compatibility tests.
- Clippy with warnings denied for all `chematic-smiles` and `chematic-py`
  targets, plus formatting checks on changed Rust files.

New tests cover the first rejected size (1025), upstream's reported
46,341/65,535/65,536/65,537 boundaries, `usize::MAX` dimension overflow,
the exact accepted Cartesian cap of 1024, and H-removal/preparation order.
The large public Rust tests use graphs without SMILES provenance, so the
specific resource-limit error also proves rejection precedes port preparation.
These tests do not allocate the oversized dense matrices.

A live check compared 30 small matrix configurations against RDKit
**2026.09.1**, covering double bonds, aromatic bonds, fragments, removed
explicit H and added H, bond-order weighting, atomic-weight diagonals and
Cartesian distances. Maximum absolute difference was **0.0** on this packet.
Both Python APIs also rejected 1025- and 65,536-atom inputs with `ValueError`.

The upstream large-input crash reproducer was not executed. The upstream
failure is attributed to its issue report; the local evidence concerns
chematic's guards and small-input numeric behavior. The source-lane manifest,
extension hash, comparisons and guard messages are recorded in
[`rdkit-matrix-limits-2026-10-14.json`](https://github.com/kent-tokyo/chematic/blob/main/validation/results/rdkit-matrix-limits-2026-10-14.json).

```sh
CARGO_INCREMENTAL=0 cargo test -p chematic-smiles --lib --offline --quiet
CARGO_INCREMENTAL=0 cargo test -p chematic-smiles --test rdkit_matrix_limits --offline --quiet
python -m pytest crates/chematic-py/tests/test_rdkit_matrix_limits.py crates/chematic-py/tests/test_rdkit_writers.py -q
CARGO_INCREMENTAL=0 cargo clippy -p chematic-smiles -p chematic-py --all-targets --offline --no-deps -- -D warnings
```

Python checks used a local debug extension without native InChI support.
This is local source validation, not a package publication or upstream change.
