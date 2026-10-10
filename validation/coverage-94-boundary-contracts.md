# Coverage boundary contracts

These tests extend the existing Rust coverage lane without removing files from
its measurement scope. The library changes only register test modules; molecular
algorithms are unchanged. Format tests check typed failures, resource boundaries,
metadata round trips, and original atom/bond indices. Mutation and numerical tests
check atom-order invariance, charge conservation, finite geometry, and analytical
relationships rather than snapshots of the implementation.

## Frozen external references

CI reads the committed fixtures and does not require Python or RDKit.
Regenerate the following fixtures with RDKit **2026.03.1**:

```sh
python scripts/generate_avalon_molblock_boundary_fixtures.py
python scripts/generate_canonical_direction_boundary_fixtures.py
python scripts/generate_potential_stereo_boundary_fixtures.py
python scripts/generate_serialization_geometry_boundary_fixtures.py
python scripts/generate_pdb_boundary_fixtures.py
python scripts/generate_mol2_boundary_fixtures.py
```

The generators check the installed RDKit version. Their references cover 14
Avalon V3000 blocks, 444 conjugated directional SMILES respellings, 29 potential
stereo cases, 2,025 SMARTS serializations, 303 deterministic depictions, 920 PDB
reader cases, and 58 MOL2 reader cases. PDB geometry cases enumerate all 24 square
planar, 120 trigonal bipyramidal, and 720 octahedral ligand permutations.

Regenerate the MMFF94 atom-type census with RDKit **2026.03.6**:

```sh
python scripts/mmff94_atom_type_oracle_tsv.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --output validation/rdkit-2026.03.6-mmff-type-census.tsv
```

The regular integration test compares all 9,796 reference-supported molecules
from the pinned 10,000-molecule corpus and requires zero atom-type differences.
The remaining rows retain the oracle's parse/unsupported status. The coverage
workflow now runs this formerly manual census and the existing CIP corpus gates.

## Explicit residuals

The depiction sample retains five layouts whose coordinates differ from RDKit
by more than 1e-9. The test requires finite, deterministic coordinates for every
case and the exact five-molecule residual set in
`crates/chematic-smiles/src/stereo_order_boundary_contract_tests.rs`.

The native approximate potential-stereo API differs on five of the 29 cases:
four dependent-ring cases and an arsenic center. The full RDKit atom-index references remain in
the fixture. The test in
`crates/chematic-chem/src/descriptor_boundary_contract_tests.rs` requires the
exact known residuals and exact agreement on the other cases. These gates do not
claim complete RDKit parity. Removing a residual requires updating its explicit
expectation after the corresponding algorithm change.

Coverage records code execution. It does not by itself establish chemical
correctness or compatibility with another toolkit.
