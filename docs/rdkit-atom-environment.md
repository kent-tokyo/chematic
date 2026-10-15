# RDKit atom-radius environments

This candidate adds Rust `rdkit_atom_environment` and Python
`Mol.rdkit_atom_environment()`, corresponding to
[RDKit's FindAtomEnvironmentOfRadiusN](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.FindAtomEnvironmentOfRadiusN).
The previous candidate had a private, diagnostics-only atom-ball BFS and
Morgan expansion traces; neither exposed ordered radius bond environments or
the complete `useHs`/`enforceSize` contract.

The reference is RDKit **2026.09.1**, tag `Release_2026_09_1`, commit
`fece8caa860bdf6c9c82bdb9253f22a0cf87150c`, specifically
`Code/GraphMol/Subgraphs/Subgraphs.cpp::findAtomEnvironmentOfRadiusN`.
The implementation retains its BSD-3-Clause attribution in
`THIRD_PARTY_NOTICES.md` and adds no RDKit binary dependency.

## Python

```python
import chematic

mol = chematic.from_smiles("CCO")
environment = mol.rdkit_atom_environment(2, 0)
assert environment == {
    "bond_indices": [0, 1],
    "atom_distances": {0: 0, 1: 1, 2: 2},
}
assert mol.rdkit_atom_environment(3, 0) == {
    "bond_indices": [], "atom_distances": {},
}
assert mol.rdkit_atom_environment(3, 0, enforce_size=False) == environment
```

The signature is:

```python
mol.rdkit_atom_environment(
    radius, rooted_at_atom, *, use_hs=False, enforce_size=True,
    max_visits=1_000_000, max_frontier=100_000,
)
```

The result contains an ordered `bond_indices` list and an `atom_distances`
dictionary mapping atom indices to minimum bond-hop radius. Unlike RDKit's
optional `atomMap` output argument, this API returns an independent dictionary
and accepts no mutable caller map. Rust uses a `BTreeMap`, so output dictionary
iteration is ascending by atom index; RDKit dictionary iteration order is not
part of the compatibility claim. Mapping values and bond-list order are exact.

## Rust

The crate root exports:

- `RdkitAtomEnvironment { bond_indices, atom_distances }`;
- `RdkitAtomEnvironmentLimits { max_visits, max_frontier }`;
- `rdkit_atom_environment(mol, radius, rooted_at_atom, use_hs, enforce_size)`;
- `rdkit_atom_environment_with_limits(..., &limits)`.

Both functions return `Result<RdkitAtomEnvironment, RdkitSmilesError>`.
The default function uses the same budgets as Python. Invalid roots,
unsupported model input and exhausted search budgets return explicit errors.
No partial environment is returned on failure. Python computation errors are
`ValueError`; unsigned argument extraction can produce a different exception
for negative/out-of-range integers.

## Exact behavior and boundaries

These functions reuse the existing sanitized RDKit model and its ordered
neighbor lists. Atom and bond IDs refer to that model, not necessarily the raw
chematic graph. Graph-H removal and explicit `rdkit_add_hydrogens()` follow
the same preparation policy as the existing graph-query/distance APIs.
For `[H]C([H])([H])[H]`, the prepared model contains one carbon atom.
After `rdkit_add_hydrogens()`, its explicit hydrogen atoms remain available.

Important boundaries covered by the implementation and oracle:

- Radius zero returns an empty bond list and `{root: 0}`.
- Positive radii expand ordered bond frontiers. Each bond is returned once,
  in discovery order. Atom distances retain their minimum value.
- With `use_hs=False`, hydrogen **neighbors** are skipped. A hydrogen root
  can still reach a heavy neighbor. Thus `[2H]C` rooted at the H includes
  its bond at radius 1, while the C root excludes the H unless `use_hs=True`.
  Hydrogen isotope labels do not bypass this neighbor filter.
- `enforce_size=True` empties both outputs when the requested frontier depth
  cannot be reached. `False` returns everything accumulated before exhaustion.
- Pending duplicate bonds are retained in the frontier even though already
  returned bonds are never emitted twice. This affects the final depth check
  in cycles. For `C1CCC1`, root 0 gives `[0, 3, 1, 2]` at both radius 2 and
  radius 3 with size enforcement, while radius 4 empties both outputs. A
  plain atom-distance ball or deduplicated queue would change this behavior.
- Bonds are traversed as graph connections, including both directions of
  dative bonds. This differs from endpoint-weighted adjacency contributions.
- Out-of-range prepared-model roots fail explicitly. Raw input indices must
  be remapped when H removal changes numbering: a raw index that remains in
  range addresses that prepared-model atom, not its former input identity.
  An empty molecule has no valid root. Query/zero bonds and other unsupported
  sanitized-model inputs are refused, as in existing named graph APIs.

Input caps are 1024 atoms and 16,384 bonds, checked before model preparation
or H removal. The visit budget counts every popped frontier entry, including
duplicates; the frontier budget caps pending entries in each individual layer.
Frontier and bond-list allocation is fallible. These budgets cover expansion;
model preparation remains subject to the input caps.

Radius is an unsigned 32-bit value. There is no radius-sized buffer or
unconditional radius-sized loop: traversal stops when the finite frontier is
empty or a budget is exceeded. Tests include `UINT32_MAX` on a 1024-atom chain.
A radius-zero call needs no frontier budget and succeeds with both budgets 0
for a supported, valid root.

The operation preserves molecule properties and serialization. Results from
separate calls can be edited independently. This delivery adds native Rust
and Python methods; it does not extend WASM or the `rdkit_compat` facade.

## Validation

`scripts/validate_rdkit_atom_environment.py` pins RDKit 2026.09.1 and verifies
the manifest and individual normalized UTF-8 SMILES hashes. Its 30 original
minimal project structures are recorded with provenance, license,
redistribution, minimization, oracle settings and expectations in
`validation/rdkit_issues/atom-environment/fixtures.json`. They are covered by
`MIT OR Apache-2.0`; no external molecule dataset is redistributed.

The gate includes every model root, both hydrogen and size flags, radii from
0 through `min(9, n+2)` plus `UINT32_MAX`, H removal/addition, and identity,
reversed and seeded shuffled atom orders where distinct. It compares against
an RDKit map seeded with unrelated values to verify the upstream map-clearing
contract. Adjacency and atom counts first establish the prepared indexed
graph, and serialization plus a sentinel property must remain unchanged.

| Local check, 2026-10-15 | Result |
|---|---|
| Original fixtures / graph variants | 30 / 162 |
| Exact environment comparisons | 43,236 |
| Ordered bond IDs compared | 165,177 |
| Atom distances compared | 176,798 |
| Radius-zero cases | 4,292 |
| Cases with changed atom order | 28,656 |
| Unexpected mismatches | 0 |
| Rust library + environment/graph/matrix integration tests | 406 passed; 1 pre-existing manual test ignored |
| Python environment/graph/matrix/renumber/FHI-aims tests | 106 passed |

Additional Rust/Python contracts cover the exact atom cap, H removal after the
input cap, a graph above 65,535 atoms, invalid roots, query-bond refusal,
visit/frontier failures, cyclic duplicate visits, detached outputs and zero
budgets for radius zero. The prepared-model and ordinary distance/graph APIs
share the existing code rather than a new chemistry interpretation.

```sh
python scripts/validate_rdkit_atom_environment.py
cargo test -p chematic-smiles --lib --test rdkit_atom_environment \
  --test rdkit_graph_queries --test rdkit_matrix_limits
```

The Python gate requires the built candidate extension and the pinned oracle.
Measurements describe this local source candidate, not a published package,
performance comparison, or general RDKit parity. The isolated candidate
inherits the preceding FHI-aims and graph-query features; its portable patch
and evidence report are delivered without committing or pushing.
