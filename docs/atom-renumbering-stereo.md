# Atom permutations with preserved stereo

This local candidate adds the missing atom-permutation workflow to chematic,
using [RDKit #9629](https://github.com/rdkit/rdkit/issues/9629) as a ring/stereo
invariance diagnostic. The native baseline already counted one stereocenter
for `C1CCN2CCCC2C1`; it had no public atom-renumbering API. This addition makes
that count stable through both identity and nontrivial atom permutations.

## Use

```python
import chematic
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import rdMolDescriptors

mol = Chem.MolFromSmiles("C1CCN2CCCC2C1")
order = list(reversed(range(mol.GetNumAtoms())))
copy = Chem.RenumberAtoms(mol, order)
assert rdMolDescriptors.CalcNumAtomStereoCenters(mol) == 1
assert rdMolDescriptors.CalcNumAtomStereoCenters(copy) == 1

native = chematic.from_smiles("[H][C@](F)(Cl)Br")
assert native.atom_count == 5
copy = native.renumber_atoms(list(reversed(range(native.atom_count))))
assert copy.smiles == native.smiles
```

Rust provides `Molecule::renumber_atoms(&[AtomIdx]) -> Result<Molecule,
RenumberError>`, with `RenumberError` re-exported by `chematic_core`.
The order convention is `new_order[new_index] = old_index`, matching
[RDKit's RenumberAtoms documentation](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RenumberAtoms).

Python's native `Mol.atom_count` reports the number of stored graph atoms,
including every explicit H. The separate `rdkit_num_atoms()` descriptor
operates on the RDKit-prepared view, which can remove ordinary explicit H;
use the stored count to build atom permutations. The compatibility wrapper
now makes `GetNumAtoms()` return the stored count by default, with
`GetNumAtoms(onlyHeavy=True)` available for the heavy-atom count.

## Preserved meaning

- Atom fields, isotope/H counts, atom maps, nonchemical tags and R labels
  follow their source atoms.
- Bond indices, bond orders, endpoint orientation and each atom's ordered
  adjacency slots stay fixed relative to their original graph entities.
  Dative donor/acceptor orientation is preserved.
- Remap atom references in tetrahedral/axial/coordination stereo frames,
  enhanced stereo groups and aromatic direction anchors. Preserve the virtual
  implicit-H sentinel and each frame's ordering. Because neighbour slots are
  unchanged, the atom chirality tag retains its meaning without a parity flip.
- Retain enhanced bond-group membership, directional-bond stashes and SMILES
  ring-closure metadata; their bond indices do not change.
- Start every derived cache empty, even for the identity permutation.
  Recompute ring and descriptor data from the copied graph. The source graph
  and its already-computed caches stay intact.
- Python copies molecule properties with their types and name, independently
  of the source property store. The compatibility wrapper also starts with
  empty atom/bond/ring wrapper caches.

The permutation must include every source atom exactly once. Wrong length,
duplicates and out-of-range atom indices return typed Rust errors or Python
`ValueError`. Python also rejects bool, noninteger, negative and oversized
indices. Invalid references in caller-set stereo/bond metadata return a typed
error rather than indexing out of bounds. An empty graph accepts an empty order.

This operation changes indices without sanitizing or assigning new chemistry.
The core graph does not own conformers or file-record arrays. Reorder external
coordinates, PQR rows and other per-atom arrays separately using the same
new-to-old order. Custom property payloads are copied verbatim; any atom
indices encoded inside those payloads remain the caller's responsibility.

## Evidence

With RDKit 2026.09.1 and legacy stereo perception enabled, the issue example
has count 1 after parsing, count 2 after identity `RenumberAtoms`, and count 1
when `GetSSSR` is called on a separately renumbered reference before counting.
The chematic copy remains at 1.

The independent validation recipe covers 15 small structures: bicyclic amines,
both tetrahedral configurations, explicit/implicit/isotopic H, atom maps,
E/Z, allene, square-planar, trigonal-bipyramidal and octahedral stereo. It uses
all permutations for at most five atoms; larger graphs use identity, reverse,
cyclic shifts and 16 seeded random shuffles, with duplicates removed.

Across 838 permutations, every native count matched both the original RDKit
count and the reference with recalculated SSSR. Every output retained the
indexed atom table and gave the original isomeric structure when read back by
RDKit. The unrepaired RDKit path changed counts in 51 of these permutations.
These are fixed diagnostic inputs, not an estimate of real-world error rates
or a claim about all RDKit stereoperception.

The fixture recipe, provenance, license, minimization, exact oracle settings,
source/module hashes, per-structure results and normalized case hash are in
`validation/results/rdkit-atom-renumbering-2026-10-14.json`. No upstream
corpus, source code or attachment was copied.

```sh
cargo test -p chematic-core -p chematic-chem -p chematic-smiles -p chematic-smarts --lib --offline
cargo test -p chematic-chem --test renumber_stereo_contract --offline
python -m pytest crates/chematic-py/tests/test_renumber_atoms.py crates/chematic-py/tests/test_rdkit_compat.py
python scripts/validation/validate_atom_renumbering.py
```

Rust passed 1,721 tests (1,709 library plus 12 new contracts), with two
pre-existing ignored tests. Python passed 190 tests (28 new, 162 existing
compatibility tests). The 838-case comparison uses RDKit 2026.09.1 with
`removeHs=False` and temporarily scoped legacy stereo perception; it restores
the process setting afterward. Clippy with warnings denied, formatting,
refactor budget, binding inventory, naming and strict documentation checks
also passed. The integration target is registered in the CI matrix.

Implementation is in an isolated source copy preserving the preceding query-H
and PQR candidates. This feature's portable patch contains only atom-numbering
changes relative to that snapshot. It is uncommitted and unpublished; these
measurements concern local source.
