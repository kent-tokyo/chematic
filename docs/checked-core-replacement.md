# Checked core removal with explicit atom maps

`replace_core_checked(mol, core, match)` removes a concrete core and places
isotope-labelled dummy atoms at its former side-chain attachments. `match[i]`
is the source atom corresponding to core atom `i`. This implements the explicit,
typed match-vector part of [RDKit #9601](https://github.com/rdkit/rdkit/issues/9601),
with a bounded native profile corresponding to `ReplaceCore(labelByIndex=True)`.

```python
import chematic as ch

source = ch.from_smiles("CCOC")
core = ch.from_smiles("CO")
result = ch.replace_core_checked(source, core, [1, 2])
assert result.original_to_output == (0, None, None, 1)
assert result.molecule.atom_count == 4
assert [(a.core_atom, a.source_core_atom, a.source_anchor, a.source_bond,
         a.anchor, a.dummy) for a in result.attachments] == [
    (0, 1, 0, 0, 0, 2), (1, 2, 3, 2, 1, 3),
]
assert source.atom_count == 4
```

`CoreReplacement` and `CoreAttachment` are frozen records. The result's `Mol`
is independent and keeps source Python properties in a copied store. Surviving
source atoms retain their relative order, attributes, H counts, atom maps and
caller tags. `original_to_output` has one entry per source atom: its new index,
or `None` if removed. Dummies follow retained atoms in source boundary-bond order.
Their isotope is the **core atom index**, so isotope zero is unlabelled.

```python
import chematic as ch

source = ch.from_smiles("Cc1ccccc1CO")
core = ch.from_smiles("c1ccccc1")
result = ch.replace_core_checked(source, core, range(1, 7))
assert [a.core_atom for a in result.attachments] == [0, 5]
assert result.original_to_output == (0, None, None, None, None, None, None, 1, 2)
```

The source bond/core/anchor indices in attachment records refer to the original
input. `anchor` and `dummy` refer to the returned graph. Reversing the matching
vector for a symmetric core may change labels; the supplied correspondence is
respected, with no matching search or canonicalization of that choice.
Disconnected concrete cores and unrelated source components are supported.
Removing the entire source returns an empty `Mol` with an all-`None` map.

## Checked profile

- Source: at most 4096 atoms and 8192 bonds. Core: 1..1024 atoms and at most
  8192 bonds. At most 1024 core/side-chain boundary bonds.
- One eager sequence of distinct integer source indices, with exactly one entry
  per core atom. Lists, tuples and ranges work; generators, text and bool indices
  are refused. Neither input nor the supplied sequence is modified.
- Corresponding atoms must have equal native element, aromatic flag, charge and
  isotope. Core bonds must exist with equal native order in the source. H counts
  and atom-map numbers are ignored **during matching**; retained source values
  survive in the result. Additional edges among mapped core atoms are permitted.
- Both graphs must be concrete and stereo-free, with no existing dummy/query
  atoms or nonstandard bond types. Only non-ring plain single bonds may cross
  the core boundary. Entire cyclic cores can be removed via exocyclic single
  attachments; cutting through a ring is refused.
- No automatic aromaticity perception is performed. A Kekulé source and an
  aromatic core can be rejected despite chemical equivalence. Provide matching
  native representations; the API does not normalize them silently.

Invalid types raise `TypeError`; invalid mappings, size limits and unsupported
graphs/boundaries raise `ValueError`. Automatic substructure matching, query
cores, alternative dummy labelling/bond policies, `replaceDummies`,
`requireDummyMatch`, stereo transfer and batch matching are outside this API.

Use `result.molecule.smiles` for fragment interchange. Existing V2000 and MolJSON
writers do not reliably preserve isotope-dummy graphs; see
[checked bond fragmentation](checked-bond-fragmentation.md). The oracle checks
SMILES graph identity and V3000 indexed atomic/bond attachment fields; V3000 dummy
query metadata and coordinates are outside that comparison.

RDKit 2026.09.1 comparisons cover 415 accepted explicit maps from 32 source
molecules and 14 core graphs, including four empty outputs. Another 12 mappings
verify the explicit Kekulé/aromatic representation refusal. Reproduction and
source/native hashes are under
`validation/results/checked-core-removal-2026-10-17/` in the repository.

Rust uses `chematic_chem::core_replacement::replace_core_checked` and returns a
`CoreReplacement` or `CoreReplacementError`. It reuses the checked single-bond
fragmenter and order-preserving native atom deletion. Python requires no RDKit.
