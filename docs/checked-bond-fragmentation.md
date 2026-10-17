# Checked bond fragmentation and cut counters

`fragment_on_bonds_checked` replaces selected bridge single bonds with paired
isotope-labelled dummy atoms. Its `cuts_per_atom: list[int]` argument is checked
before editing, then overwritten with one cut count per original atom on success.
This implements a bounded part of the mutable-output typing work tracked in
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601).

```python
import chematic as ch

source = ch.from_smiles("CCOC")
counts = [99] * source.atom_count
result = ch.fragment_on_bonds_checked(source, [2, 0], cuts_per_atom=counts)
assert counts == [1, 1, 1, 1]
assert result.cuts_per_atom == (1, 1, 1, 1)
assert result.molecule.atom_count == 8
assert source.atom_count == 4
assert [(a.source_bond, a.anchor, a.dummy, a.isotope)
        for a in result.attachments] == [
    (2, 3, 4, 2), (2, 2, 5, 3), (0, 1, 6, 0), (0, 0, 7, 1),
]
```

Indices refer to the source molecule's native atom/bond order. Original atoms,
hydrogen counts, maps, charges and isotopes remain unchanged. Each selected bond
adds two dummy atoms, in the supplied cut order. A `(begin_label, end_label)` pair
places the begin label's dummy on the original end atom, and vice versa. Default
labels are those original atom indices. Label zero means no isotope label.
Python molecule properties are copied into an independent property store.

```python
source = ch.from_smiles("CC(C)O")
result = ch.fragment_on_bonds_checked(
    source, [0, 1, 2], dummy_labels=[(7, 9), (0, 65535), (55, 55)]
)
assert result.cuts_per_atom == (1, 3, 1, 1)
assert [a.isotope for a in result.attachments] == [7, 9, 0, 65535, 55, 55]
```

The frozen `BondFragmentation` and `FragmentAttachment` records expose immutable
counts and provenance. The returned `Mol` remains an ordinary molecule object.
An optional counter must be an exact plain `list[int]`, with one entry per source
atom. Initial counter values are ignored. Invalid inputs leave the list and
source molecule unchanged; counters are replaced, never accumulated.

```python
source = ch.from_smiles("CCO")
counts = [17, 17, 17]
try:
    ch.fragment_on_bonds_checked(source, [0, 99], cuts_per_atom=counts)
except ValueError:
    pass
else:
    raise AssertionError("invalid cut was accepted")
assert counts == [17, 17, 17]
```

## Supported profile

- At most 4096 source atoms, 8192 source bonds and 1024 distinct cuts.
- Cut bonds must be plain single bonds outside all rings. Dummies always use
  single bonds. Custom isotope labels are integers from 0 through 65535.
- Nonempty edits reject stereo annotations anywhere in the molecule, existing
  dummy/query atoms, and nonstandard bond types. Remaining double, triple and
  aromatic bonds can be retained.
- Sequence arguments accept eager Python sequences such as list, tuple and
  range. Lazy iterators and text are refused; integers exclude `bool`.
- Empty cuts return a clone, zero counters and no attachments, retaining any
  source annotations. Size limits still apply. RDKit rejects empty cut lists.
- Coordinates, ring cutting, stereo transfer, alternative dummy bond orders,
  and cutting without dummies are outside this API's profile.

Use `result.molecule.smiles` for fragment interchange. SMILES graph identity and
V3000 indexed atom/bond fields were checked against RDKit 2026.09.1 for 30 source
molecules, 780 nonempty cut/label combinations and 90 defined no-op cases.
V3000 dummy query metadata is outside that comparison.

**Export limitation:** existing V2000 and MolJSON writers do not preserve these
dummy graphs reliably. V2000 can misencode dummy isotope mass offsets; MolJSON
can emit the wildcard's placeholder element as carbon. Avoid those formats for
fragment interchange. This API does not change the writers. The source, commands
and scope of the oracle are recorded in
`validation/results/bond-fragment-counters-2026-10-17/` in the repository.

The Rust entry point is
`chematic_chem::bond_fragment::fragment_single_bonds_checked`; it returns a
`SingleBondFragmentation` with the native graph, counts and attachments, or a
`BondFragmentError` without mutating its input. The Python API requires no RDKit
installation.
