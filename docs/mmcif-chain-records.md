# mmCIF model and chain records

`chematic.read_mmcif_chains(text)` provides typed, immutable atom-record
partitions keyed by `(model_num, label_asym_id)`. It addresses the chain-wise
access requested in open [Open Babel #2458](https://github.com/openbabel/openbabel/issues/2458).
The existing general residue parser retains the required identities but had no
chain-reading API. This addition uses that native parser and has no Open Babel
or Gemmi runtime dependency.

```python
from chematic import read_mmcif_chains

text = """data_example
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
C C1 LIG Alpha 1 0.125 0 0 1
O O1 HOH Solvent 1 8.5 0 0 1
C C1 LIG Alpha 1 0.25 0 0 2
"""
result = read_mmcif_chains(text)
assert [(c.model_num, c.label_asym_id) for c in result.chains] == [
    (1, "Alpha"), (1, "Solvent"), (2, "Alpha")
]
atom = result.chains[0].atoms[0]
assert atom.source_index == 0
assert atom.position_angstrom == (0.125, 0.0, 0.0)
assert atom.residue.label_comp_id == "LIG"
assert result.source_text == text
```

The grouping uses [wwPDB label asym identifiers](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_atom_site.label_asym_id.html).
Every atom must have a non-placeholder label chain ID. Author chain names are
retained independently; they never supply a missing label ID or merge different
label chains. Different models stay separate. Chains appear in first-seen
order, and each chain's atoms retain their relative source-row order even if
rows from different chains or residues are interleaved. `source_index` is the
zero-based row index within the native reader's first `_atom_site` loop.

`MmcifChainAtom` retains serial number, ATOM/HETATM field, element, both atom
name aliases, the full alternate-location ID, Cartesian angstrom coordinates,
occupancy, B factor and formal charge. Its `MmcifChainResidue` retains both
component names, both chain aliases, entity ID, both sequence namespaces,
model and the full insertion code. Missing identifiers remain `None`;
non-polymer `label_seq_id=.` does not become residue zero. Long component
names and alternate IDs are not truncated. The existing numeric reader's
defaults apply (for example, missing occupancy becomes 1 and B factor becomes 0).
Typed values do not preserve the lexical distinction between `.` and `?` or
original numeric spelling; the complete original source remains available.

The result, chain tuples, atom tuples, coordinates and residue objects are
immutable. Source-only categories, comments, cell/symmetry information and
unknown column values remain in `source_text`; `unhandled_columns` reports the
native reader's unhandled atom-site column names. No CIF text is rewritten.

This is an eager view over the existing reader's **first `_atom_site` loop**,
including its models; later loops/blocks are kept only in the source sidecar.
It is not a multi-block iterator. The source must satisfy the native residue
reader's identity checks, including integer author sequence numbers, complete
rows and unique atom-name/alternate-location pairs within each residue.
Malformed/ambiguous inputs raise errors rather than returning partial results.

Resource limits default to and cannot exceed 8 MiB UTF-8 input, 50,000 atoms,
4,096 model/label-chain pairs and 8,192 bytes per native input line. Callers
can lower these with `max_input_bytes`, `max_atoms`, `max_chains` and
`max_line_len`; limits accept integer index values and reject Python bools.
The API accepts text contents only, not paths or bytes, and returns no partial
result after failure. It adds no graph bonding or biological assembly inference:
a chain partition is not a claim that its atoms form one connected molecule.

## Validation

The new tests include all 720 permutations of six interleaved source atoms,
separate models with repeated chain names, shared author aliases, long residue
names, multiple alternate locations, insertion codes, waters, immutable outputs,
unknown-column/source retention and explicit input/resource failures.

An independent Gemmi 0.7.3 CIF parser cross-checks 96 seeded synthetic cases
with permuted columns, quoted names, wrapped rows and shuffled atoms: all
384 partitions and 4,992 atom/residue records agree, with zero observed
coordinate difference. A local Open Babel 3.2.1 control reports one record
for the authored two-model/four-partition fixture. This records its default
input granularity; it does not establish that aggregation is always a bug,
and connected-component splitting is a separate operation.

Evidence and the reproducible oracle are in
`validation/results/mmcif-chain-records-2026-10-18/`. Python tests use a verified
existing local native extension; the native mmCIF parser and residue-binding
sources are unchanged from that fixture. No new current-main native build,
published-wheel, biological assembly or large real-structure validation is
claimed.

## Python reference

::: chematic._mmcif_chains
