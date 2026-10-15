# mmCIF residue identity for non-polymers

This local source candidate addresses the residue-grouping workflow in
[Open Babel #3028](https://github.com/openbabel/openbabel/issues/3028).
The existing chematic mmCIF reader already retained per-atom charge/position
and residue fields; its separate protein hierarchy deliberately excludes
glycans, ligands and water. The new API constructs a general residue view
without using that amino-acid-only projection.

## Use

```python
import chematic

with open("structure.cif") as source:
    document = chematic.parse_mmcif_residues(source.read())

for residue in document["residues"]:
    print(residue["label_comp_id"], residue["auth_comp_id"],
          residue["auth_seq_id"], residue["atom_indices"])

# Every atom points back to a residue; indices refer to original atom_site rows.
for atom_index, residue_index in enumerate(document["atom_to_residue"]):
    assert atom_index in document["residues"][residue_index]["atom_indices"]
```

Rust exports `parse_mmcif_residues` and `parse_mmcif_residues_with_limits`,
returning `Result<MmcifResidueDocument, MmcifResidueError>` from `chematic_mol`.
The latter takes the existing `MmcifParseLimits` (input bytes, line bytes and
atom rows). Python exposes the same limits as keyword arguments and returns
`atoms`, `residues`, `atom_to_residue`, `cell`, `space_group` and
`unhandled_columns`. Invalid input raises `ValueError`.

## Identity contract

A residue identity includes model number, both label/author chain identifiers,
entity id, both sequence identifiers, full insertion code and both component
names. Every stored identity field participates in grouping. The residue list
uses first appearance order; each `atom_indices` list keeps source row order.
Rows may be interleaved. No atom row is removed or reordered.

`label_seq_id` identifies polymer sequence positions and is often absent for
non-polymers. The author's sequence identifier can still distinguish NAG,
NAG and BMA residues when all three have `label_seq_id='.'`. Absence is never
replaced by residue zero in the grouping key. A literal author identifier
`'0'` remains valid and distinct from absence; `'01'` and `'1'` remain distinct.

Both chain namespaces are retained. Two label chains with the same author
chain and residue number remain separate. Different model numbers, entity
ids and insertion codes also stay distinct. Different component identities
at a shared position form separate residue variants, retaining their names.
Alternate locations with the same component identity stay separate atom rows
within one residue; their full codes are retained. Partial aliases are not
reconciled by guessing: rows group only when their complete stored identity
matches.

The identity definitions follow the wwPDB dictionary for
[`label_seq_id`](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_atom_site.label_seq_id.html)
and [`auth_seq_id`](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_atom_site.auth_seq_id.html).

## Strict boundaries

- At least one chain, sequence, component and atom identifier must be present.
  If both sequence identifiers are missing, the parser refuses to infer a
  residue boundary. It does not claim support for identities supplied only
  by branch-scheme categories outside atom_site.
- A repeated label atom name (or author name when no label is supplied) with
  the same alternate-location code in one residue is an error. Alternate
  locations with different codes are accepted.
- Incomplete atom rows, non-finite numeric values, unknown elements and
  resource-limit violations fail explicitly. Model numbers must fit positive
  `i32`; formal charges must fit `i32` before the legacy atom projection.
- The existing numeric atom reader requires integer `auth_seq_id`. The wwPDB
  dictionary also permits nonnumeric author identifiers; those remain
  unsupported and produce an explicit parse error here.

Rust's `parsed: MmcifResult` contains the legacy Cartesian/numeric record
view. Full residue aliases/insertion strings live in `residues[].identity`;
full atom aliases/alternate-location strings live in `atom_ids`. Use these
identity tables rather than the legacy record's collapsed aliases or
single-character fields. Python's `atoms` rows include both atom aliases and
full `alt_loc`/`icode` strings alongside the legacy numeric fields.

The legacy `parse_mmcif` API keeps its existing schema and behavior. The new
view does not infer bonds, glycan connectivity or stereochemistry and does
not export PDB or provide a lossless mmCIF writer. Calling the legacy writer
with only `parsed.atoms` retains that writer's documented alias-collapse
limitations. Long component names survive this reader; PDB's three-column
name restriction is a separate conversion concern.

## Reproduction and validation

The independent six-row fixture is
`validation/openbabel_issues/mmcif/branched-residues.cif`. It contains two atom
names for each of NAG/NAG/BMA, distinct author sequence identifiers 1/2/3,
missing label sequence ids, and distinct author/label chain names. It is
synthetic metadata for a grouping test, not a chemically complete glycan.
No deposited structure or issue attachment was copied. Redistribution is
allowed under the project's `MIT OR Apache-2.0` terms.

On this host, Open Babel 3.2.1 writes all six PDB rows as NAG in one residue.
The new reader returns three residues, keeping BMA's two rows as BMA. A
control that supplies numeric label sequence ids makes Open Babel preserve
NAG/NAG/BMA, isolating the missing-label-sequence condition.

All 720 permutations of the six input rows preserve the expected residue
partition, full identity and source-row mapping. This diagnostic is not a
measurement of real-protein frequency or downstream docking quality.
The regeneration script records provenance, license, minimization, oracle
version/settings, hashes, expected values, comparison output and the
normalized permutation-case hash in
`validation/results/openbabel-mmcif-residues-2026-10-14.json`.

```sh
cargo test -p chematic-mol --lib --test mmcif_residue_contract --test mmcif_pqr_integration --offline
python -m pytest crates/chematic-py/tests/test_mmcif_residues.py crates/chematic-py/tests/test_new_formats.py crates/chematic-py/tests/test_input_errors.py
python scripts/validate_mmcif_residues.py
```

Rust passed 772 tests (753 library, 15 new contracts, four existing mmCIF/PQR
integration tests). Python passed 96 (25 new, 71 existing format/error tests).
Clippy with warnings denied, formatting, refactor budget, static binding
inventory, naming and strict MkDocs passed. The Rust integration target is
registered in the CI matrix.

Implementation is in an isolated source copy retaining the preceding local
RDKit and PQR candidates. Its portable patch contains only this feature's
changes relative to that snapshot. It remains uncommitted and unpublished;
these measurements concern local source and the pinned external engine.
