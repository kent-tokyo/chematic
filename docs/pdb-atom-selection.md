# PDB atom selection without identity changes

The isolated candidate adds `select_pdb_atoms` in Rust/Python for the workflow
reported in [Open Babel #2481](https://github.com/openbabel/openbabel/issues/2481):
removing ligand atoms must not change names or residue identity of retained
protein atoms. The issue was verified **OPEN** on 2026-10-15. A separate
original six-atom reproduction also shows this behavior in Open Babel 3.2.1.

The API selects original PDB atom records directly, reusing chematic's existing
strict coordinate parser. It performs no residue, element, bond or hydrogen
perception. This is an immutable single-model coordinate-selection operation,
not an implementation of Open Babel's mutable molecule API.
A later [coordinate-block merger](pdb-coordinate-merge.md) reuses this validator
to assemble separate sources with explicit serial and chain mappings.

## API and example

```python
import chematic

# None selects all and exposes identities in original ATOM/HETATM row order.
view = chematic.select_pdb_atoms(pdb_text, allow_metadata_loss=True)
keep = [a["source_index"] for a in view["atoms"] if a["residue_name"] != "LIG"]
selected = chematic.select_pdb_atoms(
    pdb_text, keep, allow_metadata_loss=True,
    max_atoms=99_999, max_input_bytes=64 << 20,
)
filtered_pdb = selected["pdb"]
```

Use `allow_metadata_loss=False` (the default) for a supported coordinate
block. The example explicitly consents to excluding full-entry metadata such
as HEADER/SEQRES, which cannot be safely kept stale after arbitrary selection.
Every excluded record is listed in `losses`, and the entire original text
remains in `source_text`.

Rust exports `PdbSelectionLimits`, `PdbSelectionError`, `PdbSelectedAtom`,
`PdbSelectionLoss`, `PdbSelectionReport`, `select_pdb_atoms` and
`select_pdb_atoms_with_limits`. Indices are `Option<&[usize]>`: `None` keeps
all source atoms and `Some(&[])` keeps none. Python accepts the corresponding
`None` or sequence of zero-based indices. Negative/outside/duplicate indices
are rejected. Input index-list order does not reorder the output.

The report contains:

- `pdb`: selected coordinate block, with LF and one final END.
- `source_text`, `source_atom_count`, `selected_source_indices`, `old_to_new`:
  original input and complete row-index mapping; removed rows map to `None`.
- `atoms`: selected source index, original serial, record type, raw four-column
  atom name, alternate location, residue name, chain, residue number, insertion
  code and `position_angstrom`.
- `model_num`: the optional original explicit model number.
- `losses`: source line/record name for excluded unsupported records.
- `retargeted_ter_lines`, `filtered_conect_lines`: original line numbers whose
  supported associated records were changed or filtered by the selection.

These are source ATOM/HETATM row indices. They do not refer to a sanitized
RDKit model or to atom serial numbers. No implicit hydrogen removal changes
the index range. Returned objects are independent of the input and later calls.

## Preserved and updated records

The independently written selection logic follows the
[wwPDB 3.3 coordinate layout](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)
and [CONECT layout](https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html).
No Open Babel source or original issue attachment was copied.

| Record | Behavior |
|---|---|
| ATOM/HETATM | Retain exact line contents, original serial and source order for selected atoms |
| ANISOU | Retain exact line contents only when its associated atom remains |
| TER | Omit empty segments; otherwise follow the last retained atom's residue fields and serial + 1 |
| CONECT | Retain only selected source/target atoms, preserving target order and duplicate-target multiplicity |
| MODEL/ENDMDL | Preserve one complete explicit model; omit wrapper for an empty selection |
| END | Emit one final END; missing source END is accepted |
| Other nonblank records | Reject by default; with consent, exclude and report their original line/name |

ATOM/HETATM and retained ANISOU line contents are unchanged, including
coordinate spelling, occupancy, temperature factors, record type, alternate
location, insertion code, element and charge fields. Only line endings are
normalized. Unknown fields are not inferred from atom names. Selected atom
serials are never renumbered.

TER is validated against its source segment endpoint and then updated when
needed. Its serial is derived from the last retained atom. Collisions with
retained atoms or other generated TER records are errors; the selector never
renumbers atoms to hide a collision. Segment boundaries therefore remain
explicit even when a segment's original last atom is removed.

CONECT supports a source and up to four covalent target fields per line.
Repeated targets are retained as repeated targets, and reciprocal edges are
not invented. Original lines are retained if all their referenced atoms remain;
otherwise targets are filtered and the necessary fixed-width line is emitted.
References to removed atoms are eliminated before a newly derived TER serial
could reuse their number. Unknown/self references are errors.

## Validation and boundaries

Inputs are checked before selection for printable ASCII/spaces, finite XYZ,
positive decimal atom serials, nonempty atom/residue names, integer residue
numbers and finite optional occupancy/temperature values. Atom serials must
be unique. ANISOU must follow its referenced atom, agree on the fixed-column
identity and carry six integer tensor fields. Duplicate ANISOU is refused.

The API accepts either one implicit coordinate model or one complete explicit
MODEL/ENDMDL pair. Multiple/nested/incomplete models, coordinates after
ENDMDL/CONECT, data after END, invalid references, hybrid-36 serials and legacy
CONECT fields beyond the four covalent neighbors are refused. Metadata-loss
consent does not override errors in supported records.

Hard limits are 64 MiB input/output, 1024 bytes per physical input line,
99,999 atom rows and 1,000,000 physical records, including blank lines. Rust
limits and Python keyword arguments can tighten these caps. Empty selections
produce `END\n`; the report still identifies the original source model.

This exports a coordinate block rather than a complete deposition-ready PDB
entry. It does not repair chemistry, validate a residue template, regenerate
SEQRES/LINK/SSBOND/MASTER or preserve those annotations without explicit
exclusion consent. For long identifiers and richer residue namespaces, retain
mmCIF and use the existing [residue reader](mmcif-residue-identity.md).

## Open Babel reproduction

The original `basic.pdb` fixture has two incomplete two-atom protein residue
fragments (VAL A 101, LEU B 152) and a two-atom LIG L 7, without CONECT.
It is a minimal metadata example, not a complete protein or physical model.

The public-API C++ oracle reads that PDB, snapshots `GetResidue()` identities,
calls `OBMol.DeleteAtom` in descending original index order, snapshots again
and writes PDB. The dynamically linked Open Babel reports **3.2.1**.

| Control | Changed retained atom identities |
|---|---|
| No deletion, before/after serialization | 0 |
| Delete original ligand rows 4 and 5 | 4 |

After deletion, Open Babel changes all four retained VAL/LEU residue names
to UNL, clears their chain IDs, changes their atom names and assigns new
residue numbers. The unchanged no-deletion control distinguishes the deletion
trigger from ordinary PDB serialization. The native selector retains the
original four atom lines exactly. The original issue's linked attachment was
not fetched; the claim is scoped to this independent synthetic reproduction.

## Measured checks, 2026-10-15

The richer original fixture includes ANISOU, alternate location A, a negative
residue number, insertion code I, a blank chain, charge text, TER and duplicated
CONECT targets. Six residue orders, eight within-residue orders, LF/CRLF and
implicit/explicit model variants produce 192 sources. Atom serials and
associated references are remapped consistently by the independent recipe.
All 64 subsets of each six-atom source are checked.

| Check | Result |
|---|---|
| Exact subset operations | 12,288 |
| Exact retained atom / ANISOU lines | 36,864 / 36,864 |
| Updated TER / filtered CONECT records checked | 9,216 / 55,296 |
| RDKit nonempty cross-reads / empty checks | 12,096 / 192 |
| Unexpected mismatches | 0 |

RDKit **2026.09.1** independently reads the outputs with `sanitize=False`,
`removeHs=False`, `proximityBonding=False`. Atom serial/name, residue/chain,
alternate location, insertion code, record type and coordinates match exactly.
Native checks also cover ANISOU contents, connectivity multiplicity, TER
endpoints, index maps, source immutability and idempotent reselection.

Local Rust validation passed 753 unit tests plus 37 integration tests (790
total). Python 100 tests passed, including checked mmCIF/PDB, residue and
DMol3 regressions. Scoped Clippy with warnings denied, formatting, unchanged
refactor allowances, binding inventory, product naming and strict MkDocs are
recorded in the report. No hosted CI, publication or complete upstream fix
is claimed.

## Reproduction and delivery

```sh
cargo test -p chematic-mol --lib --test pdb_selection_contract \
  --test mmcif_pdb_checked --test mmcif_residue_contract
python -m pytest crates/chematic-py/tests/test_pdb_selection.py
c++ -std=c++17 -O1 \
  validation/openbabel_issues/pdb-selection/openbabel_delete_probe.cpp \
  $(pkg-config --cflags --libs openbabel-3) -o /tmp/ob2481-probe
python scripts/validation/validate_pdb_selection.py --probe /tmp/ob2481-probe
```

Install the candidate Python extension and pinned RDKit first. The C++ helper
is an external Open Babel validation oracle, never a chematic dependency.
The fixture manifest records provenance, MIT OR Apache-2.0 licensing,
redistribution, minimization, oracle settings and normalized payload hashes.
The fixture-specific Git attributes preserve LF and allow the intentional
trailing spaces required by the fixed-column PDB records.

Source candidate: `/private/tmp/chematic-pdb-selection`, inheriting
`/private/tmp/chematic-dmol3`. Incremental patch:
`validation/openbabel_issues/pdb-selection-2026-10-15.patch`. Evidence:
`validation/results/openbabel-pdb-selection-2026-10-15.json`.
Only this document, patch and report are copied to the primary checkout.
Unrelated source work is preserved; the feature is uncommitted and unpushed.
