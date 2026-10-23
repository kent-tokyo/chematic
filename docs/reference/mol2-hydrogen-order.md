# Grouping MOL2 hydrogens with their parents

`order_mol2_hydrogens_checked(text)` places each explicit hydrogen row after
its bonded non-H parent in one MOL2 record. It adds the ordering option
requested in [Open Babel #2733](https://github.com/openbabel/openbabel/issues/2733),
verified open on 2026-10-19. Open Babel's current heavy-atoms-first output is
not inherently chemically wrong; this is an explicit alternative for workflows
that need parent/H row groups.

```python
from pathlib import Path
import chematic as ch

source = Path("input.mol2").read_text()
result = ch.order_mol2_hydrogens_checked(source)
# Inspect or explicitly save result.mol2.
print(result.new_to_old)       # Output row -> original zero-based row.
print(result.old_to_new)       # Original row -> output zero-based row.
print(result.atom_id_map)      # (original MOL2 ID, output MOL2 ID).
print(result.hydrogen_parents) # (original H ID, original parent ID).
```

The frozen `Mol2HydrogenOrder` also retains exact `source_text` and the output
`mol2`. No files are written. Non-H atoms retain their relative source order;
hydrogens sharing a parent retain their relative source order. Only explicit
H atoms are moved. The function does not add H atoms, infer bonds from distance,
rename atoms or generate NMR restraint names.

## Atom IDs and references

Moving atom rows alone can make a consumer attach bonds to the wrong atoms.
Whenever row order changes, output atom IDs become `1..N`. The function updates
all supported references together:

- BOND atom endpoints, retaining bond IDs, types and other fields;
- SUBSTRUCTURE root atom IDs (zero root remains zero);
- charge-only UNITY_ATOM_ATTR group atom IDs.

Coordinates, names, atom types, partial-charge tokens, substructure membership,
status fields and bond-row order are retained. Replaced ID tokens may change
width; other whitespace and physical line endings remain exact. Blank/comment
lines keep their physical slots and are not treated as atom annotations. The
atom rows must share LF or CRLF endings and each must be terminated, preventing
row concatenation when the final atom moves.

If row order is already correct, the exact source text and original IDs are
returned, including sparse IDs. A second call on the output is an exact no-op.
External atom-indexed data must be reordered using the returned mappings; this
API cannot update files or references outside the supplied text. It does not
promise that every consumer accepts sparse IDs in an unchanged source.

## Supported profile

The source must contain one nonempty molecule, at most 16 MiB, 16,384 atoms,
32,768 bonds and 131,072 physical lines. A physical line is at most 16,384
UTF-8 bytes. Supported sections are MOLECULE, ATOM, BOND, SUBSTRUCTURE and
UNITY_ATOM_ATTR with exactly one `charge` attribute per group. Unknown or
duplicate sections/attributes refuse because they may contain references whose
meaning is not modeled. SUBSTRUCTURE IDs and nonzero roots must resolve.

Hydrogen types are `H` and `H.*`; each requires exactly one ordinary single
bond (`1`) to a non-H atom. Isolated H, H2, bridging H and non-single H bonds
refuse. Implicit hydrogens are not materialized. Existing source-parser element,
coordinate, ID, count and formal-charge validation also applies. The function
does not establish that the supplied force-field types, valence or charges are
chemically correct.

After an edit, the native reader checks mapped atom/bond metadata, coordinates,
formal charges and header properties against the original. Wrong types raise
`TypeError`; malformed, unsupported or oversized input raises `ValueError`.
A mapped readback mismatch raises `RuntimeError`; no partial result is returned.

## Validation

134 local Python tests pass, including 33 new ordering tests. Tests cover
sparse IDs, partial/formal charge separation, root and attribute remapping,
line endings, idempotence, unsupported H/metadata, readback failure and the
actual 16,384-atom boundary.

The independent oracle uses 14 small molecules with 12 authored atom orders
each, including tetrahedral and E/Z examples. RDKit 2026.09.1 verifies all 168
mapped atom/bond tables, grouped H rows, coordinates and Standard InChI
identities. Coordinate differences after ordering are zero. Open Babel 3.2.1
rewrites all 168 outputs; subsequent independent RDKit reading retains the
mapped graph and Standard InChI, with maximum coordinate difference
`4.99e-5` Angstrom from Open Babel's four-decimal output.

Open Babel's `CCO -h -omol2` output also reproduces the heavy-first arrangement
that this optional transformation changes. Raw evidence, commands and the
local native/import boundary are in
`validation/results/mol2-hydrogen-order-2026-10-19/`.
This is a Python source-editing capability over the existing native reader,
not a general MOL2 editor, native rebuild or published-release claim.
