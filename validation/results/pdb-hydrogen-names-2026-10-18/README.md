# Distinct checked PDB hydrogen names

Base `d3cc76b828c4732c49413d9a5fcf19b29f8bee11`.
[Open Babel #2945](https://github.com/openbabel/openbabel/issues/2945) is OPEN
and requests distinct names for H atoms added to PDB. No dedicated naming API
exists in the base. This adds Python `name_pdb_hydrogens_checked`, frozen
`PdbHydrogenNaming` / `PdbHydrogenName`, and public typing/documentation.

## Behavior

- Name existing explicit element-H rows H1..H999 per residue, in source order.
  Optional zero-based ATOM/HETATM row indices select only specified H atoms.
  Reserve all unselected atom names first, including future heavy/H rows.
- Residue identity: TER segment, chain, numeric residue number, insertion code,
  residue name. Mixed ATOM/HETATM rows use the same identity. Numeric spelling
  differences do not create another residue; interleaved residue rows share names.
- Change only columns 13..16 of selected atom rows and associated ANISOU rows.
  Other bytes, record order, serials, coordinates, fields, CONECT, TER/MODEL,
  line endings and final-newline presence are preserved exactly. Repeated
  selection is idempotent. Source text and selected old/new names are retained.
- Reuse unchanged native checked PDB selection for complete source/output
  validation. Unknown metadata refuses, with no automatic exclusion consent.
  Empty/multimodel input, missing elements, selected H altloc, non-H/duplicate/
  out-of-range indices and name exhaustion refuse. No partial result is returned.
- Bounds: 16 MiB input/output, 99999 atoms, 1000000 records, 1024 bytes per line.
  No H generation, parent inference, chemical/template naming, record-kind
  reassignment, D/T naming, file I/O, Rust/WASM/CLI API or performance claim.

## Tests and independent verification

- **241 Python tests passed**, including **53 new naming cases**, 26 unchanged
  PDB-selection cases and 162 compatibility cases. Covers 9/10/99/100/999 names,
  complete exhaustion and reserved-name exhaustion, scalar/index types, input/
  reference errors, explicit models, isotope records, altloc, frozen results,
  byte/atom/name limits and exact line preservation.
- Independent authored recipe: **2048 operations**, 8 atom rows each, four row
  orders, two residue orders, LF/CRLF, implicit/explicit model and all 64 H subsets.
  Independent per-residue available-name lists validate **6144 assignments**;
  **53248 complete output lines** checked, including **4992 changed ANISOU** rows.
  All 2048 idempotence and 32 default/all-H-selection checks pass.
- **RDKit 2026.09.1 reads 16384 output atoms**, checking name/serial, elements,
  record kind, residue/chain/insertion/altloc, occupancy/temp and exact
  parsed coordinates against the same source. This synthetic corpus uses C/H;
  D/T byte preservation is tested separately, not an isotope-perception oracle.
- **1024 implicit-model bond comparisons** also match all seven declared bonds.
  In the 1024 explicit-model cases, RDKit ignores the global CONECT records after
  ENDMDL and produces zero bonds for both source and output. These are coordinate/
  identity cross-reads, not successful seven-bond topology reads. Every original
  CONECT byte and native checked reference remains verified in all cases.
- Strict mypy passes the helper and a public consumer; five invalid calls/
  mutations yield five errors. Imported legacy package stubs are followed
  silently, not checked package-wide. Source stub/cross-binding/binding/unsafe/
  refactor/name gates, strict MkDocs, one executed guide block and three API
  anchors pass. New helper syntax also parses under Python 3.9 grammar.

## Open Babel control

Installed **Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27**:

```sh
/opt/homebrew/bin/obabel -:CO -opdb -h
```

Emits six atoms, including four H rows all named ` H  `. All six rows are HETATM
in this control, so the issue's mixed ATOM/HETATM observation is not reproduced
here. This is an independently authored methanol-SMILES probe, not the original
crystallographic PDB pipeline. It does not establish meaningful 3D geometry.
Raw output is `openbabel-generated.pdb`, invocation/version in
`openbabel-control.json`. A separate explicit `select_pdb_atoms(...,
allow_metadata_loss=True)` projection excludes COMPND/AUTHOR/MASTER; all three
losses are recorded in `oracle.json`. The naming API itself excludes no records.
Its coordinate source and named result are saved separately. H1..H4 are distinct
and all atom fields outside their name columns remain unchanged.

## Native fixture and reproduction

Pure Python change, no Rust/native rebuild. A hard-linked existing extension
from **ca7e14e8** is used without modification for package import and unchanged
PDB validation. Its SHA-256 matches the prior build provenance; selected native
PDB source paths are byte-identical to this base. The fixture still has its own
recorded build revision; this is not a current-main native build, published wheel
or older-Python runtime validation. Hashes and source paths: `provenance.json`.
Python 3.13.6 / pytest 9.1.1; NumPy 2.5.3 is already a package dependency.

Set `CANDIDATE_PY` to this source package plus the verified extension and
`RDKIT_PY` to RDKit 2026.09.1.

```sh
PYTHONPATH="$CANDIDATE_PY" python3.13 -m pytest \
  crates/chematic-py/tests/test_pdb_hydrogen_names.py \
  crates/chematic-py/tests/test_pdb_selection.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONPATH="$CANDIDATE_PY:$RDKIT_PY" python3.13 \
  validation/results/pdb-hydrogen-names-2026-10-18/oracle.py
MYPYPATH="$CANDIDATE_PY" python3.13 -m mypy --no-incremental \
  --config-file validation/results/pdb-hydrogen-names-2026-10-18/mypy.ini \
  crates/chematic-py/python/chematic/_pdb_hydrogen_names.py \
  validation/results/pdb-hydrogen-names-2026-10-18/typing_valid.py
# typing_invalid.py: expect five errors, exit 1.
python3.13 -m mkdocs build --strict --site-dir /tmp/chematic-pdb-names-site
```

The fixed-column recipe follows the
[wwPDB 3.3 coordinate specification](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html).
The source checkout's four unrelated working document edits are preserved.
