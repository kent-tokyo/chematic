# Open Babel #3005: rigid PDBQT isotope boundary

Local source candidate based on `38b4c92c430d9b3bfe70e47e97293b9017ced990`.
This is not a release, docking score comparison or flexible-ligand preparation
implementation. Upstream issue: https://github.com/openbabel/openbabel/issues/3005.

## Result

- `isotope-projection.json`: 23 Python and real Node/WASM projections agree
  exactly, with 20 labelled inputs refused by the default checked writer and
  all 48 isotope annotations recorded in the explicit projection sidecar.
- Source element order and isotope annotations are checked independently by
  RDKit 2026.09.1 (`removeHs=False`). Output serials, element types, coordinates
  and supplied charges are decoded directly from fixed columns, without using
  chematic's parser as the only observer. This is not validation of docking-type
  chemistry or a lossless molecular graph round trip.
- `baseline-checked.pdbqt`: the previous checked writer accepts the mixed D/13C/T
  input and omits labels silently. The baseline runtime was built on the RDKit
  #9738 branch (`c7e348f615cf44a97fa2e19ec853eab7d2e2670d`); its native/Python
  PDBQT source files match the main base. The separate baseline JSON records
  its extension hash, rather than calling it a published main build.
- `openbabel-reproduction.json`: Open Babel 3.2.1 reproduced ordinary versus d8
  acetophenone as 9 versus 17 output atoms, 0 versus 8 `HD` atoms, `TORSDOF 1`
  versus `TORSDOF 2`, with exit status 0 and no warnings. Carbon-13-labelled
  acetophenone emits 9 atoms and `TORSDOF 1`. All input graphs have 17 explicit
  atoms with the same coordinates. SDFs were independently generated with
  RDKit AddHs and EmbedMolecule seed `0xC0FFEE`, then only isotope labels changed.

The exact input SDFs and Open Babel output blocks are retained here; stderr is
recorded in `openbabel-reproduction.json`.
They are original fixtures generated for this repository, not issue attachments.

## Reproduce

Build the candidate Python extension and a Node-target wasm-bindgen package.
Use paths to those artifacts, not an installed release:

```sh
PYTHONPATH=/path/to/candidate-python python3 scripts/validation/validate_pdbqt_isotopes.py \
  --wasm-module /path/to/node/chematic_wasm.js \
  --rdkit-python /path/to/rdkit-2026.09.1/bin/python \
  --output /tmp/pdbqt-isotope-projection.json

obabel -isdf validation/openbabel-3005-2026-10-24/acetophenone-ordinary.sdf -opdbqt
obabel -isdf validation/openbabel-3005-2026-10-24/acetophenone-d8.sdf -opdbqt
obabel -isdf validation/openbabel-3005-2026-10-24/acetophenone-C13.sdf -opdbqt
```

The candidate report records source hashes and Python/WASM binary hashes. Its
`head` is the pre-commit base; the tested dirty source is identified by those
hashes. Timestamps and host paths are provenance, not portable build requirements.

Native library tests and PQR/PDBQT contracts passed. The new native contract has
5 tests; new Python contract has 31 tests (217 including the related checked
PDBQT/PQR and atom-data suites). Both new and existing Node/WASM PDBQT contracts
passed. A paired 17-atom acetophenone/d8 check verifies identical rigid output,
eight nonpolar `H` types, unchanged supplied charges and no `TORSDOF`.
Clippy with warnings denied, formatting, stubs, public API inventory/compatibility,
cross-binding manifest, refactor budget, script layout, integration-test discovery
and strict MkDocs build passed. The strict docs build uses tracked documentation
plus the new page, excluding unrelated ignored/private workspace files.
