# Python Cleanup output parity with RDKit

On the pinned 4,072-row stress corpus, RDKit 2026.09.1 accepts 4,069
comparable rows. A source-built CPython 3.13 wheel from `origin/main`
matched RDKit Cleanup's canonical SMILES on **4,028/4,069** rows. The
candidate wheel matches **4,053/4,069**: 25 rows improve and none regress
on the same input set. The [row-level record](../validation/results/python-cleanup-stereo-output-2026-10-21.json)
contains the corpus hash, both wheel hashes, fixed indices and all 16
remaining differences.

The Python binding called the ordinary RDKit-compatible SMILES writer on the
result of `Mol.cleanup()`. That writer models a fresh parse of the output
SMILES, which discards a non-tetrahedral tag on an isolated metal. RDKit's
in-memory Cleanup result retains the tag after disconnecting the metal's
ligands. The Rust source comparison already used a Cleanup-specific writer,
so its 4,053-row count exposed this binding-only gap.

`Mol.cleanup()` now records that its result needs the Cleanup output writer.
The ordinary `from_smiles("[Pt@SP1+4]").rdkit_smiles` still yields `[Pt+4]`.
Pickle payload v3 retains the output provenance and distinguishes a present
but empty stereo-neighbor order from absent order; the decoder still accepts
v1 and v2 payloads. Focused Python tests (14) and Rust pickle tests (3)
passed on this candidate. The final source wheel was compared across all
4,072 input rows using `scripts/rdkit_cleanup_parity.py --sample-size 4072`.

The **16 remaining rows** are in the Rust Cleanup core: three acid/base
reionization differences, two canonical representation differences, two metal
charge differences, eight metal disconnection differences, and one
vacant-site metal stereo difference. The vacant-site case cannot be repaired
by simply retaining its tag: RDKit changes the square-planar permutation on
disconnection, while the current result no longer carries the original
ligand order. This change does not establish full RDKit Cleanup parity or a
published-package result.
