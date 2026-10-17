# Open Babel 3.2.1: unmarked 3D exocyclic alkene follow-up

## Paired result

RDKit 2026.09.1 embedded 150 tetrahedral and 150 alkene-stereo ChEMBL
identities, removed CTAB stereo marks, and admitted each row only when its
own read-back of both V2000 and V3000 retained the source's canonical
isomeric SMILES. The 300 identities are disjoint, under the same RDKit
canonicalization, from the 740 identities in the three earlier browser and
ChEMBL 3D selections. CheMatic before and after this change and Open Babel
3.2.1 read the **same** SDF bytes. RDKit normalized each output for exact
mapless graph and stereo comparison.

| Input | CheMatic main `db0cf128` | CheMatic `7b63f432` | Open Babel |
| --- | ---: | ---: | ---: |
| V2000, exact identity / 300 | 298 | **299** | 296 |
| V3000, exact identity / 300 | 298 | **299** | 296 |

Source index 3174 is the only CheMatic gain in each format, with no lost
identity. It is an exocyclic alkene attached to a puckered five-membered
ring. The old reader inferred E/Z from the XY projection even though the
record has genuine 3D coordinates; the projection reversed its configuration.
Unmarked, nonplanar 3D CTABs now use their full conformer for E/Z. A record
declared 2D, carrying explicit bond directions, or marking an alkene as
unspecified retains the existing interpretation. The exact source V2000
and V3000 records are regression fixtures.

Source index 591 remains an **Open Babel-only** exact identity in both
formats: CheMatic omits the tetrahedral tags in a substituted polycyclic
cage. Thus CheMatic's pass set is not a superset on this selection, despite
its higher total. Four other identities match only in CheMatic. Both readers
retain the mapless molecular graph for all 300 rows.

The candidate was also rerun on the earlier 740 development identities in
each format: CheMatic retained 740/740 and Open Babel 733/740. This is a
regression check on exposed data, not a new holdout.

## Evidence and reproduction

The [row-level report](results/unmarked-3d-chembl-disjoint-300-2026-10-17.json.gz)
contains the selection manifest, every source and normalized reader identity,
before/after exact flags, tool versions, binary and script SHA-256 values, and
the exact input hashes. Its gzip SHA-256 is
`a5eaaa7a75fdfbfffee6438b336320dd6788599d535c4bb2bde4d7ce22ff296b`.
The exact inputs are [V2000](results/unmarked-3d-chembl-disjoint-300-2026-10-17-v2000.sdf.gz)
(gzip SHA-256 `bb3aab6940cc92423738f746498bc9534ba21718e0dbed74129191641320d066`)
and [V3000](results/unmarked-3d-chembl-disjoint-300-2026-10-17-v3000.sdf.gz)
(gzip SHA-256 `06fb26dc7f775025d79af4556d210c09182421f890624b31550537f2f8f8e3e0`).
The ChEMBL source corpus SHA-256 is
`1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`.
The baseline and candidate example binaries were built with
`CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline -j1
-p chematic-mol --example sdf_v3000_accuracy_dump`; their SHA-256 values are
`cc4c32d4c9639f8fe31dd37e489c14e40c1af53e96889ace40e232c482cfb368`
and `e8a5aa569bec61fbfc02ac23449f39695e38d5f3e4ab1ebf3e477478367b4c8c`.

Extract the manifest from the report, decompress each SDF, and run
`validation/openbabel_comparison/check_unmarked_3d_stereo.py --case MANIFEST SDF`
with a source-built CheMatic example. The report also records the generator
and comparator hashes and the Open Babel binary hash.

The `chematic-mol` library's 805 tests, Clippy, format check, and the crate's
integration tests except one known main-line failure passed. The excluded
`mol2_3d_single_methyl_ring_remains_untagged` failure is reproduced on
unmodified main and is addressed separately in PR #986. The new V2000 and
V3000 exocyclic fixtures fail before this change and pass after it.

This measures source-relative identity under an RDKit gate, not independent
chemical truth or toolkit-wide superiority. These are source builds, not
published-package measurements.
