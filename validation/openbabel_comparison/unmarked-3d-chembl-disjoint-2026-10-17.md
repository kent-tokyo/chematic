# Open Babel 3.2.1: unmarked 3D alkene and cage stereo

## Paired result

RDKit 2026.09.1 embedded 150 tetrahedral and 150 alkene-stereo ChEMBL
identities, removed CTAB stereo marks, and admitted each row only when its
own read-back of both V2000 and V3000 retained the source's canonical
isomeric SMILES. The 300 identities are disjoint, under the same RDKit
canonicalization, from the 740 identities in the three earlier browser and
ChEMBL 3D selections. CheMatic before and after this change and Open Babel
3.2.1 read the **same** SDF bytes. RDKit normalized each output for exact
mapless graph and stereo comparison.

| Input | CheMatic main `db0cf128` | Alkene fix `7b63f432` | Alkene + cage `bde945d0` | Open Babel |
| --- | ---: | ---: | ---: | ---: |
| V2000, exact identity / 300 | 298 | 299 | **300** | 296 |
| V3000, exact identity / 300 | 298 | 299 | **300** | 296 |

An additional standard InChIKey comparison of the stored source and reader
SMILES gives 298/300 for the baseline, 299/300 after the alkene fix, and
**300/300 for both the final candidate and Open Babel** in each format. These
keys were generated through RDKit 2026.09.1's InChI wrapper, not by an
independent SDF reader. The four source-SMILES mismatches for Open Babel
collapse to the same standard InChIKey. Thus the exact-SMILES advantage above
does not establish an advantage under this second identity criterion.

Source index 3174 is the alkene fix's only gain in each format, with no lost
identity. It is an exocyclic alkene attached to a puckered five-membered
ring. The old reader inferred E/Z from the XY projection even though the
record has genuine 3D coordinates; the projection reversed its configuration.
Unmarked, nonplanar 3D CTABs now use their full conformer for E/Z. A record
declared 2D, carrying explicit bond directions, or marking an alkene as
unspecified retains the existing interpretation. The exact source V2000
and V3000 records are regression fixtures.

Source index 591 was an **Open Babel-only** exact identity after the alkene
fix: CheMatic omitted the tetrahedral tags in a monosubstituted adamantane
cage. The final reader recognizes only the exact subdivided K4 ring topology
with three hydrogen-bearing bridgeheads and one substituted bridgehead. Its
3D geometry then recovers all four source tags. Unsubstituted adamantane
remains untagged. The exact V2000 and V3000 source records are regression
fixtures. Four identities match the source SMILES only in CheMatic, while
both readers retain the mapless graph for all 300 rows.

The final candidate was also rerun on the earlier 740 development identities in
each format: CheMatic retained 740/740 and Open Babel 733/740. This is a
regression check on exposed data, not a new holdout.

## Evidence and reproduction

The [initial row-level report](results/unmarked-3d-chembl-disjoint-300-2026-10-17.json.gz)
contains the selection manifest, every source and normalized reader identity,
before/after exact flags, tool versions, binary and script SHA-256 values, and
the exact input hashes. Its gzip SHA-256 is
`a5eaaa7a75fdfbfffee6438b336320dd6788599d535c4bb2bde4d7ce22ff296b`.
The [final cage follow-up](results/unmarked-3d-chembl-disjoint-300-cage-2026-10-17.json.gz)
contains every final candidate and Open Babel identity, exact and standard
InChIKey flags, the source report and SDF hashes, and source/binary provenance.
Its gzip SHA-256 is
`fa577c95d6685313230aef6becdaeb32d2fbf2cd605493d3a75aef3f8b630fa6`.
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
The final candidate was built with the same command at `bde945d0`; its
example binary SHA-256 is
`e5deb74325bd664f7d8089a29aee8cde97c4e6a91bdf0d46f4a5445f4c565bd6`.

To reproduce the InChIKey counts from the row-level report, use RDKit with
InChI support:

```python
import gzip, json
from rdkit import Chem

with gzip.open("validation/openbabel_comparison/results/unmarked-3d-chembl-disjoint-300-2026-10-17.json.gz", "rt") as source:
    report = json.load(source)
key = lambda smiles: Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
for fmt, result in report["formats"].items():
    print(fmt, {engine: sum(key(row["expected"]) == key(row["actual"][engine])
                        for row in result["rows"])
                for engine in ("baseline", "candidate", "openbabel")})
with gzip.open("validation/openbabel_comparison/results/unmarked-3d-chembl-disjoint-300-cage-2026-10-17.json.gz", "rt") as source:
    final = json.load(source)
for fmt, result in final["formats"].items():
    print(fmt, "final", sum(key(row["expected"]) == key(row["candidate"])
                             for row in result["rows"]))
```

Extract the manifest from the report, decompress each SDF, and run
`validation/openbabel_comparison/check_unmarked_3d_stereo.py --case MANIFEST SDF`
with a source-built CheMatic example. The report also records the generator
and comparator hashes and the Open Babel binary hash.

The final candidate plus the independent ring-center repair from PR #986
passed all 807 `chematic-mol` library tests, integration suites, and doc tests.
Clippy, formatting, and the refactor budget passed. The PR branch by itself
still hits the pre-existing main failure
`mol2_3d_single_methyl_ring_remains_untagged`, addressed by PR #986. The new
V2000/V3000 exocyclic and cage fixtures fail before their respective fixes
and pass after them.

This measures source-relative identity under an RDKit gate, not independent
chemical truth or toolkit-wide superiority. These are source builds, not
published-package measurements.
