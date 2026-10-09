# COSMolKit parity record, round 2 — 2026-10-09

Source: branch `claude/cosmolkit-parity-2` (unreleased, on top of v1.0.38).
Reference: RDKit 2026.03.1. Competitor: COSMolKit 0.3.0 (PyPI wheel).
Previous record: [2026-10-08-cosmolkit-parity.md](2026-10-08-cosmolkit-parity.md).
Match rule as before: exact equality with RDKit; refusals, unsupported
operations and errors never count.

Corpora: ChEMBL 5k (`scripts/chembl_accuracy_corpus_4999.smi`), RDKit.js 10k
(`validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`), and a generated
non-tetrahedral corpus (9,472 SMILES with `@SP`/`@TB`/`@OH` centres; RDKit
reads 7,845).

## Corpus harness (`run_corpus.py` / `compare_corpus.py`, all ops)

Full run on a wheel built from this branch (before the alignment, reader,
non-tetrahedral and MolHash work, which were then checked op by op): chematic
matches RDKit on every row of both corpora for every operation except Ipc
(4967/5000, 9974/10000; RDKit's value is BLAS-dependent). COSMolKit matches
RDKit at least as often as chematic on no operation.

New operations in this round:

| Operation | ChEMBL 5k chematic | COSMolKit | RDKit.js 10k chematic | COSMolKit |
|---|---:|---:|---:|---:|
| `MolToSmarts` | 5000 | 4204 | 10000 | 9251 |
| `MolToCXSmarts` | 5000 | 4204 | 10000 | 9250 |
| `MolToPDBBlock` | 5000 | 363 | 10000 | 2042 |
| `FindMolChiralCenters(includeUnassigned=True)` | 5000 | 0 | 10000 | 0 |
| `GetStereoisomerCount` | 5000 | 5000 | 10000 | 10000 |
| `EnumerateStereoisomers` | 5000 | 5000 | 10000 | 10000 |
| Murcko scaffold | 5000 | 2337 | 10000 | 6329 |
| ExtendedMurcko (`rdMolHash`; COSMolKit `net_scaffold`) | 5000 | 4028 | 10000 | 8580 |
| Morgan bit-info map | 5000 | 5000 | 10000 | 10000 |

COSMolKit's `find_chiral_centers` returns every atom for every row, so it never
matches.

## Checked outside the harness

| Surface | Inputs | chematic vs RDKit |
|---|---|---|
| All 19 `rdMolHash` functions | ChEMBL 5k, RDKit.js 10k, non-tetrahedral 7,845 | every row (CXSMILES variant: every row matches or both raise) |
| `AlignMol`, `GetAlignmentTransform`, `GetBestRMS`, `CalcRMS` | heavy-atom conformers (seeds 7 vs 42), 966 ChEMBL + 950 RDKit.js rows | bit-identical |
| `MolFromPDBBlock` (heavy-atom and with-H blocks) | 1,977 ChEMBL, 1,486 RDKit.js, 60 peptides | all identical (SMILES, coordinates, failures) |
| `MolFromXYZBlock` | 1,995 blocks + 13 edge cases | all identical |
| `MolFromMol2Block` | 3,479 molecules × 3 Tripos typings, 127 edge cases, flattened / all-zero sets | all identical |
| Non-tetrahedral corpus: SMILES (all writer options), chiral centres, stereoisomers, Murcko, SMARTS/CXSMARTS, PDB | 7,845 | every row |
| Non-tetrahedral corpus: `Compute2DCoords` + `MolToMolBlock` | 7,845 | 7,725 (the other 120 are trigonal-bipyramidal centres where RDKit reads past its template array; output depends on process state) |

On the non-tetrahedral corpus COSMolKit reads 224 fewer molecules, aborts the
process on 112, and matches RDKit's SMILES and MOL block on 7,509.

## Speed (ChEMBL 5k, `bench_corpus.py`, median of 5 blocks, ms)

chematic is timed through its RDKit-exact APIs; the `molblock` operation lays
out 2D coordinates in all three engines.

| Operation | chematic | COSMolKit | RDKit |
|---|---:|---:|---:|
| parse | 36 | 913 | 810 |
| canonical SMILES | 411 | 526 | 512 |
| formula | 13 | 25 | 23 |
| mol_wt | 7 | 20 | 111 |
| TPSA | 39 | 116 | 170 |
| logP | 384 | 4602 | 1851 |
| QED | 2095 | 27568 | 9903 |
| Morgan r2 | 232 | 388 | 366 |
| MACCS | 3157 | 8512 | 4843 |
| CIP | 444 | 1771 | 514 |
| 31 SMARTS | 745 | 1297 | 1018 |
| MOL block (2D) | 1682 | 7145 | 3156 |
| InChI | 1231 | 5031 | 2471 |
| SMARTS writer | 337 | 335 | 411 |
| Murcko | 445 | 194 | 673 |
| stereoisomers | 2513 | 10009 | 9423 |

Parse, canonical SMILES, 31 SMARTS, SMARTS writer and Murcko come from the
last speed session (load about 0.8 from other jobs); the other rows from an
idle run earlier on the same branch, before those speed changes.
COSMolKit sanitizes at parse time, so its Murcko time excludes work
that chematic's Murcko time includes; parse plus Murcko is 481 ms for chematic
and 1,107 ms for COSMolKit.

## Seeded ETKDG embedding

`embed3d` harness op: `AllChem.EmbedMolecule(Chem.AddHs(m), randomSeed=42)`
coordinates (ETKDGv3 defaults) against COSMolKit
(`EmbedParameters.etkdg_v3()`, `random_seed=42`) and chematic
(`Mol.add_hydrogens().rdkit_embed(random_seed=42)`), exact float equality,
rows where RDKit's embedding succeeds:

| Corpus | COSMolKit | chematic |
|---|---:|---:|
| ChEMBL 5k | 4977/4977 | 4977/4977 |
| RDKit.js 10k | 9946/9947 | 9947/9947 |

Where RDKit's embedding fails (23 and 52 rows) chematic fails too; one
RDKit.js row raises inside RDKit (`bad lower bound`) and in chematic. The
same full run (wheel from `claude/rdkit-etkdg` 981f31ff, which includes this
branch) kept every other operation at its earlier result: 100% except Ipc.
Bit-identical coordinates are claimed for the Linux x86-64 lane recorded here;
other platforms' libm can change last bits, which the minimizations can
amplify.

## Remaining gap

None found on the compared surfaces: chematic matches RDKit at least as often
as COSMolKit on every operation in the harness.
