# COSMolKit parity record — 2026-10-08

Source: commit `a8ca103c` from branch `claude/cosmolkit-parity`, later merged
into the v1.0.38 release candidate.
Reference: RDKit 2026.03.1. Competitor: COSMolKit 0.3.0 (PyPI wheel).
Harness: `validation/cosmolkit_comparison/` (`run_corpus.py`,
`compare_corpus.py`). A row counts as a match only when the value equals
RDKit's exactly (floats compared bit-for-bit; ints and floats compared as
numbers); refusals, unsupported operations and errors never count.
Bit-identical statements below are limited to that recorded comparison lane;
macOS arm64 can follow a slightly different minimizer path, so portable unit
tests use a `1e-9` absolute tolerance for optimized force-field energies and
`1e-6` for coordinates (`1e-12` for the same-coordinate initial energy).

Corpora:

- **ChEMBL 5k** — `scripts/chembl_accuracy_corpus_4999.smi` (5,000 rows)
- **RDKit.js 10k** — `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`
  (10,000 rows)

## Corpus harness results (chematic wheel built from this branch)

Operations where both toolkits match RDKit on every row of both corpora:
aromatic ring count, canonical SMILES and its round trip, Chi0/0v/1/1v/2v/3v/4v,
CIP labels, exact and average mass, formula, atom-pair / Avalon / layered /
pattern / RDKit / torsion fingerprints, Morgan r2 and r3 (2048 bits), Morgan
with `includeChirality`, Morgan `countSimulation`, Fsp3, Hall–Kier alpha,
HBA, HBD, heavy atoms, InChI, InChIKey, InChI read-back (`MolFromInchi`),
Kappa1–3, Labute ASA, Crippen logP and MR, MACCS, MQN, ring counts (all
aliphatic/aromatic/saturated/hetero variants), amide, bridgehead, spiro and
heteroatom counts, QED, rotatable bonds, SlogP_VSA, SMR_VSA, stereocentre
counts, TPSA, SMILES writer options (Kekulé, non-isomeric, all bonds and Hs
explicit), SMILES after `AddHs`, and the 31-pattern SMARTS gate
(155,000/155,000 and 310,000/310,000 pattern-row pairs).

Operations where the two differ:

| Operation | Corpus | COSMolKit 0.3.0 | chematic |
|---|---|---:|---:|
| BalabanJ, BertzCT, PEOE_VSA, hashed atom-pair counts | both | unsupported | 100% |
| Ipc | ChEMBL 5k | unsupported | 4967/5000 |
| Ipc | RDKit.js 10k | unsupported | 9974/10000 |
| MOL block round trip | ChEMBL 5k | 4993/5000 | 5000/5000 |
| MOL block round trip | RDKit.js 10k | 9992/10000 | 10000/10000 |
| Murcko scaffold | ChEMBL 5k | 2337/5000 | 4998/5000 |
| Murcko scaffold | RDKit.js 10k | 6329/10000 | 9932/10000 |
| EnumerateStereoisomers (defaults) | ChEMBL 5k | 5000/5000 | 4997/5000 → 5000/5000* |
| EnumerateStereoisomers (defaults) | RDKit.js 10k | 10000/10000 | 9972/10000 → 9993/10000* |

\* The harness run predates RDKit's seeded random sampling for molecules with
more than 1024 isomers; the sampling port was then checked with a separate
script on the affected rows (3/3 and 21/22). Remaining stereoisomer
differences: one double-bond reference-atom choice that needs RDKit's
`FindPotentialStereo` ranking, five aromatic `[s+]([O-])` rows and one
ferrocene row.

Ipc differences come from RDKit's eigenvalue path (BLAS-dependent last-bit
results), not from the graph.

Murcko differences: RDKit's `MurckoDecompose` recomputes implicit hydrogens on
hypervalent atoms that lost substituents (P, S+, Sb, Ti, Ni), and one
bridgehead-N stereo label differs.

## Other RDKit-exact surfaces checked outside the harness

| Surface | Scope | Result |
|---|---|---|
| `Compute2DCoords` coordinates (`Mol.rdkit_2d_coords`) | both corpora | 5000/5000, 10000/10000 bit-identical on the recorded lane; portable tests use `1e-12` coordinate tolerance |
| `MolToMolBlock` after `Compute2DCoords` (`Mol.rdkit_mol_block_2d`) | both corpora | 5000/5000, 9999/10000 (V3000 output not written yet) |
| MMFF94 / MMFF94s energy, gradient, `MMFFOptimizeMolecule(maxIters=200)` | first 1000 rows of each corpus (990 + 994 compared) | all bit-identical on the recorded lane |
| UFF energy, gradient, `UFFOptimizeMolecule(maxIters=200)` | same | all bit-identical on the recorded lane |
| `MolToSmiles` option matrix (32 switch combinations × 3 roots) | both corpora | 1,440,000 / 1,440,000 |

COSMolKit 0.3.0 also reproduces `Compute2DCoords` (40/40 on a pilot) and
`MMFFOptimizeMolecule` (99/99 on a pilot).

## Known gaps where COSMolKit is ahead

- **Seeded ETKDG embedding.** COSMolKit reproduces
  `EmbedMolecule(AddHs(m), randomSeed=42)` coordinates exactly (pilot 99/99).
  chematic's RDKit embedding port is in progress (distance-geometry core exact
  on a 60-molecule pilot with RDKit's bounds matrix; bounds-matrix builder and
  ETKDG torsion terms not yet ported) and is not part of this branch.
- **SMARTS writer** (`MolToSmarts` / `MolToCXSmarts`), **PDB writer**
  (`MolToPDBBlock`), `GetStereoisomerCount` and `FindMolChiralCenters`
  equivalents are not yet offered with RDKit-exact output.

## Speed

Not re-measured after these correctness changes. The last measurement
(`bench_corpus.py`, ChEMBL 5k, median of 5 rotating blocks, before the parity
work) had chematic faster than COSMolKit on every timed operation except the
SMARTS gate (0.6×). Re-run `bench_corpus.py` before quoting speed for a
release.
