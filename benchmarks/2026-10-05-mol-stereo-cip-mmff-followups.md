# MOL stereo writer, CIP phosphorus and lone pairs, MMFF typing (after v1.0.34)

Source-branch results for ten open roadmap items, measured on one Linux
x86-64 host against RDKit 2026.03.6. Published-package rows are labelled as
such; everything else is unreleased source. Corpora: the exposed 10k lane
(`validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`, 1,687 stereo
rows) and `scripts/chembl_accuracy_corpus_4999.smi` ("ChEMBL 5k", 1,670
stereo rows).

## MOL/SDF stereo

> Superseded for the writer by
> [2026-10-05-mol-stereo-loss-and-clean-a6.md](2026-10-05-mol-stereo-loss-and-clean-a6.md)
> (1,681 / 1,687 and 1,663 / 1,670, with a loss report). The numbers below
> are from `ea09a51a`.

**Writer.** `write_mol` / `write_mol_v3000` now lay a molecule with stereo
out (`chematic_mol::stereo_depiction`), reflect one side of each E/Z bond
whose drawn geometry disagrees with its declared markers, and draw one wedge
or hash per tetrahedral centre from the centre, kept only when re-reading it
gives the declared parity. A centre is left undrawn rather than drawn
ambiguously: two of its bonds closer than 15°, or a four-bond centre whose
wedge sits between plain bonds that span less than a half-plane. A one-bond
neighbour of a centre (methyl, OH, halogen) is moved into the widest gap
first. Ring E/Z bonds that reflection cannot set, and stereogenic double
bonds without declared E/Z, are written "either" (V2000 stereo 3, V3000
`CFG=2`).

RDKit `MolFromMolBlock` on the written blocks (no coordinates given):

| Corpus | Read back as input | E/Z written "either" | Tetrahedral centre lost | Wrong label |
|---|---|---|---|---|
| exposed 10k, 1,687 stereo rows | **1,581** | 69 | 37 | 0 |
| ChEMBL 5k, 1,670 stereo rows | **1,593** | 14 | 63 | 0 |
| v1.0.34 `to_mol_block`, 10k | 0 (170 unreadable) | | | |

Lost centres sit in bridged or cage systems where the layout puts atoms on
top of each other or bonds in one direction. The coordinates path also wrote
the atom-map column one character late (RDKit read every atom as mapped 0).

**Reader.** RDKit-written blocks (`Compute2DCoords`) read by `parse_mol`,
compared on RDKit's canonical SMILES:

| Corpus | v1.0.34 | Source |
|---|---|---|
| exposed 10k stereo rows | 1,548 / 1,687 | **1,666 / 1,687** |
| ChEMBL 5k stereo rows | 1,596 / 1,670 | **1,647 / 1,670** |

Changes: two double bonds whose carriers disagree on a shared bond are
reconciled by flipping one (only an inconsistent cycle rejects), and an end
whose two substituents both carry markers needs them on opposite sides; a
wedge bond can carry an adjacent E/Z direction in `bond_direction`, which
every E/Z reader now prefers over the literal order; double bonds in rings
of fewer than eight atoms get no E/Z (RDKit's `shouldDetectDoubleBondStereo`).
Remaining differences: one wedge reading on an artemisinin-type cage, RDKit
reading its own unspecified C=N geometry as E, and Kekulé exocyclic imines on
a four-membered ring.

**Canonical SMILES.** A marker that specifies no E/Z is dropped before
canonicalization, as RDKit drops it on parse, and ends of non-stereogenic
double bonds (`=C(C#N)2`, `=CD2`) no longer steer carrier choice. On 15,000
rows: 0 unstable, 0 order-variant (two random spellings), RDKit reads every
output as the input; 8 strings differ from v1.0.34, all carrier moves.

## Hybridization

`hybridization_per_atom` uses `chematic_smarts::rdkit_hybridization` on the
RDKit aromaticity view. Agreement with RDKit `GetHybridization`:

| Corpus | v1.0.34 | Source |
|---|---|---|
| exposed 10k atoms | 199,617 / 220,015 | **220,003 / 220,015** |
| ChEMBL 5k atoms | 126,658 / 138,655 | **138,654 / 138,655** |

SMARTS `^n` through `find_matches_perceived` (the Python `smarts_find` path)
on Kekulé-written copies of both corpora, nine queries: 134,990 / 135,000
cells exact; the 10 differ on Zn/Hg dithiocarbamate complexes.

## MMFF94 atom types

`scripts/mmff94_atom_type_census.py` (heavy atoms typed differently from
RDKit `MMFFGetMMFFAtomType`):

| Corpus | Before this batch | Source |
|---|---|---|
| exposed 10k | 451 | **89** (78 in one fullerene cage, row 3837) |
| ChEMBL 5k | — | **103** (cyclophosphazenes, porphyrin, indolocarbazole) |

Ported from RDKit's `setMMFFHeavyAtomType`: aliphatic carbon by total degree
(CR4E 30, CO2M 41, CNN+ 57; an MMFF-aromatic seven-ring carbon was typed 1),
nitroso N 46, NSO N 48, N-oxides 67/68, N=N-N / N=N-C N 10, =N= 53,
isonitrile N 61, phosphorus by degree, O 49/51/70, perchlorate Cl 77.
[10k](../validation/results/mmff94-atom-type-census-source-10k-2026-10-05.json),
[5k](../validation/results/mmff94-atom-type-census-source-chembl5k-2026-10-05.json).

## CIP

RDKit's CIPLabeler does not expand bond orders at the root, so a P=O oxygen
ranks like the O⁻ of `[P+][O-]`; the accurate engine now does the same.
Acyclic phosphorus is labelled (the four phosphonamidate rows of the 10k lane
agree with RDKit, which is stable over 30 atom orders and the charge-separated
spelling); P on an unsaturated ring stays `OracleUnstable` (cyclophosphazene
Kekulé respellings flip RDKit's label). A lone pair outside rings is a
lowest-priority phantom ligand: 540 / 540 sulfoxide, sulfinamide,
sulfoximine, phosphine and selenoxide centres agree with RDKit over random
atom orders. Ring lone-pair centres (bridgehead amines) stay
`LonePairCenter`: RDKit's label depends on how a ring-closure digit on the
centre was written, which `stereo_neighbor_order` does not keep.

| Corpus | Agree with RDKit | Abstain | Differ |
|---|---|---|---|
| exposed 10k stereo rows | 4,394 | 1 (was 5) | 0 |
| ChEMBL 5k stereo rows | 4,166 | 20 | 0 |

All 1,687 exposed stereo rows keep identical labels across five random
spellings, an H round trip and a canonical reparse.

## Native SMIRKS: centre gaining a bond

`[C;$(CO):1]>>[C:1]C(=O)C` on `C[C@H](N)O`: the native runner now keeps the
centre, with the new group in the H's place (geometric retention). RDKit
copies the raw tag onto a bond list with the new bond appended, which is the
opposite arrangement (`CC(=O)[C@](C)(N)O` versus native
`CC(=O)[C@@](C)(N)O`); `run_smirks_checked(..., rdkit_compat=True)` gives
RDKit's. Unit test: `centre_gaining_a_neighbour_in_place_of_its_h_keeps_stereo`.

## Published 310k opt-in SMARTS gate

`scripts/check_python_smarts_parity_310k.py --archived-oracle` checks the
pinned RDKit 2026.03.6 oracle by SHA-256 instead of importing RDKit, so it
runs on the PyPI v1.0.34 Linux wheel (CPython 3.9 only): **309,982 exact,
18 typed refusals (the order-dependent `[R1]`/`[R2]`/`[R3]` cells), 0
unexpected**. [Report](../validation/results/v1.0.34-published-pypi-linux-cp39-python-smarts-parity-310k.json).

## A6 external scorer on the libm build

> This row was measured on a wheel built from `4823938f` plus an
> uncommitted diff, not on a commit. The clean-commit rerun is in
> [2026-10-05-mol-stereo-loss-and-clean-a6.md](2026-10-05-mol-stereo-loss-and-clean-a6.md).

The stereo-safe MMFF94 arm rerun with this batch's source wheel (libm math
and the new MMFF typing): 265 / 265 successes, 102 gradient-converged.
`crates/chematic-3d/examples/pipeline_v2_vs_rdkit_common_scorer.rs` on those
rows: 265 / 265 scored, independently sound, stereo-clean, no gross-clash
rows (RDKit ETKDGv3 + MMFF94: 264 / 264).
[Summary](../validation/results/a6-libm-source-external-scorer-2026-10-05.json).
macOS was not rerun.

## Performance

Release builds of v1.0.34 and the branch on the same host, median of three
interleaved runs over the exposed 10k corpus:

| Operation | v1.0.34 | Branch | Ratio |
|---|---|---|---|
| canonical SMILES | 552.5 ms | 587.4 ms | 1.06 |
| add + remove H | 117.2 ms | ~126 ms | 1.08 |
| SMARTS `[C^2]` | 28.6 ms | ~44 ms | 1.5 |
| SMARTS `[R2]` | 101.6 ms | 102.3 ms | 1.01 |
| kekulize | 11.0 ms | 11.1 ms | 1.01 |
| `write_mol` | 11.2 ms | 190.1 ms | 17x slower (superseded: 9.5x, see the follow-up record) |
| MMFF94 typing | 573.6 ms | 549.4 ms | 0.96 |
| UFF 3D, 249 molecules | ~127 s | ~127 s | 1.0 |
| 13 reaction cases | 409 µs/iter | 437 µs/iter | 1.07 (Kekulé reactants aromatized) |

These are measurements of the named artifacts on one host; they do not
claim general RDKit parity.
