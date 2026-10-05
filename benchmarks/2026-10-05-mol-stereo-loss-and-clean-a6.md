# MOL stereo loss, write_mol speed, A6 on a clean commit

Follow-up to [the MOL stereo / CIP / MMFF record](2026-10-05-mol-stereo-cip-mmff-followups.md),
whose MOL writer numbers, `write_mol` timing and A6 row this record replaces.
Source branch only, one Linux x86-64 host, RDKit 2026.03.6. Code revision:
`b014d970` (the evidence commit that follows changes documents and result
files only).

## MOL writer: what is lost, and saying so

RDKit `MolFromMolBlock` on chematic's blocks (no coordinates given),
compared with the input on RDKit's canonical isomeric SMILES.
"Lost" means the block carries less than the input; no row in either
revision reads back with an inverted centre or a wrong E/Z.

| Corpus (stereo rows) | Read back as input | Centre lost | E/Z written "either" | Previous draft (`ea09a51a`) |
|---|---|---|---|---|
| exposed 10k (1,687) | **1,681** | 6 | 0 | 1,581 (37 centre rows, 81 E/Z rows) |
| ChEMBL 5k (1,670) | **1,663** | 6 | 1 | 1,593 (58 centre rows, 26 E/Z rows, 5 pseudo-asymmetric) |

V2000 and V3000 give the same RDKit SMILES on every row. chematic's own
reader agrees: after aromatizing both sides, only the same rows differ
(10k: 6, the other 43 differences are aromaticity spellings RDKit
reads as one molecule; 5k: 7).
[Summary](../validation/results/mol-writer-stereo-readback-source-2026-10-05.json).

What changed:

- **Macrocycle E/Z** (all 117 lost E/Z bonds on 10k, 49 of 50 on 5k, rings
  of 16–29 atoms): a ring drawn as a regular polygon makes every ring double
  bond cis. One end of the double bond, with its own substituents, is
  reflected across the line from its other ring neighbour to the far end
  (both bond lengths kept), when that sets the bond, changes no bond already
  drawn as declared and lands on no atom.
- **2D layout bug:** the fused-ring loop's bound was computed from the
  shrinking list of unplaced rings, so the last ring of pentacyclic
  triterpenes fell to an origin-centred fallback (bonds of 6–8 Å). The bound
  is now fixed at the start. A ring fused on two atoms goes to the side
  where it lands on fewer placed atoms. Molecules with any bond outside
  0.9–1.1 bond lengths: 233 / 10,000 (v1.0.34: 272), ChEMBL 5k 242 / 5,000
  (297). SVG depictions of such molecules change.
- **Bridges:** a centre with two bonds in nearly one direction gets a second
  layout (atoms in reverse order; the better drawing is kept), and a
  two-bond neighbour (a bridge atom) can be moved into the widest gap.
- **Wedges between centres:** a wedge also shows at its wide end, so a wedge
  that ends on a centre that already has its wedge is kept only if that
  centre still reads as declared. The "next to a double bond" penalty in
  wedge choice applies only to stereo double bonds (a carboxyl C=O pushed
  wedges onto bonds between two centres).

Remaining: one cage centre in each of four trichothecene macrolides
(roridin/verrucarin type, in both corpora), a bridged lactone and a
quinuclidine spiro-oxazoline (10k), an aphidicolane and an ellagitannin
(5k); one E/Z bond of a bryostatin (5k).

**Reporting.** These losses were detected but not exposed. Now:

| Binding | Report | Refuse |
|---|---|---|
| Rust | `write_mol_with_stereo_report`, `write_mol_v3000_with_stereo_report` → `MolStereoLoss` | caller checks `is_empty()` |
| Python | `Mol.to_mol_block_with_report()` → `(block, dict)` | `strict=True` on `to_mol_block`, `to_mol_block_2d`, `to_mol_v3000` (`ValueError`) |
| WASM | `mol_block_stereo_loss_json` | `to_mol_block_strict` |

`MolStereoLoss` lists unwedged centres, E/Z bonds written "either",
square-planar centres (no 2D field) and whether V2000 dropped enhanced
stereo groups. The default writers still return the lossy block; whether
they should refuse is open. WASM `to_mol_block` passed plain layout
coordinates to the writer, which then kept them as given and did not set
E/Z; it now uses the writer's stereo layout.

## write_mol speed

Release builds, exposed 10k (10,000 molecules), best of five calls per run,
median of five interleaved runs:

| | v1.0.34 | Previous draft | Source |
|---|---|---|---|
| `write_mol`, all rows | 8.2 ms | 190.1 ms | **77.9 ms** |

The source is about **9.5x slower** than v1.0.34, which wrote zero
coordinates and lost all stereo. Of the 77.9 ms: layout of the 1,325
molecules with stereo ≈ 21 ms, wedge choice and E/Z ≈ 21 ms, the 8,675
molecules without stereo ≈ 20 ms (mostly Kekulé bonds chosen in canonical
atom order where type 4 would lose an aromatic H, which v1.0.34 did not
write), the rest formatting. Removed from the draft: a molecule rebuild for
every wedge check (wedges are now toggled on one unmarked copy), recomputed
angles in a sort comparator, topological classes computed for every
molecule, the general float formatter (a fixed-point formatter with a
formatter fallback near ties; a unit test compares 400,000 values), and a
molecule clone for the aromatic-H check. Output is byte-identical to the
draft for the optimisations alone (checked on both corpora before the
loss-reduction changes). `write_mol_with_coords` with given coordinates
skips the layout.

## A6 on a clean commit

The 265 / 265 row in the previous record came from a wheel built at
`4823938f` with an uncommitted diff (sha256 `43b6d248…`), before later MMFF
changes. Rerun from a clean worktree of `b014d970` (`git diff` empty, sha256
of the empty diff `e3b0c442…`), CPython 3.11 release wheel `ca04127d…`:

| | Successes | Gradient-converged | Externally sound | Stereo-clean | Gross clash |
|---|---|---|---|---|---|
| chematic stereo-safe MMFF94 | 265 / 265 | 102 | 265 | 265 | 0 |
| RDKit ETKDGv3 + MMFF94 (pinned rows) | 264 / 264 | — | 264 | 264 | 0 |

The scorer (`pipeline_v2_vs_rdkit_common_scorer`) was built from the same
commit. [Summary](../validation/results/a6-clean-source-external-scorer-2026-10-05.json),
[rows](../validation/results/a6-clean-source-external-scorer-rows-2026-10-05.jsonl).
macOS and Windows were not run here; CI now has a Windows source-wheel job
(pytest, reaction 83 rows, SMARTS 310k, A6) beside the existing macOS one.
WASM could not be built on this host (the wasm32 target download is
blocked); `test-wasm` in CI builds it and runs the Node tests.
