# Stereo integrity, SMARTS, A6 and CIP follow-ups (after v1.0.34)

Source-branch results for ten open roadmap items, measured on one Linux
x86-64 host against RDKit 2026.03.6. Published-package rows are labelled as
such; everything else is unreleased source.

## Silent stereo and chemistry changes found by an order-invariance audit

Every Molecule-rebuilding function was run on the 3,357 stereo rows of the
exposed 10k corpus and `scripts/chembl_accuracy_corpus_4999.smi`, once on the
input and once on each of two random-SMILES spellings; a function whose
canonical output depends on input atom order, or that drops stereo where it
should not, is a finding.

| Function | Before | After | Cause |
|---|---|---|---|
| `add_hydrogens` → `remove_hydrogens` | 652 / 10,000 canonical changes (45 after the earlier H-count fix) | 0 / 10,000 and 0 / 5,000 | `add_hydrogens` dropped the bond-direction side table (aromatic-ring `/` of exocyclic imines), stereo groups and R-group labels |
| `reionize` | 1,956 / 6,714 order-variant | 0 | an ether/ester O was charged with two bonds, picking the first C neighbour |
| `uncharge` | radicals from nitro / N-oxide | 0 changes | forced every charge to 0 |
| `canonical_tautomer` (and `standardize`) | 105 / 6,714 order-variant; 216 / 10,000 rows turned an amide/urea/carbamate into an iminol | 0 order-variant; 0 / 10,000 and 0 / 5,000 iminols | a methoxy group sent the molecule through the aromatic O/C selection, whose O-H preference picked iminols; tetrahedral tags were kept on atoms that became sp2 |
| MOL V2000/V3000 read | benzene read as `C:1:C:C:C:C:C:1` | aromatic atoms | type-4 bonds did not flag their atoms aromatic |
| MOL V2000/V3000 write | `c1cc[nH]c1` read back as `c1ccnc1` | Kekulé bonds when type 4 would lose an aromatic H | |

`uncharge` now follows RDKit's `Uncharger` (cations keep their charge when
they have no H to lose; charge-separated groups stay; balancing anions are
chosen in canonical order). `reionize` deprotonates only O–H of carboxylic
acids and phenols and protonates only aliphatic primary/secondary amines.
`CanonicalMode::Backbone` keeps the old force-neutral topology key.

**Still open:** `write_mol` without coordinates writes no tetrahedral stereo
(there are no wedges to derive it from), and writes E/Z as wedge codes that
only chematic reads back. A coordinate-and-wedge writer is a separate task.
Kekulé fallback output depends on atom order (123 / 6,714 spellings), as any
Kekulé SMILES does.

## SMARTS and SMIRKS

* `^n` now follows RDKit's hybridization model (orbital count = total degree
  + lone pairs; four-orbital atoms with a conjugated bond and fewer than four
  neighbours are SP2; RDKit's conjugation candidates and transition-metal
  electron counts). Per-atom agreement with RDKit on the exposed 10k corpus:
  **219,883 / 220,015** (was 199,617; 116 of the 132 left are Kekulé-written
  thiophene/thiazole S, which RDKit perceives as aromatic and chematic keeps
  as written); ChEMBL 5k: 138,654 / 138,655.
* `k` alone means "in a ring" and `k0` matches every atom, as in RDKit.
* Bracket atoms read two-digit charges and H counts (`[C+10]`, `[CH10]`), so
  product templates `[C;+10:1]` / `[C;H15:1]` give RDKit's result (product
  kept / dropped) instead of a parse error.
* An aromatic N cation already at valence 4 kekulizes (`C[n+]1(C)cccc1`);
  a template that rewrites an aromatic ring bond as `=` between aromatic atoms
  gives no product, as RDKit's sanitize does for pyridine.
* SMARTS-query reactant templates on aromatic and Kekulé inputs: 244 / 247
  template × molecule cases give RDKit's sanitized products; the other three
  rewrite an aromatic ring bond with `~` → `=` (RDKit keeps two oddly kekulized
  products chematic drops).
* xsmarts-autoconf tweak run (3,000 examples, seed 1): 2 differing kinds
  (was 16 on v1.0.34, 4 before this batch); both are decided policies.

## Independent-corpus SMARTS rerun

[`scripts/check_smarts_independent_corpus.py`](../scripts/check_smarts_independent_corpus.py)
ran the 31 pinned queries on the 4,625 rows of
`scripts/chembl_accuracy_corpus_4999.smi` that are not in the exposed 10k
corpus (143,375 cells), against RDKit 2026.03.6 `GetSubstructMatches`.

| Artifact | Opt-in `find_matches_rdkit_parity` | Native `smarts_find` |
|---|---|---|
| PyPI v1.0.34 Linux wheel (CPython 3.9) | 143,375 exact | 143,213 exact, 162 differ |
| Source branch | 143,375 exact | 143,213 exact, 162 differ |

The 162 native differences are `[R1]`/`[R2]`/`[R3]` (161) and one `[k5]`,
all the documented SSSR ring basis. Reports:
[published](../validation/results/v1.0.34-published-pypi-linux-cp39-smarts-independent-chembl-4625.json),
[source](../validation/results/source-smarts-independent-chembl-4625.json).

## The 18 opt-in SMARTS refusals

[`scripts/check_rdkit_ring_order_dependence.py`](../scripts/check_rdkit_ring_order_dependence.py)
renumbered every exposed-10k molecule with three or more rings (4,328) six
times and re-ran eleven ring queries in RDKit. Exactly six rows and three
queries change their RDKit match set with atom order: rows 9, 23, 28, 29, 30
and 34 × `[R1]`/`[R2]`/`[R3]` — the 18 cells the opt-in profile refuses as
`ring_model_ambiguous`. Over 40 orders each row gives 10–32 distinct RDKit
answers. There is no order-independent RDKit result to match, so the typed
refusal is the decided outcome. [Report](../validation/results/rdkit-ring-order-dependence-10k.json).

## A6: Linux stereo failures (#739)

Force-field and distance-geometry code (chematic-3d, chematic-ff) now takes
`sin`, `cos`, `acos`, `atan2`, `exp`, `ln`, `pow`, … from the pure-Rust
`libm` crate instead of the host C library. Host libm results differ in the
last bit between glibc and Apple's libm, and the minimiser amplified that into
different trajectories; `wasm32` builds already used `libm`, so native and
WASM now compute the same values. On Linux the stereo-safe MMFF94 arm of the
265-row A6 set now returns **265 / 265** successes (published v1.0.31/v1.0.34
Linux wheels: 263, failing rows 53 and 246 with `FinalStereoViolation`;
macOS: 265): 102 gradient-converged, 162 iteration-limit, 1
constraint-rejected fallback (output SHA-256 `9a9afe19…`). macOS was not
rerun here, and the external geometry/stereo/clash scorer was not rerun.

## CIP stability

On the 1,687 stereo rows of the exposed 10k corpus, Accurate-mode
tetrahedral labels, the five typed abstentions (one lone-pair centre, four
phosphorus) and bond-keyed E/Z labels are identical across five random
SMILES spellings, an `add_hydrogens`/`remove_hydrogens` round trip and a
canonical-SMILES reparse (1,687 / 1,687 rows). The atom-keyed E/Z entries of
`assign_cip` sit on the double bond's first atom by index, as documented;
`assign_ez_bonds` is the order-independent form.

These are measurements of the named artifacts on one host; they do not
claim general RDKit parity.
