# RFC: closing the chematic-3d / RDKit ETKDG gap (diagnosis + planning only — no implementation in this round)

Agent: "Agent G" (diag/etkdg-3d-gap). Scope: `chematic-3d` (ETKDG-style conformer
generation, MD, SASA, WHIM) and `chematic-ff` (MMFF94 / DREIDING atom typing),
per `CLAUDE.md`'s crate table. This round is **diagnosis and planning only** —
no conformer-generation algorithm code was written or modified, and nothing
under `crates/*/src/**` was touched. Reference oracle: RDKit `2026.03.3`
ETKDGv3 (`AllChem.EmbedMolecule`/`EmbedMultipleConfs`), via the repo's
`.venv/bin/python`.

## Files touched by this round

- `docs/etkdg_3d_gap_rfc.md` — this document (new).
- `scripts/etkdg_vs_rdkit_gap.py` — re-runnable comparison harness (new).
- `validation/results/etkdg_vs_rdkit_summary.json` — aggregate measurements (new).
- `validation/results/etkdg_vs_rdkit_rows.jsonl` — per-molecule raw rows, nothing
  silently dropped (new).

**Out of scope / not touched**: anything under `crates/*/src/**` (no algorithm
code, not even a one-line fix — candidate fixes are called out below as
follow-ups for a future, separately-authorized round); `feat/io-mrv`,
`feat/io-tdt`, `feat/io-smiles-supplier-writer`, `fix/smiles-bracket-implicit-h`,
`diag/stereo-reader-integration-boundary`, `feat/stereo2d-local-parity`, or any
other agent's branch/files (2D stereo, aromaticity, descriptor census,
canonical-SMILES — separate concurrent diagnosis streams).

## Done condition for this round

A PR against `main` containing the four files above, with the measurements in
this document traceable to a `git diff`-visible script and a committed JSON/
JSONL artifact — not narrative claims. No src changes, no merge.

---

## 1. What the current code actually does (read end to end, not guessed)

`chematic-3d`'s "ETKDG" entry point (`generate_coords_etkdg` /
`generate_coords_etkdg_with_noise`, `crates/chematic-3d/src/etkdg.rs`) is:

```
dg::generate_coords(mol)                    // deterministic DFS bond/angle/dihedral placement
  -> apply_torsion_preferences_with_noise    // snap ~180 hand-coded/SMARTS torsion rules
  -> build_constraints + satisfy_constraints // 3 iterations of local bond/angle projection
  -> snap_amide_torsions                     // post-hoc amide-planarity fixup
```

`dg::generate_coords` (`crates/chematic-3d/src/dg.rs`) is explicitly documented
as "not a full distance-geometry solver" — it is a **deterministic DFS
bond-angle-dihedral builder**: rings go on regular-polygon/chair/envelope/crown
templates, chain atoms are placed by walking the bond graph with ideal lengths
and 120°-staggered dihedrals. There is **no distance bounds matrix, no
triangle-inequality smoothing, no random embedding, and no eigen/Gram-matrix
step** anywhere in this path — i.e. it is not distance geometry in the sense
RDKit's ETKDG (or the `dgeom`/`ChemBase` literature) uses the term.

A **second, complete, but entirely orphaned implementation of exactly what's
missing already exists** in the same crate:
`crates/chematic-3d/src/dg_fft.rs` (`#![allow(dead_code)]`, ~840 lines,
tested, has a real commit history across multiple sprints —
`4923818 feat: A3 — ETKDGv3 distance geometry with eigenvalue decomposition`).
It builds a bond/angle/VdW bound matrix, applies Floyd-Warshall
triangle-inequality smoothing (`smooth_bounds`), converts to a Gram matrix via
classical MDS double-centering, does a Jacobi eigendecomposition, and refines
with SHAKE-like bound-driven relaxation. `generate_coords_dg`, the module's
only public entry point, is declared `pub` in `lib.rs` but **has zero callers
anywhere in the workspace** (verified: `grep -rn "generate_coords_dg\b"`
across every crate returns only its own definition and its own tests). Phases
1-2 below are therefore "audit and wire up (or replace) existing dead code,"
not "build from scratch."

Chirality/parity is **never read** during coordinate generation: `grep -n
"chiral\|parity\|Chirality" crates/chematic-3d/src/{dg,etkdg,constraints,
dg_fft,etkdg_knowledge}.rs` returns zero matches. The declared `@`/`@@`/`/`/`\`
stereo descriptor from the input molecule has no path into `generate_coords`,
`apply_torsion_preferences_with_noise`, or `satisfy_constraints`. Phase 3 does
not exist in any form, not even a stub.

`Coords3D` (`crates/chematic-3d/src/coords.rs`) is explicitly "3D coordinates
for all **heavy** atoms" — there is no implicit-hydrogen materialization
anywhere in `chematic-3d` or `chematic-ff` (`grep -rn "implicit_hcount"` in
both crates: zero hits). This is a real, load-bearing scope difference from
RDKit, which embeds all-atom (`Chem.AddHs`) by default, and it is the direct
cause of the chirality-coverage gap measured in §3.

`chematic-ff` separately ships a **complete** MMFF94 implementation —
`mmff94_energy::{bond, angle, oop_stbn (stretch-bend + out-of-plane),
torsion, vdw}` plus its own L-BFGS minimizer (`mmff94_minimizer::
minimize_mmff94_lbfgs`, `minimize_mmff94_full`, `mmff94_energy_breakdown`) —
and a separate UFF module (`uff::{minimize_uff, uff_total_energy}`). None of
this is used by `chematic-3d`'s own minimizer
(`crates/chematic-3d/src/minimize.rs`), which hand-rolls a **much smaller**
energy function (bond + angle + VdW, +electrostatic for MMFF94 only — **no
torsion term, no stretch-bend, no out-of-plane term**) and a naive
fixed-step finite-difference gradient descent, and even `minimize_uff` in
`chematic-3d` is just an alias for the generic bond/angle/VdW minimizer, not a
call into `chematic_ff::uff`. This is confirmed by the Python surface too:
`Mol.mmff94_energy_breakdown()` already calls the *complete* chematic-ff
energy function directly and is used as the "fair yardstick" in §2 below —
but the conformer that gets *minimized* never sees the torsion/oop/stretch-bend
terms that same function would score it on. Phase 8 is therefore substantially
"wire the existing full implementation into the existing minimization
pipeline," not "implement MMFF94."

Ensemble generation, RMSD (Kabsch, `conformer.rs`), duplicate detection
(`is_duplicate`), greedy-leader clustering (`cluster_conformers_by_rms`), and
USR-based diversity (`conformer_diversity_usr`) all **already exist** and are
exercised by real tests — Phases 6-7 are about *fixing the geometry feeding
into them* (§2), not building the ensemble/pruning machinery itself.

The RNG (`crates/chematic-3d/src/prng.rs`) is a process-global atomic-counter
xorshift64 with **no public seed parameter** — ensemble runs are not
bit-reproducible across runs (noted in the JSON `meta.reproducibility_note`;
the `noise_sigma_deg = 0.0` single-conformer path used for every non-ensemble
metric below *is* deterministic, and was re-run twice during this diagnosis
with identical output).

---

## 2. Measured current state — Phase 0: an existing correctness defect, found before any of the roadmap phases were even reachable

Running `scripts/etkdg_vs_rdkit_gap.py` against a 58-molecule curated corpus
(rigid/fused rings, flexible chains, macrocycles, sp3 stereocenters split by
heavy-neighbor-count, alkene E/Z pairs, drug-like/steroid/stress cases — see
the `CORPUS` list in the script) surfaced something more urgent than the
missing-distance-geometry gap itself: **the pipeline that exists today
produces geometrically torn structures for most of the corpus.**

| metric | value |
|---|---|
| `raw_embed_returned_rate` (chematic) | 93.1% (RDKit: 100%) — the embedder "returns coordinates," but see next row |
| **`geometrically_valid_rate`** (chematic) | **25.9%** (15/58) — fraction with no bond >50% off an external covalent-radius reference |
| bond-length-blowup bucket | **39/58 (67%)** |
| gross-clash bucket (chematic) | 4/58 |

"Geometrically valid" here means *no bond exceeds 50% relative error* against
an external, engine-independent reference (RDKit's own periodic-table
covalent radii × a bond-order scale factor — **not** chematic's own internal
tables, so this isn't circular). Real force-field relaxation artifacts are a
few percent even for strained rings; 50%+ is unambiguously a torn molecule,
not imprecision. This bucket is reported separately (`status ==
"chematic_bond_length_blowup"`), not silently folded into "success," per the
task's no-silent-drops requirement — every one of the 58 rows in
`validation/results/etkdg_vs_rdkit_rows.jsonl` carries an explicit named
`status`.

**Isolation.** Bypassing the ETKDG-specific step (`Mol.generate_3d()`, which
is `dg::generate_coords` + `minimize_dreiding` with **no** call into
`etkdg.rs`) gives **0 blow-ups on the same corpus**. Going through the actual
ETKDG path (`Mol.conformer_ensemble(1, 0.0, force_field, 0.0)`, i.e.
`etkdg::generate_coords_etkdg_with_noise`) reproduces the blow-up under
**both** `"dreiding"` and `"mmff94"` force-field choices. This pins the defect
to `etkdg.rs`'s torsion-preference/constraint-repair stage, not `dg.rs`'s base
placement or the minimizer alone. Minimal repros:

```python
import chematic, numpy as np
m = chematic.from_smiles("CCCCCCCCCC")  # decane — no ring, no exotic atoms
def worst_bond(coords):
    a = np.array(coords)
    return max(np.linalg.norm(a[i] - a[i + 1]) for i in range(len(a) - 1))

worst_bond(m.generate_3d())                                   # ~1.58 Å — fine
worst_bond(m.conformer_ensemble(1, 0.0, "dreiding", 0.0)[0])  # ~11.3 Å — torn
```

Source inspection localizes (at least) three contributing mechanisms, found
by reading the code, not by fixing it — **no fix was attempted**, per this
round's diagnosis-only mandate:

1. **Under-iterated constraint repair.** `etkdg.rs` calls
   `satisfy_constraints(&coords, mol, &constraints, 3)` — exactly 3
   iterations — immediately after torsion preferences have potentially
   rotated an entire subtree by up to ~180°. `constraints.rs`'s own module
   doc says "Fast convergence: 5–10 iterations typical," and its projection
   is local/non-rigid (each bond or angle constraint moves only its own 2-3
   atoms, ignoring every other constraint those atoms participate in) — not a
   convergent global solver. 3 iterations undershoots the module's own stated
   minimum for exactly the class of large perturbation the torsion step just
   introduced.
2. **Ring-unaware torsion rotation.** `get_torsion_preference`
   (`etkdg_knowledge.rs`) classifies atoms purely by hybridization/aromaticity
   type, with no ring-membership check on the B–C bond (only the separate
   SMARTS-rule path consults `ring_bond_set`, and only to gate *which SMARTS
   rules* apply). When it fires on a bond that is part of a ring,
   `set_dihedral` → `find_rotated_atoms` (`mol_transforms.rs`) computes the
   "rotatable subtree" via a BFS that excludes only the literal parent atom
   index, not the (parent, target) edge or cycle membership — for a ring bond
   this walks nearly the entire ring back around to the far side of the fixed
   pivot, and rotating that while the pivot stays put tears the ring-closing
   bond on the other side. This is consistent with the blow-up set skewing
   heavily toward ring-containing molecules (naphthalene, thiophene,
   adamantane, cubane, indole, purine, quinoline, anthracene, pyrene, and
   every polycyclic/aromatic drug-like case in the corpus).
3. **Silent missing-parameter zero-gradient under MMFF94.**
   `bond_energy_mmff94`/`angle_energy_mmff94` (`chematic-3d/src/minimize.rs`)
   use `if let Some(params) = mmff94_bond_params(...) { ... }` — when the
   atom-type pair isn't covered by `chematic-ff`'s MMFF94 tables, that
   internal coordinate silently contributes **zero energy and zero gradient**,
   i.e. no restoring force, while VdW repulsion still pushes the atoms apart.
   Reproduced directly: `[C@H](F)(Cl)Br` gives a worst bond of 2.8 Å under
   `etkdg+dreiding` (DREIDING has full generic element coverage) vs **24.3 Å**
   under `etkdg+mmff94` on the identical starting geometry.

The acyclic, non-exotic residual (decane, hexadecane, but-2-ene,
2-chlorobutane, l-serine/l-threonine) blows up under *both* force fields, so
mechanism (1) — or something else not yet isolated in the local/non-rigid
constraint-projection path — is the best-supported explanation for that
subset; this is flagged as an open residual, not fully attributed, and is
listed as a candidate follow-up below rather than fixed here.

**Recommendation:** treat this as **Phase 0**, ahead of the numbered roadmap.
None of Phases 1-8's value can be measured reliably while 67% of a
representative corpus produces torn geometry from the *existing* pipeline —
any Phase 1-2 (real distance geometry) that replaces the current pathway
likely fixes this as a side effect (a real bounds-matrix embedding never
"partially rotates a ring"), but that should be verified, not assumed, once
Phase 1-2 lands. If Phase 1-2 is deferred, Phase 0 should be fixed or the
torsion/constraint-repair step should be disabled by default first.

All numeric metrics below (§3) that depend on geometry validity (RMSD,
energy) are computed **only** over the 15/58 rows that did not hit the
blow-up bucket, and are explicitly caveated as measuring a non-representative
"lucky" subset of the corpus. Chirality (§3) is reported **both** ways —
rolled up over every row that reached that stage (coverage is a structural
fact independent of geometry, so this is valid) **and** restricted to the
geometrically clean subset (the only valid source for a match-rate, since a
signed-volume/dihedral-sign read on a torn molecule is not evidence about
stereo enforcement) — see §3's chirality subsection for why this split
matters before citing either number.

---

## 3. Measured current state — the rest of the requested metrics

Full JSON: `validation/results/etkdg_vs_rdkit_summary.json`. Full per-molecule
rows (nothing dropped): `validation/results/etkdg_vs_rdkit_rows.jsonl`.
Corpus: 58 molecules, `scripts/etkdg_vs_rdkit_gap.py::CORPUS`.

### Chirality retention

**The unconfounded evidence, first.** §1 already showed by source inspection
that chirality/parity is never read anywhere in the generation path
(`grep -n "chiral\|parity\|Chirality"` across `dg.rs`, `etkdg.rs`,
`constraints.rs`, `dg_fft.rs`, `etkdg_knowledge.rs` returns zero matches) —
this alone is proof that declared stereo cannot survive generation by
anything other than accident, independent of any geometry-quality confound.
Two further structural, geometry-independent facts corroborate it:

- **0.0% coverage for implicit-H stereocenters** (30/30 declared,
  `chirality_all_geometry.tetrahedral_implicit_h_lt4_heavy_neighbors`) is a
  hard structural fact, not a sample of a match-rate: `assign_stereo_from_3d`
  (`crates/chematic-3d/src/stereo3d.rs`) requires
  `mol.neighbors(idx).count() == 4`, and `chematic-3d` never materializes
  hydrogens (§1), so any stereocenter with one implicit H — alanine-type
  amino acids, ibuprofen, naproxen, menthol, most drug-like sp3 stereocenters
  — is structurally invisible to chirality verification regardless of
  geometry quality. This is also already flagged in the codebase's own
  docstring (`crates/chematic-py/src/mol_methods.rs`, `stereo_from_coords`:
  *"Only assigns R/S for atoms with four heavy-atom neighbours... Chiral
  centres with an implicit H... are not currently assigned"*).
- **One clean-geometry E/Z inversion**: of the 8 declared E/Z pairs, only
  `chloropropene_E`/`chloropropene_Z` avoided Phase 0's blow-up bucket. On
  that clean pair, `chloropropene_Z` (declared Z) re-perceives as **E** from
  chematic's own valid conformer — a real, geometry-unconfounded mismatch
  (`chirality_clean_geometry.alkene_ez`: match-given-covered 50% on n=2).

**Why the aggregate match-rate numbers below must not be read at face
value.** A signed-volume (tetrahedral) or dihedral-sign (E/Z) read is
mathematically defined on *any* coordinates, including torn ones — so it was
tempting to report a single "match rate given covered" over the whole
corpus. Checking which rows actually populate each bucket shows why that
would overstate the finding: **all 7** of the "4 heavy neighbors" tetrahedral
centers (quaternary_1/2, testosterone, cholesterol) and **6 of the 8** E/Z
pairs fall inside the `chematic_bond_length_blowup` bucket (§2) — a
signed-volume computed on a stretched-to-11-Å bond is not evidence about
stereo *enforcement*, it's evidence about a torn molecule. The table below
therefore reports **both** rollups side by side; only the clean-geometry
column (small n, by construction, until Phase 0 is addressed) is an
unconfounded enforcement signal, and it agrees directionally with the
source-grep proof above.

| stereo class | n declared | chematic coverage (all) | match-given-covered (**all geometry, confounded — see above**) | match-given-covered (**clean geometry only**) | RDKit coverage / match (sanity ceiling) |
|---|---|---|---|---|---|
| tetrahedral, 4 heavy neighbors (assessable) | 7 | 100% | 57.1% | n=0 (all 7 in blow-up bucket) | 100% / 100% |
| tetrahedral, <4 heavy neighbors (implicit H) | 30 | **0.0%** (structural, unconfounded) | n/a | n/a | 100% / 100% |
| alkene E/Z | 8 | 100% | 50.0% | **50.0% (n=2, unconfounded)** | 100% / 100% |

RDKit's 100%/100% column is the sanity-check ceiling confirming the *test
methodology* is sound (RDKit's own embedder reliably reproduces its own
declared stereo on every corpus molecule), not a claim that chematic's
number should also be 100%. Full breakdown:
`chirality_all_geometry` and `chirality_clean_geometry` in the JSON, each
row's `chirality.geometry_clean` flag in the JSONL.

### RMSD, energy, bond/angle geometry (n=15, the non-blown-up subset only)

| metric | value |
|---|---|
| Kabsch-aligned heavy-atom RMSD vs RDKit | mean 1.46 Å, median 1.35 Å, max 3.13 Å |
| "fair" MMFF94 energy delta (see below) | mean **+37.6**, median **+4.4** kcal/mol |
| bond-length violation rate (ext. covalent-radius ref, ±15%) | chematic 0.0%, RDKit 0.0% |
| bond-angle violation rate (ext. hybridization-angle ref, ±12°) | chematic **16.3%**, RDKit 3.1% |

"Fair" energy delta: `chematic`'s own **full** `mmff94_energy_breakdown()`
(bond+angle+stretch-bend+torsion+oop+vdw+electrostatic — the complete
`chematic-ff` implementation, §1) applied to (a) chematic's own conformer and
(b) RDKit's conformer's heavy atoms remapped onto the *same* chematic
`Molecule` topology (verified index correspondence, see Methodology below).
One topology, one energy function, two geometries — an apples-to-apples
geometry-quality signal, *not* a comparison of absolute energies (RDKit's own
native energy is all-atom; chematic's conformer is heavy-atom-only, so those
two absolute numbers are never diffed against each other). Positive means
chematic's own geometry scores worse under chematic's own complete energy
function than RDKit's geometry does — even though chematic's *minimizer*
never actually optimizes the torsion/oop/stretch-bend terms it's being scored
on (§1), which is very likely the dominant cause of the median +4.4 kcal/mol
gap on top of any bond/angle residual.

Ring planarity (`chematic_ring_planarity_rms`/`rdkit_ring_planarity_rms` per
row, RMS deviation of aromatic-ring atoms from their best-fit plane) and
angle-hybridization violations are recorded per molecule in the rows file;
the angle-violation rate above (16.3% vs 3.1%) is consistent with §1's
"no torsion energy term" finding — angles relax reasonably under a spring
potential even without a torsion term, but not as tightly as RDKit's. Ring
planarity specifically is thin this round: most fused-aromatic corpus members
(naphthalene, indole, quinoline, anthracene, pyrene, caffeine, ...) fall in
the §2 blow-up bucket, so the surviving `ok`-bucket aromatic-ring sample is
small (benzene, pyridine, furan) — worth re-measuring once Phase 0 is
addressed rather than trusted as representative today.

### Runtime (single-conformer, noise_sigma=0, deterministic path)

| | chematic | RDKit |
|---|---|---|
| mean | 125.7 ms | 12.5 ms |
| median | 46.4 ms | 3.7 ms |

Wall-clock, so expect run-to-run variance of a similar magnitude to the gap
itself; treat this as "same order of magnitude slower," not a precise ratio.
chematic is roughly an order of magnitude slower per embed on this corpus
(mostly small/medium molecules), on top of producing torn geometry 67% of the
time. Not a primary finding of this RFC (correctness dominates), but relevant context for Phase
1-2 sizing: an O(n³) Floyd-Warshall smoothing step (which the orphaned
`dg_fft.rs` already implements, `DG_MAX_ATOMS = 500`) is not free, and should
be budgeted against this baseline.

### Ensemble generation (n=8 flexible/drug-like molecules, 20 conformers requested each)

Full detail: `ensemble_generation` array in the JSON. Headline: `chematic`'s
retained-conformer count is erratic relative to what was requested —
`ibuprofen` kept **1/20**, `aspirin` kept **4/20** (RDKit: 18/20 and 2/20 on
the same two, respectively) — while `hexadecane` kept 20/20 for both engines.
Diversity (mean pairwise RMSD across the kept ensemble) is directionally
*higher* for chematic than RDKit on every molecule in this subset (e.g.
hexadecane: 4.39 Å vs 2.06 Å), consistent with §1's undirected torsion noise
model rather than a genuine multi-basin conformational search; duplicate rate
within the kept set is ~0 for both engines (the pruning mechanism itself
works — `cluster_conformers_by_rms`/`is_duplicate` do what they say — the
*input* to pruning is the problem).

### Methodology notes (read before citing any number above)

- **Atom correspondence**: both engines parse the *same literal SMILES
  string* with no internal canonicalization on parse, so heavy-atom index i
  in chematic corresponds to heavy-atom index i in RDKit (`Chem.AddHs`
  appends new H atoms after the existing heavy-atom indices, preserving
  their order). This is verified **per molecule, per index**, by an explicit
  element-symbol check (`chematic_syms != rdkit_syms` → `status =
  "atom_correspondence_mismatch"`, fails loud rather than silently trusting
  the mapping); it did not fire on this corpus, but the corpus was
  deliberately hand-written rather than pulled from an external file for
  exactly this reason.
- **RDKit reference**: ETKDGv3 (`AllChem.ETKDGv3()`), `randomSeed=0xF00D`,
  one retry with `useRandomCoords=True` on `-1`; `AllChem.MMFFOptimizeMolecule`
  post-embed. Version `2026.03.3`, via `.venv/bin/python` (symlinked into
  this worktree from the shared repo venv; rebuilt via `maturin develop
  --release` at the start of this round so the installed `chematic` package
  is known to match this exact commit, mitigating the shared-venv race
  against other concurrently-running diagnosis agents).
- **External geometry reference**: RDKit's own periodic-table covalent radii
  × a bond-order scale factor, and RDKit's per-atom hybridization for ideal
  angles — applied identically to both engines' coordinates. This is *not*
  chematic's own internal ideal-bond-length/angle table (`dg.rs`,
  `constraints.rs`, `minimize.rs` each have their own, slightly different,
  copy of such a table — a minor duplication noted but out of scope to fix
  here), so "how close to ideal" is not circular.

---

## 4. Phased roadmap (planning only — no code in this round)

Each phase lists: current state (from §1-3), gap, proposed scope, and an
acceptance signal expressed in terms of *this* diagnosis's own metrics/script
(re-run `scripts/etkdg_vs_rdkit_gap.py`, compare `validation/results/
etkdg_vs_rdkit_summary.json` before/after). None of these are started; this
section is the plan, not a commitment to timeline.

### Phase 0 — fix or bypass the existing blow-up defect (prerequisite, see §2)
**Current**: 67% of the corpus produces a >50%-off-reference bond under the
existing `etkdg.rs` pipeline, in both force fields. **Proposed**: either (a)
increase `satisfy_constraints` iteration budget to what `constraints.rs`'s
own doc says is typical (5-10+) and add ring-membership gating to
`get_torsion_preference`'s call sites, with a regression test asserting
`geometrically_valid_rate == 1.0` on this round's corpus, or (b) if Phase 1-2
is scheduled immediately, verify Phase 1-2's replacement naturally avoids
this class of defect (a proper bounds-matrix embedding has no "rotate one
side of a ring" step to get wrong) before assuming it's moot. **Acceptance**:
`geometrically_valid_rate` (chematic) on the frozen corpus reaches 1.0, or the
replacement path is shown not to reproduce any of the three named mechanisms.

### Phase 1 — bounds matrix
**Current**: no bounds matrix in the live path; a complete, tested
implementation (`build_bound_matrix`, `dg_fft.rs`) exists but has zero
callers. **Proposed**: audit `dg_fft.rs`'s bound construction against
RDKit's ETKDG bounds-generation logic (1-2, 1-3, 1-4 bonded distances; VdW
floor for non-bonded pairs) for completeness/correctness, decide
wire-up-as-is vs. rewrite, and replace `dg::generate_coords` as the base
placement stage of `generate_coords_etkdg`. **Acceptance**: `dg_fft.rs` (or
its replacement) has a caller inside `etkdg.rs`; `#![allow(dead_code)]`
removed.

### Phase 2 — smoothing / triangle-inequality constraints
**Current**: `dg_fft.rs::smooth_bounds` already implements Floyd-Warshall
smoothing (O(n³), gated at `DG_MAX_ATOMS = 500`). **Proposed**: once wired
per Phase 1, benchmark its cost against this round's runtime baseline
(chematic single-conformer mean ~125 ms/embed vs RDKit ~13 ms, §3) on the
corpus and on a larger stress set;
decide whether `DG_MAX_ATOMS` is the right cutoff or whether a sparser
smoothing pass (RDKit uses a bounds-smoothing variant, not naive
all-pairs Floyd-Warshall, for large molecules) is needed. **Acceptance**:
smoothed bounds satisfy the triangle inequality by construction (property
test), runtime on the corpus stays within an agreed multiple of the Phase 0
baseline.

### Phase 3 — stereo constraints
**Current**: does not exist in any form (§1 grep evidence — the strongest,
unconfounded proof); measured 0% coverage for implicit-H stereocenters
(structural, unconfounded) and a 50% match rate on the one clean-geometry E/Z
pair available this round (§3). **Proposed**: add a chiral-volume
sign constraint (or explicit improper-torsion constraint) per declared
tetrahedral/double-bond stereocenter into `build_constraints`/the embedding
objective, keyed off `chematic_core`'s existing CIP/parity data — this
requires deciding how (or whether) to materialize implicit H first, since
today's `Coords3D` is heavy-atom-only and the 4th substituent of most
stereocenters is exactly the implicit H (§1, §3). **Acceptance**: re-running
this round's chirality metric shows >90% match-given-covered on both the
4-heavy-neighbor and (once H are addressed) implicit-H rows, holding
coverage at or above where it started.

### Phase 3 status update (2026-08-11) — supersedes the "Proposed" text above for planning purposes; original left intact for history

**What changed since this RFC was written (2026-07-22):** Phase 0-2's underlying
technical goals were achieved, but through a different code path than this
document originally envisioned. This RFC assumed wiring `dg_fft.rs`'s
`build_bound_matrix`/`smooth_bounds` into the legacy `etkdg.rs`/
`generate_coords_etkdg` pipeline (Phase 1's acceptance criterion literally
says "has a caller inside `etkdg.rs`"). Instead, the Wave 2/3 program's
"Agent C" work (`distance_geometry_v2.rs`, first scaffolded 2026-07-26, four
days after this RFC) built a new module that became the actual production
embedding path (`embed_pipeline_v2` in `pipeline_v2.rs`), and it *does* call
`dg_fft::build_bound_matrix`/`smooth_bounds` directly (confirmed via grep,
2026-08-11) — so those two functions have real production callers now, just
not through `etkdg.rs` as this document assumed. Treat Phase 1/2's
"Current"/"Acceptance" text above as historical, not a live gap: their
technical substance is done, just via `distance_geometry_v2.rs` instead of
`etkdg.rs`.

**Phase 3 itself remains genuinely open**, confirmed independently by two
separate 2026-08 investigations (issue #285, filed for a 2/265-corpus E/Z
sign error; issue #210, re-investigated for a UFF-rescue stereo problem):
`distance_geometry_v2.rs`'s `build_bound_matrix` has zero awareness of
`Atom.chirality` or `BondOrder::Up/Down` — confirmed by direct code
inspection, not just absence-of-evidence. The `stereo_constraints.rs` module
(built in the Wave 2 program, after this RFC) is entirely post-hoc
verify/repair — `verify_stereo`/`repair_stereo` — never a pre-embedding
constraint generator. This RFC's original Phase 3 diagnosis ("does not exist
in any form") is still accurate for the *embedding* side; it's just now also
confirmed true of the production module specifically, not only the legacy
one this RFC was scoped to.

**Committed scope for this phase** (per `ROADMAP.md`'s v0.14.0 plan, S-tier
item 1, added 2026-08-11 — a committed plan, not just a proposal): build
`SMILES stereo -> StereoConstraintSet -> bounds/chiral-volume/dihedral
constraints -> DG -> verify`, covering **both** E/Z and tetrahedral chirality
in one constraint framework, not just the E/Z case #285 happened to surface
first:
- A `StereoConstraintSet` derived from declared stereochemistry
  (`Atom.chirality`, `BondOrder::Up/Down`, and/or CIP data), materialized
  *before* embedding, not after.
- Chiral-volume sign constraints for tetrahedral centers, dihedral/improper-
  torsion sign constraints for E/Z double bonds — injected into
  `build_bound_matrix`'s bounds (or a new constraint stage between bounds
  construction and `smooth_bounds`/refinement), not bolted on as a post-hoc
  rejection filter the way `stereo_constraints.rs` currently is.
- Must handle the implicit-H case this RFC's original Phase 3 flagged (0%
  coverage for implicit-H stereocenters, §3 above) — `Coords3D` is
  heavy-atom-only, so a stereocenter whose 4th substituent is an implicit H
  needs either H materialization or an equivalent geometric proxy.
- Directly reusable by issue #210's UFF-rescue path: 2 of #210's 4
  reproducing molecules (acyclic-bridge substituents) are already fixable
  today with the *existing* `stereo_constraints::repair_stereo` (empirically
  confirmed, no new machinery needed); the other 2 (ring-fused stereocenters,
  e.g. testosterone/cholesterol) have no acyclic bridge for that repair
  strategy and need this new embedding-constraint machinery specifically.
- Motivating evidence beyond #285/#210's own named cases, now measured on
  **both** paths on the identical 58-molecule corpus (29 stereo-bearing,
  re-measured 2026-08-11, filed as issue #291): the legacy `generate_coords`
  first-attempt path is silently stereo-wrong on 13/29 (44.8%); production
  `embed_pipeline_v2` (`UffOnly`+`Ignore`, same "no repair" policy) is
  silently stereo-wrong on **18/29 (62.1%) — higher than legacy, not lower**.
  All 58/58 molecules reported pipeline success either way; `verify_stereo`
  never abstained. This is strong, production-measured (not just directional)
  evidence the missing-constraint mechanism causes real, majority-of-cases
  harm under this policy when nothing catches it — see issue #291 for the
  full breakdown and caveats (single-seed measurement, policy-level not
  variable-controlled comparison between the two structurally different
  embedders).

**Acceptance, more concrete than the original ">90% match-given-covered"**:
re-running this RFC's chirality metric (`scripts/etkdg_vs_rdkit_gap.py` or
the newer `pipeline_v2_vs_rdkit_*` harnesses, whichever is live at
implementation time) shows the constraint-generation stage produces
embeddings that satisfy `verify_stereo` on ≥99% of first-attempt embeds for
molecules with declared stereochemistry, without needing `RepairAndVerify`'s
post-hoc repair pass as a crutch — repair should become the rare exception,
not the primary mechanism, for stereo satisfaction going forward.

**Correction (2026-08-11) — the "injected into `build_bound_matrix`'s
bounds" design above is mathematically impossible, do not attempt it as
written.** A pairwise distance matrix is reflection-invariant: a molecule
and its full mirror image have an *identical* distance matrix (every
pairwise distance is preserved by a reflection). This is not an
implementation gap, it's a property of Euclidean distance geometry itself —
no amount of clever `DistanceBoundAdjustment`-style pairwise bound narrowing
(the mechanism `pipeline_v2.rs` already uses for macrocycle 1-4 relaxation,
`crates/chematic-3d/src/pipeline_v2.rs:629`) can encode which of the two
mirror-image arrangements a bound-matrix-only embedder should produce,
because the bounds themselves cannot distinguish them. Anyone re-deriving a
design from the macrocycle-adjustment precedent will naturally reach for
this approach and should stop here instead.

The design space that actually works operates on **real 3D coordinates**,
not just pairwise distances, since only actual coordinates (not a distance
matrix) can carry a chirality sign at all:
- **Check-and-repair after embedding**: what `stereo_constraints.rs`
  (`verify_stereo`/`repair_stereo`) already does — reflect/rotate a
  bridge-eligible substituent group post-embedding. Structurally cannot fix
  ring-fused stereocenters (no bridge-eligible substituent to move without
  corrupting the ring — issue #210's testosterone/cholesterol case).
- **Check-and-retry with a new stochastic seed**: cheap to wire (the
  existing `max_attempts` retry loop in `embed_distance_geometry_v2_with_
  adjustments` already retries for other failure causes), but success
  probability per attempt is roughly independent-per-stereocenter, so for a
  molecule with *k* declared centers the naive per-attempt success rate is
  roughly on the order of 2^-k — likely still near-zero within a handful of
  attempts for a molecule like cholesterol (8 declared centers). Needs
  measuring before being presented as a fix for the ring-fused case.
- **A genuine chiral-volume/improper-torsion penalty term inside
  `refine_coords`'s SHAKE-like iteration**: since this operates on real
  coordinates each pass (not a static distance matrix), it *can* discriminate
  mirror images the way a real force field's improper-torsion term does.
  Not yet attempted; the harder but mathematically sound path for the
  ring-fused case specifically.

Which of these (or which combination) is the right next step is a policy/
scope decision, not a pure implementation one — see `ROADMAP.md`'s S-tier
item 1 status notes for the live decision state.

### Phase 4 — ring handling
**Current**: `dg.rs` has hand-picked chair/envelope/crown templates by ring
size (6/5/≥8), correctly tested for those cases in isolation, but Phase 0's
mechanism 2 shows the *torsion-preference* layer applied on top is
ring-unaware and actively destroys these templates for polycyclic/fused
cases. **Proposed**: once Phase 1 supplies a real embedding, re-evaluate
whether hand-picked templates are still needed as an initial guess (classical
MDS from a good bound matrix + smoothing typically produces reasonable ring
puckering on its own) or should remain as a fallback/seed. **Acceptance**:
ring-planarity metric (already recorded per molecule) on fused/polycyclic
corpus members stays within RDKit's range; no ring-tearing regressions
(Phase 0's acceptance criterion is a subset of this).

### Phase 5 — experimental torsion preferences
**Current**: `etkdg_knowledge.rs` (1385 lines) already implements a
substantial SMARTS + atom-type torsion-preference library — this is real,
working infrastructure, not a gap — but it is the same code identified in
Phase 0 as ring-unaware and interacting badly with an under-iterated
constraint-repair pass. **Proposed**: after Phase 0/1 land, re-audit
`etkdg_knowledge.rs`'s rule coverage against RDKit's actual ETKDG torsion
library (the CSD-derived preferences RDKit ships) for coverage gaps, now that
applying a rule won't risk tearing the molecule. **Acceptance**: no new
regressions in `geometrically_valid_rate`; measurable reduction in this
round's bond-angle-violation gap (16.3% vs 3.1%) attributable to better
torsion seeding.

### Phase 6 — ensemble generation
**Current**: `ConformerEnsemble`, `generate_conformer_ensemble_with_config`,
Gaussian torsion noise, force-field choice — all exist and work as
documented; the gap is entirely in what's fed into them (Phase 0-5), not the
ensemble machinery itself. **Proposed**: no new machinery; re-measure this
round's ensemble metrics (kept-count erraticism, diversity vs RDKit) after
Phase 0-2 land, since a real embedding + repaired torsion step should also
fix the "keeps 1/20 for ibuprofen" symptom (currently indistinguishable from
"most attempts are torn and discarded by RMSD-pruning coincidence" vs. a
genuine diversity problem). **Acceptance**: kept-conformer count on the
ensemble subset stops varying more than ~2× the requested-vs-kept ratio RDKit
shows on the same molecules.

### Phase 7 — pruning
**Current**: `cluster_conformers_by_rms` (greedy leader-linkage) and
`is_duplicate` (Kabsch-RMSD threshold) both exist, are tested, and do what
they claim — duplicate rate within kept ensembles is already ~0% for both
engines in this round's measurement. **Proposed**: no functional gap
identified; revisit only if Phase 6's re-measurement (after Phase 0-2)
reveals the greedy leader-linkage algorithm (order-dependent by construction)
produces a materially different kept-set than RDKit's pruning on the same
inputs. **Acceptance**: not currently blocking; defer unless Phase 6 surfaces
a concrete problem.

### Phase 8 — force-field refinement
**Current**: `chematic-3d/src/minimize.rs` reimplements a smaller
bond+angle+VdW(+electrostatic) energy function and a fixed-step
finite-difference gradient descent, duplicating (worse) what already exists
as a complete, separately-shipped implementation in `chematic-ff`
(`mmff94_energy::*` with torsion/oop/stretch-bend, `mmff94_minimizer::
minimize_mmff94_lbfgs` with proper L-BFGS + Armijo backtracking, and
`uff::minimize_uff`/`uff_total_energy`). **Proposed**: replace
`chematic-3d`'s `minimize_mmff94`/`minimize_uff` call sites with calls into
`chematic_ff::minimize_mmff94_lbfgs` / `chematic_ff::uff::minimize_uff`
respectively; keep the existing DREIDING path (chematic-ff has no DREIDING
minimizer) but consider whether DREIDING should still be the *default*
(`ConformerForceField::Dreiding`) once a correct, complete MMFF94 minimizer is
one line away. **Acceptance**: this round's "fair" MMFF94 energy delta
(median +4.4 kcal/mol) drops materially, ideally toward parity, since the
minimizer would now optimize the exact terms it's scored on.

---

## 5. Non-goals for this round

- No conformer-generation, minimization, or stereo algorithm code was written,
  even for the three Phase-0 mechanisms pinned down precisely enough to
  reproduce. Fixing them is explicitly deferred to a future, separately
  authorized round.
- No changes to `chematic-ff`'s MMFF94 parameter tables (the missing-coverage
  mechanism in §2.3 is a `chematic-ff` gap, not wired into this RFC's phase
  list, which is scoped to `chematic-3d`'s conformer-generation pipeline per
  the assignment — flagged here for whoever owns `chematic-ff`'s parameter
  completeness).
- No attempt to reconcile this RFC's chirality methodology with the separate,
  already-tracked CIP-correctness workstream (`docs/cip_accurate_rfc.md` and
  the `cip_*` memory entries) — this round only measures "does the declared
  label survive 3D generation," never "is the declared label itself right."
- No production benchmark gate wiring (this is a diagnosis script, not a CI
  check); no scope creep into `chematic-depict`, `chematic-mol`, or any other
  crate.

## 6. Reproduce

```bash
.venv/bin/python -c "import rdkit; print(rdkit.__version__)"   # 2026.03.3
.venv/bin/maturin develop --release -m crates/chematic-py/Cargo.toml
.venv/bin/python scripts/etkdg_vs_rdkit_gap.py
# writes validation/results/etkdg_vs_rdkit_summary.json (aggregate)
#    and validation/results/etkdg_vs_rdkit_rows.jsonl   (every molecule, every status)
```
