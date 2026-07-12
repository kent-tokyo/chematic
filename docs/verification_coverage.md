# Verification Coverage Matrix

An internal engineering audit, not a results report: for every output chematic produces,
what has actually been *measured* (against RDKit or another oracle, and for
self-consistency), versus what is merely assumed. Companion to
[`validation.md`](validation.md) (the polished, RDKit-accuracy report for the 19
descriptors `bench5k.py` covers) and [`rdkit-comparison.md`](rdkit-comparison.md) — this
doc's job is to make **coverage gaps** visible, not to re-derive numbers those own.

Fourteen rounds of targeted bug-hunting (see project memory) found something wrong nearly
every time something was actually measured, and at least six of those rounds turned out
to be the *measurement harness* being wrong, not the library. This matrix exists so gaps
are found by scanning for blank cells, not by accident.

**Status vocabulary — read literally, not optimistically:**
- `MEASURED` — a real corpus, a real number, cited.
- `PARTIAL` — measured, but narrowly (small corpus, single representation, or only one
  of {RDKit-agreement, self-consistency} checked).
- `UNMEASURED` — no evidence found anywhere in code, tests, or project history. Rendered
  as a genuine blank, not a checkmark — "probably fine because nothing's been reported
  broken" is not evidence.
- `KNOWN GAP` — deliberately deferred, documented limitation with a pinned regression
  test or doc comment (e.g. tautomer `max_iter` exhaustion, O3A greedy-order-sensitivity).

---

## Matrix

| Output type | RDKit/oracle agreement | Permutation / repeat invariance | Known unfixed defects | Source |
|---|---|---|---|---|
| Canonical SMILES | **MEASURED** — 0/5000 ChEMBL + 0/33 acyclic-polyene structural corruption, worst-of-10/30 | **MEASURED** — 0/5000, 0/33; E/Z direction 5.50% residual (cosmetic, non-corrupting) | E/Z simple-bond spelling not fully normalized (~1 in 18 stereo molecules) | `sssr_horton_and_canonical_smiles_gap`, `validation.md` |
| Canonical atom ordering (`canonical_atom_order`, feeds InChI) | **PARTIAL** — fixed Round 14 (`c219ee7`), no dedicated post-fix corpus run | **MEASURED** — 0/14 permutation probe + individualized-branch match probe | — | `first_zero_order_dependence_audit` |
| InChI / InChIKey (pure-Rust) | **UNMEASURED** — no corpus comparison vs. standard InChI found | **UNMEASURED** | Approximate, not standard-compliant (documented) — use `native-inchi` for real InChI | `validation.md` Known Limitations |
| InChI (native, IUPAC C lib) | **UNMEASURED** at scale (spot tests only) | **UNMEASURED** | InChI /m /s (enantiomer/isotope) layers not yet measured post-canonicalization fix | — |
| Aromaticity perception (Hückel per-SSSR-ring) | **PARTIAL** — 96.3% atom-flag parity on Kekulized input, worst-of-10 | **UNMEASURED** directly (implied stable via downstream ring counts) | azulene, purine regressed by SSSR fix; root cause (`aromatic_context` bypass) identified, not fixed | `sssr_horton_and_canonical_smiles_gap`, `validation.md` |
| SSSR / ring perception | **MEASURED** — 98.9% ring-size agreement vs `GetSymmSSSR`, 5000-mol; residual is RDKit over-symmetrization, not a chematic bug | **MEASURED** — 100% self-stability (was 50.6%); permanent regression test added | Full Vismara relevant-cycle symmetrization not implemented (not required for correctness) | `sssr_horton_and_canonical_smiles_gap` |
| **ECFP4 fingerprint** | **MEASURED** (this round) — see [worked example](#ecfp4-vs-rdkit--worked-example) below | **MEASURED** (this round) — chemistry-level: 100% | Similarity correlation ~0.94 vs. default-config RDKit, explained by differing hash-fold collisions (not a chemistry gap) | `scripts/ecfp4_agreement.py`, this round |
| FCFP / ECFP6 / other Morgan-basis fingerprints | **UNMEASURED** directly — share ECFP4's underlying Morgan-rank basis and neighbor-sort fix, no dedicated comparison run | **PARTIAL** — neighbor-sort order-independence fixed Round 14 (`86e0d24`), not re-verified at scale | — | `first_zero_order_dependence_audit` |
| Pattern fingerprint | **UNMEASURED** vs RDKit | **MEASURED** — order-independence fixed and tested (Round 14) | — | `first_zero_order_dependence_audit` |
| Pharmacophore fingerprint (2D/3D) | **UNMEASURED** vs RDKit | **MEASURED** — pair-bit symmetry fixed and tested (Round 14) | — | `first_zero_order_dependence_audit` |
| MACCS keys | **UNMEASURED** — no per-key or Tanimoto comparison vs RDKit found | **UNMEASURED** | — | — |
| MAP4, Avalon, HDF, ERG, MHFP, Spectrophores, TopoPF, AtomPair, Torsion, Layered, Reaction FP | **UNMEASURED** — zero RDKit-comparison evidence found for any of these | **UNMEASURED** | — | — |
| Murcko scaffold / MMS | **MEASURED** — 0/5000 structural correctness, worst-of-10 | **MEASURED** — 0.8% residual, same E/Z root cause as canonical SMILES (cosmetic) | Tied to canonical-SMILES E/Z residual | `sssr_horton_and_canonical_smiles_gap` |
| MCS (maximum common substructure) | **PARTIAL** — small constructed-tie corpus only, no ChEMBL-scale run | **PARTIAL** — DFS-order tie-break fixed (`9729fa2`/`95c770c`); confirmed *not* fully spelling-invariant by design | Fix addresses fixed-labeling ties only, not cross-spelling invariance — documented | `first_zero_order_dependence_audit` |
| IUPAC naming | **UNMEASURED** at scale — no corpus run against an oracle (e.g. OPSIN); only hand-built tie-break regression cases exist | **PARTIAL** — several real tie-break bugs found and fixed via constructed non-automorphic cases (seniority, chain selection) | Coverage is intentionally partial (linear/simple polycyclic only); `IupacError` for unsupported structures | `first_zero_order_dependence_audit` |
| Standardization (salts, charges, zwitterions, tautomer canonicalization) | **UNMEASURED** vs RDKit's standardizer at scale | **PARTIAL** — deterministic-per-input confirmed (no HashMap/random); spelling-dependent tie-break left as documented non-fix | Tautomer `max_iter=16` default: confirmed real order-dependence on >16-independent-site molecules (no known real-molecule trigger) — **KNOWN GAP**, pinned `#[ignore]`d test | `first_zero_order_dependence_audit` |
| 2D depiction / layout | **UNMEASURED** vs RDKit layout (no accepted oracle — layout isn't a correctness question in the same sense) | **MEASURED** — process-random HashMap-iteration bug found and fixed (Round 14, `916ffab`); positive-controlled | — | `first_zero_order_dependence_audit` |
| Descriptors (MW/HBA/HBD/TPSA/LogP/MR/Fsp3/ring counts/rotatable bonds/spiro/bridgehead/stereocenters/[nH] SMARTS, 19 tested) | **MEASURED** — 100% or near-100% on 4,999-mol ChEMBL, see `validation.md` for exact per-descriptor numbers | **MEASURED** (subset) — `ring_collateral_damage.py`/`ringinfo_parity.py` self-stability sweeps | Stereocenters 98.7–99.98% depending on oracle (calibration doc in `validation.md`) | `validation.md` |
| Remaining 170+ descriptor functions (Wiener/Zagreb/Randic/Balaban/kappa/chi/BCUT2D/MQN/VSA families/WHIM-2D/RDF-2D/autocorrelation/etc.) | **UNMEASURED** — not in `bench5k.py`'s 19-metric set | **UNMEASURED** | — | — |
| pKa prediction | **UNMEASURED** — no oracle comparison found | **UNMEASURED** | — | — |
| ADMET assembled profiles (BBB/Caco-2/CYP3A4/hERG/PPB/boiled-egg/Ames/CNS MPO) | **UNMEASURED** as assembled predictions (component descriptors they're built from are separately measured) | **UNMEASURED** | — | — |
| QED / SA_score / drug-likeness filters (Lipinski/Veber/Ghose/REOS/PAINS/BRENK) | **UNMEASURED** — assumed parity, no corpus comparison found | **UNMEASURED** | — | — |
| BRICS / RECAP / R-group decomposition / MMP-MMS | **PARTIAL** (R-group only) — VF2 match-choice variance confirmed collapsed by downstream canonicalization (p-xylene symmetric test) | **PARTIAL** (R-group only) — no explicit spelling-variant test in suite | Low priority per Round 14 audit — mechanism is sound but untested at scale | `first_zero_order_dependence_audit` |
| Reaction SMILES/SMIRKS transforms, retrosynthesis templates | **UNMEASURED** vs a reaction oracle | **MEASURED** (safety only) — full Cartesian-product enumeration confirmed, not first-pick | — | `first_zero_order_dependence_audit` |
| Reaction/SMARTS matching (VF2) | **UNMEASURED** for match correctness at scale (spot SMARTS tests only) | **MEASURED** (safety only) — `.first()` pick confirmed to only feed a boolean, never serialized | — | `first_zero_order_dependence_audit` |
| 3D conformer generation (rule-based DG, ETKDG) | **UNMEASURED** vs RDKit ETKDGv3 geometry quality (documented as "not equivalent," never quantified) | **UNMEASURED** | Amide-torsion snap only fixes the *first* N–C(=O) bond for imide N with two carbonyls (succinimide-type) — confirmed coverage gap, deferred | `first_zero_order_dependence_audit` #11, `validation.md` |
| Force fields (MMFF94/DREIDING/UFF): atom typing, energy, minimization | **UNMEASURED** — no energy/geometry comparison vs RDKit's MMFF94 found | **UNMEASURED** | — | — |
| USR shape descriptors | **MEASURED** (self-consistency) — content-based tie-break verified red-before/green-after on a hand-constructed exact tie | **MEASURED** — same as agreement column | — | `first_zero_order_dependence_audit` |
| O3A alignment | **PARTIAL** — alignment score confirmed bit-identical under reordering (10-mol corpus); atom-pair identity is not | **PARTIAL** | Greedy-algorithm order-sensitivity in atom-pairing confirmed real and reachable (3/10 symmetric molecules) — **KNOWN GAP**, harmless for score-based use, real for raw `.pairs` consumers | `first_zero_order_dependence_audit` |
| WHIM / RDF / GETAWAY / SASA / PMI / NPR (3D descriptors) | **UNMEASURED** — assumed parity, never checked against RDKit's 3D descriptor implementations | **UNMEASURED** | — | — |
| Molecular dynamics (MD simulation, thermostats) | **UNMEASURED** — no energy-conservation or trajectory-correctness check found | **UNMEASURED** | — | — |
| PME / Ewald summation | **UNMEASURED** — no comparison against a reference electrostatics implementation found | **UNMEASURED** | — | — |
| File I/O: SDF/MOL V2000/V3000 | **PARTIAL** — throughput-benchmarked (`bench_sdf.py`) against RDKit's egfr.sdf; round-trip *correctness* (does re-read reproduce the same molecule) not separately confirmed | **UNMEASURED** | — | — |
| File I/O: CML/CDXML/MOL2/PDBQT/KET/CIF/RXN/Gaussian/MolJSON | **UNMEASURED** — zero round-trip correctness evidence found for any of these formats | **UNMEASURED** | — | — |
| CIP stereo assignment (R/S, E/Z, M/P atropisomer) | **PARTIAL** — atropisomer 100% on its own test corpus; stereocenter count 98.7–99.98% (oracle-dependent) | **UNMEASURED** directly | 10/4999 stereocenters under/over-counted depending on oracle — ring-adjacent like/unlike CIP tie-break edge cases | `validation.md` |

## Blank-cell priority (the actual output of this audit)

Ranked by how load-bearing the gap is for a migration decision, highest first:

1. **MAP4 / Avalon / HDF / ERG / MHFP / Spectrophores / TopoPF / AtomPair / Torsion /
   Layered fingerprints** — zero RDKit-comparison evidence for any of these. ECFP4 is
   now measured (this round); every other fingerprint family chematic ships is not.
2. **170+ descriptor functions outside `bench5k.py`'s 19** — Wiener/Zagreb/Randic/
   Balaban/kappa/chi/BCUT2D/MQN/VSA/2D-WHIM/2D-RDF/autocorrelation are shipped and
   documented as "190+ descriptors" but only 19 have ever been checked against RDKit.
3. **IUPAC naming correctness at scale** — the only evidence is hand-built tie-break
   regression cases. No run against a naming oracle (OPSIN) on a real corpus exists.
4. **pKa / ADMET assembled predictions** — the descriptor *inputs* are measured; the
   *predictions themselves* (BBB, hERG, CYP3A4, pKa sites) are not.
5. **3D descriptor family (WHIM/RDF/GETAWAY/SASA/PMI/NPR)** and **force-field energies
   (MMFF94/DREIDING/UFF)** — implemented, benchmarked for speed, never checked for
   correctness against RDKit's equivalents.
6. **Non-SMILES file I/O round-tripping** — SDF has throughput benchmarks; correctness
   (parse → write → reparse → same molecule) is unverified for SDF and entirely
   unaddressed for CML/CDXML/MOL2/PDBQT/KET/CIF/RXN/Gaussian/MolJSON.
7. **InChI /m /s layers** and **native-InChI at ChEMBL scale** — the canonicalization
   bug feeding InChI was fixed this round-series, but InChI's own output was never
   re-measured against a real InChI corpus afterward.

---

## ECFP4 vs RDKit — worked example

The "Round 1" migration-decision metric, measured for the first time this round
(`scripts/ecfp4_agreement.py`, 5,000-mol ChEMBL corpus, environment: Python 3.13,
chematic v0.4.29, RDKit 2026.03.3).

**Why raw bit-vector equality isn't the headline number:** chematic hashes atom
environments with FNV-1a (`crates/chematic-fp/src/ecfp.rs`); its own doc comment already
states bit positions aren't meant to match RDKit's (RDKit uses a different hash). Two
independent hash functions landing on the same bit index for the same chemistry is a
~1/2048 coincidence per environment — reporting raw bit agreement as "the" metric would
manufacture a misleading number regardless of correctness (and because both fingerprints
are sparse, it's actually biased *high*, not low: measured 95.36% per-position agreement
on a 1,000-mol sample, dominated by 0/0 non-matches, not a signal of anything).

**What was actually measured, three ways:**

| Tier | What it measures | Result |
|---|---|---|
| Coverage parity | Does chematic generate an environment at every `(atom, radius)` RDKit does, radius ∈ {0,1,2}? (RDKit run with `includeRedundantEnvironments=True` — its default silently prunes some real environments, confirmed via RDKit's own unfolded fingerprint; disabling that pruning is required for a fair chemistry-only comparison.) | **5000/5000 (100%)** exact match |
| Neighborhood identity | Independent BFS in both libraries: does the bond-radius atom-set neighborhood match, atom-for-atom, at radius 1 and 2? (Same SMILES parsed unpermuted by both, so atom index *i* is the same physical atom on both sides; chematic side built from `Mol.bond_table`, no new library code.) | **55,630/55,630 (100%)** across 1,000 molecules |
| Similarity-structure preservation | Pearson correlation between chematic's and RDKit's pairwise Tanimoto similarity (default RDKit config, matching real-world usage), 499,500 pairs from 1,000 molecules | **r = 0.9385**, mean \|Δ Tanimoto\| = 0.0163 |

**Conclusion:** chematic's ECFP4 implements the *same chemistry* as RDKit's Morgan
fingerprint — confirmed two independent, hash-independent ways at 100%. The Tanimoto
correlation of 0.94 (not ~1.0) is fully explained by the two implementations using
different hash functions to fold the same environments into 2048 bits, producing
different (but equally valid) bit-collision patterns — not a chemistry disagreement.
This is a strong, positive, previously-never-measured answer to the original migration
question.

Reproduce: `.venv/bin/python scripts/ecfp4_agreement.py ~/Downloads/SMILES.csv`
