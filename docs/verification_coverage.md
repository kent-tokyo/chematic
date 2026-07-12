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
| **ECFP4 fingerprint** | **MEASURED** (this round) — see [worked example](#ecfp4-vs-rdkit--worked-example) below. Not the standard Rogers-Hahn/RDKit invariant set (includes aromaticity); ~77% invariant-partition match, r=0.94 similarity correlation | **MEASURED** (this round) — self-consistency bug found: 92% of molecules get a different `ecfp4()` for Kekulé vs. aromatic spelling unless `apply_aromaticity()` is called first; ~13% still differ even after calling it | **KNOWN GAP** — 92% naive mismatch is a documented-contract footgun (call `apply_aromaticity()` first); the post-mitigation ~13% residual splits into ~1/3 consistent with the known `aromatic_context` perception bug and ~2/3 **confirmed new defect** (identical aromatic-atom/bond assignment multiset — not just counts — yet still a different fingerprint, verified via 3 independent checks) | `scripts/ecfp4_agreement.py`, this round |
| FCFP4 / ECFP6 (share `initial_atom_id` with ECFP4) | **UNMEASURED** vs RDKit directly, but almost certainly inherit ECFP4's aromaticity-invariant deviation — same seed function, same `atom.aromatic` byte | **MEASURED** (this round, 300-mol spot check) — inherit the SAME self-consistency defect: ECFP6 94.0% and FCFP4 94.0% naive Kekulé-vs-aromatic mismatch (vs. ECFP4's 92.2%) — not ECFP4-specific, a shared-seed-function issue | **KNOWN GAP** — same root cause and same remediation (`apply_aromaticity()`) as ECFP4; post-mitigation residual not separately measured for these two | `scripts/ecfp4_agreement.py`, this round |
| MACCS / Pattern / other Morgan-adjacent fingerprints not sharing `initial_atom_id` | **UNMEASURED** — no dedicated comparison run; not confirmed to share or avoid the aromaticity-representation-dependence defect | **UNMEASURED** | — | `first_zero_order_dependence_audit` (neighbor-sort fix only, Round 14) |
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

## New defect found this round (not a blank cell — a confirmed bug)

Measuring ECFP4 surfaced a real, previously-unknown, unfixed correctness issue, not just
a coverage gap: **`ecfp4()` (and `ecfp6()`/`fcfp4()`, which share the same invariant seed
function) is representation-dependent for ~13% of molecules even after following the
documented `apply_aromaticity()` contract** (see the worked example below). Of that
residual, roughly a third (41/130) is still explained by an aromaticity-perception
disagreement between the two spellings (consistent with the known `aromatic_context`
bug). The other two-thirds (89/130) is **not** — verified three independent ways at
increasing granularity (ring count, then atom/bond counts, then the full
order-independent atom/bond assignment multiset), all converging on the same 89 — a
confirmed, new, unattributed defect with no existing tracking. Not fixed this round
(scope was measurement, not remediation) — flagged here so it doesn't get lost.

---

## ECFP4 vs RDKit — worked example

The "Round 1" migration-decision metric, measured for the first time this round
(`scripts/ecfp4_agreement.py`, 5,000-mol ChEMBL corpus, environment: Python 3.13,
chematic v0.4.29, RDKit 2026.03.3). An earlier pass at this measurement wrongly reported
"100% chemistry-level agreement" based on two checks that didn't actually test what they
claimed to (see [[ecfp4_agreement_methodology]] for the full self-correction — a
connectivity-only BFS check that couldn't detect an invariant difference by construction,
and an unverified assumption about what caused the Tanimoto correlation gap). This
section reflects the corrected, source-verified result.

**Why raw bit-vector equality isn't the headline number:** chematic hashes atom
environments with FNV-1a (`crates/chematic-fp/src/ecfp.rs`); its own doc comment already
states bit positions aren't meant to match RDKit's (RDKit uses a different hash). Two
independent hash functions landing on the same bit index for the same chemistry is a
~1/2048 coincidence per environment — reporting raw bit agreement as "the" metric would
manufacture a misleading number regardless of correctness (and because both fingerprints
are sparse, it's actually biased *high*, not low: measured 95.36% per-position agreement
on a 1,000-mol sample, dominated by 0/0 non-matches, not a signal of anything).

**Finding A — chematic's ECFP4 is a related-but-different fingerprint, not the standard
one.** Its radius-0 invariant explicitly includes `atom.aromatic`
(`crates/chematic-fp/src/ecfp.rs:initial_atom_id`, source-read, not inferred); RDKit's
default Morgan invariant (atomic number, degree, H-count, charge, isotope delta, ring
membership) does not. This is a legitimate design choice — closer in spirit to FCFP than
to standard ECFP — not a bug, but it means chematic's ECFP4 is **not bit-compatible with
RDKit's by design**, not merely by hash-function accident:

| Tier | What it measures | Result |
|---|---|---|
| Coverage parity | Does chematic generate an environment at every `(atom, radius)` RDKit does, radius ∈ {0,1,2}? (RDKit run with `includeRedundantEnvironments=True` to disable its default silent pruning, for a fair comparison.) | **5000/5000 (100%)** exact match — same emission slots |
| Invariant partition agreement | Within each implementation, do environments that hash identically (i.e. "this implementation considers these chemically identical") form the same grouping structure on both sides? Hash-*value*-independent — isolates genuine invariant-encoding disagreement. | **3849/5000 (76.98%)** exact match. Root cause confirmed, not just correlated: aromatic-ring-free molecules match **100%** (363/363); molecules with ≥1 aromatic ring match **75.18%** (3486/4637) — the aromaticity-in-invariant difference fully accounts for the gap, no residual mismatch exists outside it |
| Similarity-structure preservation | Pearson correlation between chematic's and RDKit's pairwise Tanimoto similarity (default RDKit config, matching real-world usage), 499,500 pairs from 1,000 molecules | **r = 0.9385**, mean \|Δ Tanimoto\| = 0.0163 — primarily consistent with (not fully decomposed to) the invariant difference above |
| Connectivity sanity check (auxiliary) | Independently-run BFS in both libraries: does the bond-radius atom-set neighborhood match, atom-for-atom? This checks *parser* agreement only — it never touches fingerprint invariant code and cannot by itself detect an invariant-encoding difference. | **55,630/55,630 (100%)** across 1,000 molecules — parsers agree; not evidence the fingerprints agree |

**Migration-decision answer:** not bit-compatible with RDKit (expected — different hash),
**not the standard ECFP4 definition** (aromaticity-augmented invariant, a deliberate
deviation), similarity correlates strongly but not perfectly (r≈0.94) — **RDKit-trained
similarity thresholds or ML models should not be assumed to transfer without
re-validation.**

**Finding B — a real self-consistency defect, found by testing whether chematic agrees
with *itself* across two spellings of the same molecule** (not an RDKit comparison):
because `atom.aromatic` feeds the invariant (Finding A) and is not auto-perceived for
Kekulé-written SMILES (chematic requires an explicit `apply_aromaticity()` call, unlike
RDKit's auto-sanitize-on-parse), the same molecule can get two different fingerprints
depending on which valid spelling was used to construct it:

| Check | Result |
|---|---|
| Naive (no `apply_aromaticity()`) — aromatic vs. Kekulé spelling of the same molecule, 1,000-mol sample | **922/1000 (92.2%)** get a different `ecfp4()` |
| After calling `apply_aromaticity()` as documented | **130/1000 (13.0%)** *still* mismatch |
| Of that residual: still explained by aromaticity-perception disagreement between the two spellings | **41/130 (~32%)** — consistent with the known `aromatic_context` bug or an extension of it |
| Of that residual: **not** explained by aromaticity perception at any granularity checked | **89/130 (~68%)** — **confirmed new, unattributed defect** |

The 41/89 split was re-derived three times at increasing granularity — ring count, then
aromatic-atom/bond *counts*, then the full order-independent aromatic-atom
`(element, aromatic, degree)` and aromatic-bond `(element, element, bond_type)`
*assignment* multiset (the finest check: counts alone don't rule out two spellings
assigning aromaticity to a *different* set of atoms/bonds while preserving the total,
which the multiset does catch) — all three converged on the identical 89, which is what
earns "confirmed" rather than "not yet explained by the checks tried so far." Full
escalation history: [[ecfp4_agreement_methodology]].

**Blast radius beyond ECFP4**: `initial_atom_id` (the invariant seed carrying the
aromaticity byte) is shared by ECFP6 and FCFP4, and a 300-mol spot check confirms both
inherit the same naive representation-dependence (ECFP6 94.0%, FCFP4 94.0%, vs. ECFP4's
92.2%) — this is a shared-seed-function defect, not ECFP4-specific. Not separately
measured post-`apply_aromaticity()` for ECFP6/FCFP4.

The 92.2% naive case is arguably working as documented — CLAUDE.md/README already state
Kekulé input needs `apply_aromaticity()` first, and this is the first time that contract
has been checked specifically against fingerprint output. The 13.0% residual is not
excused by that contract: a caller following the documented procedure still gets a
different fingerprint for the same molecule. About a third of that residual matches an
already-known, already-deferred limitation; about two-thirds does not and has no existing
tracking — see [New defect found this round](#new-defect-found-this-round-not-a-blank-cell--a-confirmed-bug)
above.

Reproduce: `.venv/bin/python scripts/ecfp4_agreement.py ~/Downloads/SMILES.csv`
