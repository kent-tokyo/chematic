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
| Canonical atom ordering (`canonical_atom_order`, feeds InChI) | **PARTIAL** — fixed Round 14 (`c219ee7`), no dedicated post-fix corpus run | **PARTIAL** — only a 14-case hand-built permutation probe + individualized-branch match probe exists; no full-ChEMBL-scale worst-of-N run of its own. Shares the same individualize-refine code path as `canonical_smiles` (which does have full-corpus worst-of-10/30 coverage — `c219ee7` extracted a shared `winning_individualized_ranks` helper both functions call), so this is *partially*, not *fully*, uncovered — but that inheritance argument itself was never independently verified at scale for this specific function, same "verify analogous fixes independently" gap as the o3a.rs incident | — | `first_zero_order_dependence_audit` |
| InChI / InChIKey (pure-Rust) | **UNMEASURED** — no corpus comparison vs. standard InChI found | **KNOWN GAP** (this round) — 13.4% mismatch on a pure atom-order-only control (seeded, same aromatic origin, no Kekulization), confirmed mostly disjoint (122/134) from the shared Kekulé-origin residual below — a real, independent, previously-unknown order-sensitivity defect, not yet root-caused | Approximate, not standard-compliant (documented) — use `native-inchi` for real InChI; **plus** the order-sensitivity defect above | `validation.md` Known Limitations, `scripts/ecfp4_agreement.py` (tier 6), this round |
| InChI (native, IUPAC C lib) | **UNMEASURED** at scale (spot tests only) | **UNMEASURED** | InChI /m /s (enantiomer/isotope) layers not yet measured post-canonicalization fix | — |
| Aromaticity perception (Hückel per-SSSR-ring) | **PARTIAL** — 96.3% atom-flag parity on Kekulized input, worst-of-10 | **UNMEASURED** directly (implied stable via downstream ring counts) | azulene, purine regressed by SSSR fix; root cause (`aromatic_context` bypass) identified, not fixed | `sssr_horton_and_canonical_smiles_gap`, `validation.md` |
| SSSR / ring perception | **MEASURED** — 98.9% ring-size agreement vs `GetSymmSSSR`, 5000-mol; residual is RDKit over-symmetrization, not a chematic bug | **MEASURED** — 100% self-stability (was 50.6%); permanent regression test added | Full Vismara relevant-cycle symmetrization not implemented (not required for correctness) | `sssr_horton_and_canonical_smiles_gap` |
| **ECFP4 fingerprint** — Layer 1 (definition difference, not a bug): chematic's invariant includes aromaticity, RDKit's default doesn't; structural ceiling, not fixable without an RDKit-compat mode | **MEASURED** (this round) — see [worked example](#ecfp4-vs-rdkit--worked-example) below. ~77% invariant-partition match, r=0.94 similarity correlation — cannot reach ~100% regardless of any future fix, the definitions differ | n/a — single-representation input, ~98.8% unaffected by Layer 2 (1.18% overlap measured, see worked example) | **N/A — design choice**, not a defect. RDKit-compat mode is a feature request, not a fix | `scripts/ecfp4_agreement.py`, this round |
| **ECFP4 fingerprint** — Layer 2 (real bug, independent of Layer 1): representation-dependence | n/a — this is a self-consistency question, not an RDKit comparison | **MEASURED** (this round) — 92% of molecules get a different `ecfp4()` for Kekulé vs. aromatic spelling unless `apply_aromaticity()` is called first (likely the same apply_aromaticity-bypass pattern as Round 8–12's canonical-SMILES/InChI bugs — see worked example); ~13% still differ even after calling it as documented | **KNOWN GAP, corrected scope this round**: the ~13% residual is **not ECFP4-specific** — `canonical_smiles()` (13.2%) and InChI (13.4%) diverge on the *same* molecules (100% overlap with ECFP4's residual set both ways, n=1000), so this is one shared, systemic defect touching ≥3 core output functions, not three coincidentally-similar ones. Mechanism is **unidentified**: not aromatic-atom/bond flag assignment for 89/130 (identical multiset, still differs), not SSSR ring decomposition for any of the 130 (identical ring-size multiset), and — for ECFP4 specifically — not plain atom-order-dependence (0.0% mismatch on an order-only control, vs. 13.0% residual). `canonical_smiles` (3.8%) and, more substantially, InChI (13.4%) each carry their *own separate* order-sensitivity defect — confirmed via molecule-set intersection, not assumed from magnitude alone, to be mostly disjoint from this shared residual (37/38 and 122/134 outside, respectively). InChI's is large enough to be its own new, unattributed finding — see InChI row above and worked example | `scripts/ecfp4_agreement.py` (tier 6), this round |
| FCFP4 / ECFP6 (share `initial_atom_id` with ECFP4) | **UNMEASURED** vs RDKit directly, but almost certainly inherit ECFP4's aromaticity-invariant deviation — same seed function, same `atom.aromatic` byte | **PARTIAL** — this round's 300-mol spot check confirms they inherit ECFP4's Layer-2 representation-dependence defect (ECFP6 94.0%, FCFP4 94.0% naive Kekulé-vs-aromatic mismatch). But that's a *different axis* from Round 14's neighbor-sort order-independence fix (`86e0d24`, originally applied to ECFP4/pattern-fp) — **no dedicated Rust test confirms FCFP4/ECFP6 inherited that fix too**; the only FCFP/ECFP6-specific test found (`ecfp6_vs_ecfp4_benzene_differ`) checks bit-count difference, not order-independence. "MEASURED — order-independence fixed and tested" was too strong for this axis; downgraded to PARTIAL | **KNOWN GAP** (Layer-2 representation-dependence) — same root cause and same remediation (`apply_aromaticity()`) as ECFP4, post-mitigation residual not separately measured; **UNMEASURED** (Round-14 neighbor-sort inheritance) — plausible but unverified | `scripts/ecfp4_agreement.py`, this round; `first_zero_order_dependence_audit` |
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
| Descriptors (MW/HBA/HBD/TPSA/LogP/MR/Fsp3/ring counts/rotatable bonds/spiro/bridgehead/stereocenters/[nH] SMARTS, 19 tested) | **MEASURED** — 100% or near-100% on 4,999-mol ChEMBL, see `validation.md` for exact per-descriptor numbers. **Single-representation only** — `bench5k.py` parses each corpus SMILES once (no `doRandom`/worst-of-N in the RDKit-agreement comparison itself; confirmed by reading the file) | **PARTIAL** (9/19, not all) — `ring_collateral_damage.py`/`ringinfo_parity.py` worst-of-10 self-stability sweeps cover mol_wt/tpsa/hba/hbd/ring_count/num_aromatic_rings/num_saturated_rings/num_aliphatic_rings/logp/mr/scaffold; **zero self-stability coverage despite being perception-sensitive**: num spiro atoms, num bridgehead atoms, num amide bonds, aromatic/aliphatic heterocycle counts, stereocenters | Stereocenters 98.7–99.98% depending on oracle (calibration doc in `validation.md`); the 6 zero-coverage descriptors above are a real, unaudited gap, not a known-and-accepted one | `validation.md`, `scripts/bench5k.py`, `scripts/ring_collateral_damage.py` |
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

## MEASURED-column audit (2026-07-12 continuation)

Round 15's ECFP4 measurement went through four rounds of self-correction (a tautological
check, then two successively-too-coarse residual classifiers) before its "MEASURED"
claim was actually earned. That raised an obvious question: do *other* rows' MEASURED
labels — especially ones set in early rounds, before this project's positive-control
discipline existed — hold up to the same scrutiny? Three parallel research passes
audited the matrix's MEASURED/PARTIAL rows against the actual scripts, tests, and git
history behind each claim.

**A meta-finding worth stating plainly**: several proposed corrections from that research
were themselves wrong, caught only by directly re-checking README.md's own prose and the
repo (not by trusting the research's summary). Canonical SMILES, Murcko scaffold, and
Aromaticity perception were each flagged as under-evidenced; all three turned out to have
solid, already-documented positive controls or dedicated scripts the research had missed
(a pre-fix-vs-post-fix historical comparison for canonical SMILES; a caught-and-fixed
measurement-harness bug for Murcko scaffold; a dedicated `aromaticity_atom_parity.py` for
aromaticity) — this is the exact discipline the audit itself exists to apply, applied one
level up, to the audit's own output.

**Confirmed solid, no change**: Canonical SMILES, Murcko scaffold, Aromaticity
perception, SSSR/ring perception, Pattern FP, Pharmacophore FP, 2D depiction/layout, USR,
Standardization (common case), IUPAC naming, MCS, InChI (both variants) — existing labels
already reflect the real evidence (genuine non-automorphic test cases, git-verified
red-before-green, or accurately UNMEASURED/thin where the evidence really is thin).

**Genuine gaps found and corrected in the matrix rows above**: the Descriptors row
implied uniform worst-of-N coverage across all 19 `bench5k.py` metrics — the RDKit
comparison itself is single-representation only, and the self-stability companion
scripts cover 9 of 19, leaving 6 perception-sensitive descriptors (spiro/bridgehead
atoms, amide bonds, aromatic/aliphatic heterocycles, stereocenters) with zero
self-stability evidence despite being exactly the kind of descriptor most likely to be
representation-sensitive. Canonical atom ordering and FCFP4/ECFP6 were both downgraded
from overstated MEASURED/PARTIAL claims to precisely-scoped PARTIAL, each for a distinct
reason (see their rows above).

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
8. **Self-stability for 6 perception-sensitive descriptors** (found this audit pass) —
   num spiro atoms, num bridgehead atoms, num amide bonds, aromatic/aliphatic heterocycle
   counts, and stereocenters have RDKit-agreement numbers but zero worst-of-N
   self-stability evidence, unlike the other 13 of `bench5k.py`'s 19 descriptors.

## New defect found this round (not a blank cell — a confirmed bug)

Measuring ECFP4 surfaced a real, previously-unknown, unfixed correctness issue — this is
**Layer 2 only** (see the worked example below for the full Layer 1/Layer 2 split; Layer
1, chematic's ECFP4 including aromaticity in its invariant, is a design choice, not a
bug, and is unrelated to what follows): **`ecfp4()` (and `ecfp6()`/`fcfp4()`, which share
the same invariant seed function) is representation-dependent for ~13% of molecules even
after following the documented `apply_aromaticity()` contract**.

**Corrected this round: this is not an ECFP4-specific defect.** A follow-up
cross-consumer check (`scripts/ecfp4_agreement.py` tier 6) found `canonical_smiles()`
and InChI diverge on the *same* molecules ECFP4 does, post-`apply_aromaticity()`, at
near-identical rates (n=1000: ECFP4 130, canonical_smiles 132, InChI 134 residual
mismatches; 100% of ECFP4's residual set is contained in both of the others'). The
earlier framing of "89/130 confirmed new, unattributed ECFP4 defect" undersold what this
actually is — a single shared root cause affecting at least three independent core
output functions identically, not a fingerprint-specific quirk.

Of the residual, roughly a third (41/130) still shows an aromatic-atom/bond
flag-assignment disagreement between the two spellings — the same *class* of bug as the
known `aromatic_context` regression (azulene, purine), though whether it's literally the
same code path is **unverified**: `aromatic_context` was previously characterized via
only 2 molecules, and a rate this much higher (≈4% of the whole corpus, not 2 isolated
cases) suggests this may be a broader flag-perception issue than what was previously
scoped, not confirmed either way this round. The other two-thirds (89/130) has an
**identical** aromatic-atom/bond assignment multiset yet still diverges — verified three
independent ways at increasing granularity (ring count, then atom/bond counts, then the
full order-independent assignment multiset), all converging on the same 89.

This round additionally ruled out two more candidate mechanisms for the full 130-molecule
residual: **SSSR ring decomposition** (0/130 differ in ring-size multiset between the two
spellings — the rings found are identical, so this isn't a ring-finding bug) and, for
ECFP4 specifically, **plain atom-order-dependence** (two aromatic-preserving, non-Kekulized
respellings of the same molecule — differing only in traversal order, not origin — mismatch
in ECFP4 0.0% of the time, vs. the 13.0% residual; this rules out "it's just an order
artifact of the Kekulé round-trip" as ECFP4's explanation and supports that the divergence
tracks Kekulé-vs-aromatic *origin* specifically).

**A second, previously-unknown defect surfaced while checking this, and is reported
separately rather than folded into the number above**: `canonical_smiles` and InChI each
carry their own independent order-sensitivity, confirmed (not assumed) disjoint from the
shared 13% residual by intersecting the actual molecule sets rather than just comparing
percentages. `canonical_smiles`'s is small (3.8%, 37/38 outside its own residual set).
InChI's is **large — 13.4%, coincidentally the same magnitude as InChI's residual count,
but confirmed to be a mostly different set of molecules (only 12/134 overlap)**. A first
draft of this section nearly reported the `canonical_smiles` case as "self-evidently
separate" purely from the 3.8 << 13.2 arithmetic, without checking set membership — that
would have missed the far more consequential InChI case, where the matching magnitude
(134 = 134) would have supported the wrong conclusion ("InChI's residual is just
order-dependence") if the sets had gone unchecked. **InChI's 13.4% order-sensitivity is a
new, real, unattributed defect in its own right** — logged as a new blank cell in the
InChI row above, not yet root-caused; a plausible but unverified hypothesis is that it
inherits and amplifies `canonical_smiles`'s ordering sensitivity through InChI's
additional numbering/stereo-layer steps on chematic's shared canonical-atom-order code
path.

**Net: the shared residual's mechanism remains unidentified.** Per-atom flag assignment,
per-bond *type* assignment (the multiset keys on bond type, not just element pair), and
SSSR ring topology are all ruled out as order-independent-multiset-level explanations for
the 89/130 majority; order-dependence is ruled out for ECFP4's share of the whole 130.
Two candidates remain, **not distinguished this round, and not ranked** — a one-molecule
spot check attempted to localize the divergence (via `morgan_fp_counts` hash-set diff,
folded to a bit index, resolved to specific atom indices) but atom-table indices don't
correspond 1:1 between the two spellings (the Kekulé-origin respelling has a different
internal atom order), so per-index atom comparison couldn't distinguish the candidates
without a real atom-correspondence step (e.g. a substructure match), which wasn't
attempted: (a) a *positional* bond-structure difference below the multiset's resolution
(same type counts, different placement), or (b) a symmetric aromaticity-flag swap between
two atoms of identical `(atomic_number, degree)` — the multiset's own documented blind
spot, which would point back at aromaticity perception rather than bond order. Do not
treat either as confirmed. Not fixed this round (scope was measurement/characterization,
not remediation) —
flagged here so it doesn't get lost, with
the corrected (broader, shared) scope so it isn't mistaken for an ECFP4-only issue in a
future round.

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

**Correcting the premise, not just the number.** Round 1's original instruction — "measure
ECFP4-vs-RDKit bit-agreement; that number is the migration-decision signal" — assumed
chematic's ECFP4 shared RDKit's exact invariant definition. That assumption was never
checked before being issued, and it was wrong (see Layer 1 below). The two things this
round actually found need to be reported as **two independent layers, not one blended
number** — conflating them (as the first draft of this section did) makes a design choice
look like a bug backlog, and hides a real bug behind a design choice:

- **Layer 1 (definition difference — not a bug)**: chematic's ECFP4 deliberately includes
  aromaticity in its invariant; RDKit's default doesn't. Different definitions, so
  RDKit-parity was never structurally reachable — no future fix moves Tier 1/2/3's numbers
  toward 100%, only a new RDKit-compatible mode would (a feature request, not a bug fix).
- **Layer 2 (real bug — representation-dependence)**: `ecfp4()` gives a different value for
  the Kekulé vs. aromatic spelling of the identical molecule. Independent of Layer 1, and a
  genuine fingerprint-contract violation regardless of which invariant definition is used.

**These two layers are ~98.8% cleanly separated by construction, confirmed empirically, not
assumed**: Tier 1/2/3 below use one fixed representation per molecule (the corpus's raw,
as-given SMILES) and never vary it — Layer 2's mechanism (Kekulé input skipping
`apply_aromaticity()`) requires *two* representations to manifest, so it should barely touch
Tier 1/2/3's numbers if the corpus is already aromaticity-perceived on parse. Checked
directly: forcing `apply_aromaticity()` on the raw corpus input changes `ecfp4()` output for
only **59/5000 molecules (1.18%)**. Tier 1/2/3's 76.98%/r=0.94 numbers are therefore
**~98.8% purely Layer 1** — the small remainder is a negligible overlap, not evidence the
two layers were conflated in the measurement itself, only in how the first draft described it.

**Why raw bit-vector equality isn't the headline number:** chematic hashes atom
environments with FNV-1a (`crates/chematic-fp/src/ecfp.rs`); its own doc comment already
states bit positions aren't meant to match RDKit's (RDKit uses a different hash). Two
independent hash functions landing on the same bit index for the same chemistry is a
~1/2048 coincidence per environment — reporting raw bit agreement as "the" metric would
manufacture a misleading number regardless of correctness (and because both fingerprints
are sparse, it's actually biased *high*, not low: measured 95.36% per-position agreement
on a 1,000-mol sample, dominated by 0/0 non-matches, not a signal of anything).

**Layer 1 — chematic's ECFP4 is a related-but-different fingerprint, not the standard
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

**Layer 2 — a real self-consistency defect, found by testing whether chematic agrees
with *itself* across two spellings of the same molecule** (not an RDKit comparison, and
independent of Layer 1's invariant-definition question): because `atom.aromatic` feeds the
invariant and is not auto-perceived for Kekulé-written SMILES (chematic requires an
explicit `apply_aromaticity()` call, unlike RDKit's auto-sanitize-on-parse), the same
molecule can get two different fingerprints depending on which valid spelling was used to
construct it:

**Likely the same root cause as Round 8–12's canonical-SMILES/InChI bugs, for the naive
92% case specifically.** Those historical bugs were code paths that didn't route through
`apply_aromaticity()` normalization before comparison/canonicalization. The naive-mismatch
number here (92.2%, before any mitigation) drops to 13.0% once `apply_aromaticity()` is
called first — a 79-point reduction directly consistent with "most of this is the exact
same normalization-bypass pattern, recurring in a new code path (`ecfp4()`) instead of a
new mechanism." Not yet verified beyond this consistency check (no direct code-path
comparison against the Round 8–12 fix was done this round) — the 13.0% residual is
explicitly **not** covered by this hypothesis, since `apply_aromaticity()` was already
called for those cases and they still mismatch:

| Check | Result |
|---|---|
| Naive (no `apply_aromaticity()`) — aromatic vs. Kekulé spelling of the same molecule, 1,000-mol sample | **922/1000 (92.2%)** get a different `ecfp4()` |
| After calling `apply_aromaticity()` as documented | **130/1000 (13.0%)** *still* mismatch |
| Of that residual: still explained by aromaticity-perception disagreement between the two spellings | **41/130 (~32%)** — same *class* of bug as the known `aromatic_context` regression; literal same-code-path attribution unverified (see below) |
| Of that residual: **not** explained by aromaticity perception at any granularity checked | **89/130 (~68%)** — confirmed via 3 independent checks, mechanism still unidentified (see tier 6) |

The 41/89 split was re-derived three times at increasing granularity — ring count, then
aromatic-atom/bond *counts*, then the full order-independent aromatic-atom
`(element, aromatic, degree)` and aromatic-bond `(element, element, bond_type)`
*assignment* multiset (the finest check: counts alone don't rule out two spellings
assigning aromaticity to a *different* set of atoms/bonds while preserving the total,
which the multiset does catch) — all three converged on the identical 89, which is what
earns "confirmed" rather than "not yet explained by the checks tried so far." Full
escalation history: [[ecfp4_agreement_methodology]].

**Tier 6 follow-up — is this ECFP4-specific, and is it really about aromaticity?** Two
questions the 41/89 split above doesn't answer on its own: does this residual hit other
`apply_aromaticity()`-consuming functions on the same molecules (shared root cause vs.
ECFP4-specific), and is "Kekulé-vs-aromatic origin" really the operative variable, or is
it confounded with plain atom-order-dependence (RDKit's non-canonical Kekulé respelling
changes *both* the origin and the atom traversal order at once — a confound the naive
check above never controlled for)?

| Check | Result |
|---|---|
| `canonical_smiles()` residual (post-`apply_aromaticity()`), same 1,000-mol sample | **132/1000 (13.2%)** |
| InChI residual (post-`apply_aromaticity()`), same sample | **134/1000 (13.4%)** |
| Overlap: ECFP4's residual set ∩ `canonical_smiles`'s residual set | **130/130 (100%)** of ECFP4's residual molecules also mismatch in `canonical_smiles` |
| Overlap: ECFP4's residual set ∩ InChI's residual set | **130/130 (100%)** of ECFP4's residual molecules also mismatch in InChI |
| SSSR ring-size multiset, restricted to ECFP4's 130-molecule residual set | **0/130 differ** — ring decomposition is identical between the two spellings for every residual molecule |
| Order-only control (two aromatic-preserving, non-Kekulized respellings of the same molecule, seeded — origin held fixed, only atom traversal order varies), ECFP4 | **0/1000 (0.0%)** mismatch |
| Same order-only control, `canonical_smiles` | **38/1000 (3.8%)** mismatch — of which only **1/38** falls inside `canonical_smiles`'s own residual set, 37/38 outside it |
| Same order-only control, InChI | **134/1000 (13.4%)** mismatch — numerically equal to InChI's residual count, but only **12/134** falls inside InChI's own residual set, 122/134 outside it (a coincidence of magnitude, not the same molecule set — checked, not assumed) |

Four findings from this table, each independently load-bearing:

1. **This is a shared, systemic defect, not an ECFP4 quirk.** 100% of ECFP4's residual
   set reappears in both `canonical_smiles`'s and InChI's own residual sets, at
   near-identical magnitude (13.0% / 13.2% / 13.4%). One root cause, three consumers —
   the earlier "confirmed new ECFP4 defect" framing was too narrow.
2. **SSSR ring decomposition is ruled out** for the entire 130-molecule residual, not
   just the 89-subset that already had matching flags: every residual molecule has an
   identical ring-size multiset between the two spellings. Whatever differs, it isn't
   which rings get found.
3. **Atom-order-dependence is ruled out as ECFP4's explanation.** If the 13.0% residual
   were mostly a side effect of RDKit's non-canonical Kekulé respelling also reshuffling
   atom order (rather than the Kekulé-vs-aromatic origin itself), the order-only control
   — which reshuffles order without ever touching Kekulé/aromatic origin — should show a
   comparable mismatch rate. It doesn't (0.0%), so the residual really does track origin,
   not incidental ordering, for ECFP4.
4. **`canonical_smiles` and InChI each carry their own separate order-sensitivity defect,
   confirmed disjoint from the shared residual, not assumed.** An earlier draft of this
   section asserted `canonical_smiles`'s 3.8% order-only mismatch was "layered on top, not
   competing for the same explanatory share" purely from the numbers (3.8 << 13.2) —
   without checking whether those 38 molecules were even inside the 132-molecule residual
   set. Checked directly: only 1/38 is. The independence claim holds for
   `canonical_smiles`, but it holds *because it was verified*, not because the arithmetic
   made it self-evident. The same check on InChI caught a real near-miss: InChI's
   order-only mismatch count (134) is numerically identical to its residual count (134),
   which would have supported a much stronger (and wrong) claim — "InChI's residual is
   pure order-dependence" — if left unchecked. It isn't: only 12/134 overlap. InChI has a
   **large, independent order-sensitivity defect (13.4%, larger than `canonical_smiles`'s
   3.8%)**, disjoint from the shared Kekulé-origin residual — a fourth, previously-unknown
   finding, not yet root-caused. A plausible but unverified hypothesis is that InChI
   inherits and amplifies `canonical_smiles`'s ordering sensitivity through additional
   numbering/stereo-layer steps built on top of chematic's shared canonical-atom-order code
   path — consistent with, not confirmed by, [[feedback_canonicalization_theory]]'s
   tie-break concerns.

**Mechanism of the shared residual: still unidentified.** Per-atom aromatic-flag
assignment, per-bond *type* assignment, and SSSR ring topology are all ruled out for
89/130 (the multiset keys on bond type, not just element pair — it rules out more than a
loose "counts matched" reading would suggest); atom-order-dependence is ruled out for
ECFP4's share of all 130. A one-molecule spot check (`morgan_fp_counts` hash-set diff,
folded to a bit index, resolved to atom indices via `ecfp_bitinfo`) attempted to localize
the divergence further but stalled on atom-index correspondence — the Kekulé-origin
respelling has a different internal atom order, so raw per-index atom-table comparison
isn't meaningful without a real correspondence step (e.g. a substructure match), which
wasn't attempted this round. **Two candidates remain, deliberately left unranked**: (a) a
*positional* bond-structure difference below the multiset's resolution (same type counts,
different placement — would require reading `crates/chematic-perception`'s
`apply_aromaticity()` and `crates/chematic-smiles/src/canonical.rs`'s bond-order handling
directly), or (b) a symmetric aromaticity-flag swap between two atoms of identical
`(atomic_number, degree)` — the multiset's own documented blind spot, which would point
back at aromaticity perception, not bond order. An earlier draft of this section named (a)
as "the leading candidate" with no evidence distinguishing it from (b) — corrected here;
neither is confirmed. InChI's own 13.4% order-sensitivity defect (finding 4) is a
**separate, untouched blank cell**, not part
of this mechanism question.

**Blast radius beyond ECFP4**: `initial_atom_id` (the invariant seed carrying the
aromaticity byte) is shared by ECFP6 and FCFP4, and a 300-mol spot check confirms both
inherit the same naive representation-dependence (ECFP6 94.0%, FCFP4 94.0%, vs. ECFP4's
92.2%) — this is a shared-seed-function defect, not ECFP4-specific. Not separately
measured post-`apply_aromaticity()` for ECFP6/FCFP4.

The 92.2% naive case is arguably working as documented — CLAUDE.md/README already state
Kekulé input needs `apply_aromaticity()` first, and this is the first time that contract
has been checked specifically against fingerprint output. The 13.0% residual is not
excused by that contract: a caller following the documented procedure still gets a
different fingerprint for the same molecule — and, per the tier 6 findings above, the same
molecule also canonicalizes and InChI-generates differently. About a third of that residual
matches the same *class* of bug as an already-known, already-deferred limitation; about
two-thirds does not and has no existing tracking, mechanism unidentified — see
[New defect found this round](#new-defect-found-this-round-not-a-blank-cell--a-confirmed-bug)
above.

Reproduce: `.venv/bin/python scripts/ecfp4_agreement.py ~/Downloads/SMILES.csv --aromaticity-sample 1000 --layer2-sample 1000`
