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
| InChI / InChIKey (pure-Rust) | **UNMEASURED** for the InChI *string* vs. standard InChI (numbering/`/h`-layer/`/m`-layer remain structurally non-standard, documented, unchanged this round). **MEASURED** for the underlying R/S it encodes, newly this round: `assign_cip()` vs RDKit's per-atom CIP oracle, 4163 stereocenters, 5000-mol corpus — 96.83% (was 76.22% before this round's fix; see "Defect found and FIXED" section below) | **MEASURED, FIXED this round** — order-only mismatch 13.4%→3.5% (n=1000); InChI-*specific* excess over `canonical_smiles`'s own baseline order-sensitivity: 74%→**0%** (residual 3.5% fully explained by shared ordering-layer churn, no separate InChI-specific mechanism detected at this sample size) | Approximate, not standard-compliant (documented) — use `native-inchi` for real InChI; **plus** ~132-155/4163-4186 stereocenters (3.2-3.7%, oracle-dependent) still wrong vs RDKit — attempted and reverted (not just deferred): the residual is at least 4 distinct mechanisms, and a follow-up fix for one (triple-bond duplication) went net negative, revealing the comparator itself (`cip_branch_spheres`, shell-multiset pooling, not true recursive CIP) is the real limiter, not a missing per-bond-type rule — design RFC for the proper fix at `docs/cip_accurate_rfc.md`, frozen 155-case corpus at `validation/cip_label_corpus.jsonl`, see "Remaining, explicitly out of scope" below | `validation.md` Known Limitations, `scripts/ecfp4_agreement.py` (tier 6), "Defect found and FIXED this round: `assign_cip()`" section below, `docs/cip_accurate_rfc.md` |
| InChI (native, IUPAC C lib) | **UNMEASURED** at scale (spot tests only) for the InChI string itself. Its stereo *input* is fixed this round: `crates/chematic-inchi/src/native/convert.rs` calls `tetrahedral_stereo_neighbors()` directly, the same function fixed above, so native-inchi's R/S input inherits the same 76.22%→96.83% correctness improvement — not independently re-measured against the C library's own output, only inferred from the shared code path | **UNMEASURED** | InChI /m /s (enantiomer/isotope) layers not yet measured post-canonicalization fix | — |
| Aromaticity perception (Hückel per-SSSR-ring) | **PARTIAL** — 96.3% atom-flag parity on Kekulized input, worst-of-10 | **UNMEASURED** directly (implied stable via downstream ring counts) | azulene, purine regressed by SSSR fix; root cause (`aromatic_context` bypass) identified, not fixed | `sssr_horton_and_canonical_smiles_gap`, `validation.md` |
| SSSR / ring perception | **MEASURED** — 98.9% ring-size agreement vs `GetSymmSSSR`, 5000-mol; residual is RDKit over-symmetrization, not a chematic bug | **MEASURED** — 100% self-stability (was 50.6%); permanent regression test added | Full Vismara relevant-cycle symmetrization not implemented (not required for correctness) | `sssr_horton_and_canonical_smiles_gap` |
| **ECFP4 fingerprint** — Layer 1 (definition difference, not a bug): chematic's invariant includes aromaticity, RDKit's default doesn't; structural ceiling, not fixable without an RDKit-compat mode | **MEASURED** (this round) — see [worked example](#ecfp4-vs-rdkit--worked-example) below. ~77% invariant-partition match, r=0.94 similarity correlation — cannot reach ~100% regardless of any future fix, the definitions differ | n/a — single-representation input, ~98.8% unaffected by Layer 2 (1.18% overlap measured, see worked example) | **N/A — design choice**, not a defect. RDKit-compat mode is a feature request, not a fix | `scripts/ecfp4_agreement.py`, this round |
| **ECFP4 fingerprint** — Layer 2 (real bug, independent of Layer 1): representation-dependence | n/a — this is a self-consistency question, not an RDKit comparison | **MEASURED, 89/130 FIXED this round** — naive Kekulé-vs-aromatic mismatch (92%) drops to **41/1000 (4.1%)** post-`apply_aromaticity()`, down from 130/1000 (13.0%). Root cause: `implicit_hcount()`'s aromatic-path heuristic can't distinguish pyrrole-type from pyridine-type heteroatoms once `apply_aromaticity_ex()` normalizes bond order — fixed by freezing the pre-normalization H count for atoms where it would otherwise change (`crates/chematic-perception/src/aromaticity.rs`). Verified causally: `canonical_smiles`/InChI residuals dropped the *same* 89 molecules simultaneously (130→41, 132→43, 134→45); MW verified against RDKit at full-corpus scale (4998/5000, the 2 exceptions pre-existing/unrelated), not just internal agreement. Remaining 41/1000 is the pre-existing, distinct `aromatic_context`-class flag-assignment defect (untouched, out of scope) — see [defect writeup](#defect-found-and-fixed-this-round-implicit-h-count-lost-across-the-kekuléaromatic-boundary) | **KNOWN GAP (41/1000 remaining)** — genuine aromatic-atom/bond flag-assignment disagreement, same class as `aromatic_context` (azulene, purine), same-code-path attribution still unverified. `canonical_smiles` (3.8%) and InChI (13.4%) each separately carry their own order-sensitivity defect, confirmed disjoint from this mechanism by molecule-set intersection — untouched by this fix, still open | `scripts/ecfp4_agreement.py` (tier 6), `scripts/aromaticity_mechanism_probe.py` (new), this round |
| FCFP4 / ECFP6 (share `initial_atom_id` with ECFP4) | **UNMEASURED** vs RDKit directly, but almost certainly inherit ECFP4's aromaticity-invariant deviation — same seed function, same `atom.aromatic` byte | **PARTIAL, 89/130 FIXED this round (verified independently, not assumed)** — this round's spot check confirms both inherit ECFP4's Layer-2 representation-dependence defect (naive Kekulé-vs-aromatic mismatch) AND independently confirms both inherited the `implicit_h` fix: post-`apply_aromaticity()` residual dropped from ECFP6 94.0%/FCFP4 94.0% (naive) to **41/1000 (4.1%) for both**, matching ECFP4's own post-fix number exactly (checked directly via `ecfp6()`/`fcfp4()` calls, not inferred from ECFP4's result). Round-14 neighbor-sort order-independence inheritance (`86e0d24`) is a *different axis*, still **no dedicated Rust test** confirming it for FCFP4/ECFP6 specifically — remains PARTIAL for that reason | **KNOWN GAP** (Layer-2 representation-dependence) — same root cause and same remediation as ECFP4, remaining 41/1000 shared with the still-open `aromatic_context`-class defect; **UNMEASURED** (Round-14 neighbor-sort inheritance) — plausible but unverified | `scripts/ecfp4_agreement.py`, `scripts/aromaticity_mechanism_probe.py`, this round; `first_zero_order_dependence_audit` |
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
| CIP stereo assignment (R/S, E/Z, M/P atropisomer) | **PARTIAL** — atropisomer 100% on its own test corpus; stereocenter count 98.7–99.98% (oracle-dependent) | **UNMEASURED** directly | 10/4999 stereocenters under/over-counted depending on oracle — ring-adjacent like/unlike CIP tie-break edge cases; **`assign_ez`/`substituent_is_up` never read the `bond_direction` side channel** — an exocyclic double bond (e.g. imine `/N=c...`) whose geometry is anchored on an aromatic ring bond is invisible to CIP on every parse, even though an independent RDKit check confirms it's real, stable stereochemistry (`STEREOZ`), not a chematic-specific convention. Found while fixing Canonical-Stereo-D0 (a parser inconsistency that let this same information silently become visible-but-wrong after a canonical round trip); D0 fixed the round-trip inconsistency (5,000/5,000 corpus-stable) without closing this recall gap. Quantification and fix tracked as follow-up EZ-A0 (diagnosis) / EZ-S1 (production) | `validation.md`, Canonical-Stereo-D0 PR, `scripts/ez_stash_gap_report.py` (EZ-S1 full-corpus before/after accounting) |

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

## Defect found and FIXED this round: implicit-H count lost across the Kekulé→aromatic boundary

**Aromaticity-derived output contamination — three distinct items, kept separate so a
future fix's causal signal is never contaminated by mixing with another:**

```
芳香族性由来の出力汚染:
  ✅ implicit_h(ピロール/ピリジンN の型損失): 修正済み(8d0b992/9c50c08)
     — ECFP4/canonical/InChI の 89分子を同時に解消、因果を介入で確認
  ⬜ aromatic_context 系: 41/1000 未解決、既知、次スコープ
  ✅ InChI 固有 13.4% 順序依存 → assign_cip() 自体の正誤バグと判明、修正済み
     — 「自己一貫性のみ」の予定だったが根本原因は自己一貫性の問題ですら
       なかった: RDKit の per-atom CIP と比較したところ全立体中心の
       76.22%が(順序に関係なく、初回パースの時点で)誤っていた。
       独立した2バグ: (1) 環開き結合の隣接順序ずれ、
       (2) CIP球展開の二重結合重複原子が片側にしか追加されていなかった。
       RDKit正解率 76.22%→96.83%(4163中132残)。残り132は着手・撤回済み:
       三重結合にも同型の重複ルールを拡張したところ正味で悪化(16件新規誤り
       vs 1件修正)、原因はルール不足ではなく比較器自体(cip_branch_spheres
       がBFS深さごとに多重集合でプールする近似比較で、真の再帰的CIP有向
       グラフ比較ではない)と判明、コミットせずrevert。芳香環(96件、最大
       グループ)は未着手、P中心(15件)・三重結合(1件)・その他(43件)も
       未着手(件数は凍結コーパス`validation/cip_label_corpus.jsonl`基準、
       current oracle準拠 — legacy oracle基準の暫定値「三重結合10件・P6件」
       は上書き)。legacy版RDKitとrdCIPLabeler(現行)は43件で不一致があり、
       現行オラクル基準では 96.30%(155/4186残)がより正確な数値。
       設計RFC `docs/cip_accurate_rfc.md`(実装コードなし、chematic-cip
       クレート新設・真の再帰的digraph比較器・4マイルストーン計画)を
       この回で追加。InChI固有の順序超過は 74%→0% に解消(残り3.5%は
       canonical_smiles と共有の既知ベースライン順序依存で説明可能)。
       詳細は下記セクション。
```

Measuring ECFP4 surfaced a real correctness issue that took two rounds to fully
root-cause and fix — this is **Layer 2 only** (see the worked example below for the full
Layer 1/Layer 2 split; Layer 1, chematic's ECFP4 including aromaticity in its invariant,
is a design choice, not a bug, and is unrelated to what follows).

**Round 1 (measurement)**: `ecfp4()` (and `ecfp6()`/`fcfp4()`, which share the same
invariant seed function) was representation-dependent for ~13% of molecules even after
following the documented `apply_aromaticity()` contract. A cross-consumer check
(`scripts/ecfp4_agreement.py` tier 6) found `canonical_smiles()` and InChI diverging on
the *exact same* molecules, post-`apply_aromaticity()` (n=1000: 130/132/134 residual
mismatches, 100% overlap both ways) — a single shared root cause across ≥3 core output
functions, not an ECFP4-specific quirk. Of that residual, 41/130 was attributable to a
genuine aromatic-atom/bond flag-assignment disagreement (same *class* as the pre-existing
`aromatic_context` limitation, not confirmed the same code path). The other 89/130 had an
**identical** flag-assignment multiset yet still diverged, with SSSR ring decomposition
and (for ECFP4) atom-order-dependence both ruled out as explanations — mechanism left
unidentified, two unranked candidates.

**Round 2 (root cause + fix, `scripts/aromaticity_mechanism_probe.py`, new)**: built a
genuine per-*physical*-atom/bond correspondence between the two representations using
RDKit's canonical atom ranking as an independent oracle (chematic parses a given SMILES
string in the same atom order RDKit does — verified directly, not assumed — so rank-based
correspondence gives a true physical-atom mapping, unlike raw index comparison, which is
meaningless across differently-ordered respellings). **Positive control** (run first, per
this project's standing discipline): the correspondence tool correctly detected a
difference on 10/10 molecules from the already-known 41-set before any result on the
89-set was trusted. Result on the 89: **0 atom-flag diffs, 0 bond diffs, 0 SSSR-membership
diffs under true correspondence** — both candidates from Round 1 were false leads.

The actual divergence, found via a full 7-field `atom_table` correspondence diff: **the
`implicit_h` count** (not `aromatic`, not `bond.order`) differed for a specific
ring atom in **89/89** cases — 100% explained, no residual "mystery" cases. Root cause,
confirmed by reading `crates/chematic-perception/src/aromaticity.rs` and
`crates/chematic-core/src/valence.rs`:

- `apply_aromaticity_ex()` normalizes every aromatic-model ring bond to
  `BondOrder::Aromatic` uniformly, discarding the original Kekulé Single/Double pattern
  (confirmed by source read, not inferred).
- `implicit_hcount()`'s aromatic-path heuristic (`floor(n_aromatic_bonds × 1.5)`) is
  *correct* for its documented input — direct aromatic-written SMILES, where a
  lone-pair-donating "pyrrole-type" heteroatom is always written with an explicit bracket
  H (`[nH]`) and a "pyridine-type" one never is (OpenSMILES convention). But once
  `apply_aromaticity_ex()` erases the Single/Double pattern, a pyrrole-type N (2 ring
  single bonds pre-normalization, needs 1 H) and a pyridine-type N (1 single + 1 double,
  needs 0 H) become **byte-identical** post-normalization (aromatic, 2 aromatic-order ring
  bonds, no bracket H, neutral) — there is no local signal left to tell them apart, so the
  heuristic silently returns 0 for both. An atom reached via direct aromatic-written
  SMILES gets the right answer for free (its bracket H count is read directly, bypassing
  the heuristic entirely); the identical atom reached via Kekulé-then-perceive does not.

**Fix** (`crates/chematic-perception/src/aromaticity.rs`, `apply_aromaticity_ex`):
capture `implicit_hcount()`'s value on each non-bracket atom **before** bond-order
normalization (while the Kekulé pattern is still intact and the heuristic still gives the
correct answer), then compare it against what the heuristic would compute *after*
normalization; only atoms where the two disagree get an explicit `hydrogen_count` frozen
onto the output atom. This is fully surgical: atoms where pre- and post-normalization
already agree (ordinary aromatic CH, pyridine-type N) are untouched — no spurious bracket
notation introduced anywhere. Verified by two new Rust unit tests
(`test_apply_aromaticity_preserves_pyrrole_nh_implicit_hydrogen`,
`test_apply_aromaticity_does_not_add_h_to_pyridine_type_n`) plus the full existing
`cargo test --workspace` suite (all pre-existing tests, 671+, still pass unchanged —
notably the *existing* `pyrrole_kekule()` test fixture never caught this because it
manually set `hydrogen_count` on construction, sidestepping the exact gap between
synthetic test fixtures and real SMILES-parser output that let this slip through 15
rounds).

**Causal verification, not just correlation** (the standard this round's plan explicitly
required): re-running `scripts/ecfp4_agreement.py` post-fix shows the residual drop from
130/132/134 to **41/43/45 simultaneously across all three consumers** — each losing
*exactly* the 89 molecules that had matching flags/bonds/rings, leaving the
pre-existing, distinct 41-set (real flag-assignment disagreement, untouched by this fix)
exactly where it was. `scripts/aromaticity_mechanism_probe.py`'s 89-set is now empty
(0/1000), while its positive control still correctly detects all 10/10 of the remaining
41-set. Fixing one root cause zeroed the identical 89 molecules in three independent
output functions at once — the intervention-based proof that this was genuinely one
shared defect, not three coincidentally-similar ones.

**Ground-truth verified at full-corpus scale, not just internally consistent, and not
just a handful of spot checks** (agreement between the two chematic representations
alone would not rule out "both now agree on the same wrong answer" — an early check that
compared `implicit_h` only where chematic's *own* two representations still disagreed was
vacuous, 0-of-0, since the fix already made them agree by construction; caught before it
was mistaken for the real evidence). The decisive check: `mol.mw` vs RDKit's
`Descriptors.MolWt()` for the **full 5,000-molecule corpus**, on the exact Kekulé-origin
→ `apply_aromaticity()` path the fix touches — **4998/5000 (99.96%) match RDKit's MW
exactly** (tolerance 0.02 Da). The 2 remaining mismatches are pre-existing and unrelated
to this fix (an isotope-labeled `[3H]` atom and a `[Te]`-containing molecule, both
already off by the identical amount in a *baseline* run with no Kekulé round-trip at
all — a standard-atomic-weight-table gap, not an aromaticity bug; logged as a new,
separate, low-priority blank cell rather than investigated further this round). No
regression on the ~4,910 molecules unaffected by the original bug. **Severity was larger
than the fingerprint/canonicalization framing suggested**: `implicit_hcount()` feeds
molecular weight and formula directly, so both were silently wrong (short ~1.008 Da) for
any pyrrole-type NH heterocycle (imidazole, indole, pyrrole, purine substructures)
reached via the extremely common "parse Kekulé SDF/MOL input → `apply_aromaticity()` →
compute descriptors" pipeline — a core-property correctness bug, not merely a
fingerprint/canonicalization quirk.

**Remaining, explicitly out of scope for this fix**: the 41-molecule genuine
flag-assignment disagreement (same class as `aromatic_context`, azulene/purine) is
unchanged and still open. `canonical_smiles`'s (3.8%) and InChI's (13.4%, confirmed
disjoint via molecule-set intersection, not assumed) own separate order-sensitivity
defects are also unchanged — neither was ever part of this mechanism, both remain
unattributed blank cells for a future round.

---

## Defect found and FIXED this round: `assign_cip()` R/S was wrong on ~24% of stereocenters, not just order-unstable

**Scope changed mid-investigation, and the change is the finding.** The plan going in was
narrow and explicitly bounded by the user: fix InChI's own 13.4% order-sensitivity
(previous round), self-consistency only, not full RDKit/standard-InChI compliance —
because chematic's InChI numbering was already known to be non-standard, so "does the
`/t`-layer parity sign stay the *same* across two respellings of the same molecule" was
treated as the only tractable, well-scoped target. Root-causing that instability
(`crates/chematic-chem/src/cip.rs::assign_tetrahedral`) found it wasn't a self-consistency
problem at all: `assign_cip()`'s R/S codes are absolute (unlike InChI's numbering, there is
no legitimate "chematic's own convention" for CIP R/S), and measuring them against RDKit's
per-atom `_CIPCode` oracle (`Chem.AssignStereochemistry(cleanIt=True, force=True)`) — for
the very first time this project has done so — showed **76.22% correct (3173/4163
stereocenters)** even on first-parse, non-respelled SMILES from the 5,000-molecule corpus.
The instability chasing InChI specifically was a symptom; the disease was upstream and
touched every consumer of `assign_cip`/`tetrahedral_stereo_neighbors`: pure-Rust InChI,
**`native-inchi`'s C-library conversion path** (`crates/chematic-inchi/src/native/convert.rs`
calls `tetrahedral_stereo_neighbors` directly — the exact "use `native-inchi` for real
InChI" escape hatch documented elsewhere in this file was *also* silently affected), and
any other R/S consumer. `canonical_smiles`'s own stereo output was unaffected — it reads
chirality through a separate, already-correct code path
(`crates/chematic-smiles/src/canonical.rs`, see below).

**Root cause 1 — wrong SMILES-encounter order, discarding an existing correct mechanism.**
`assign_tetrahedral()` rebuilt a chiral atom's 4 substituents from raw
`Molecule::neighbors()` adjacency order plus a hand-rolled heuristic to place the bracket-H
slot. That heuristic assumed adjacency order always matches SMILES textual order, which is
false for exactly one case: a stereocenter that **opens** a ring before its other
substituents (`[C@@H]1N...`, ring digit written before the continuation atom). A
ring-*opening* bond's partner is unknown until the matching closing digit is reached later
in the string, so `Molecule::neighbors()` only materializes it then — by which point a
branch/continuation atom that appears later in the text, but has nothing to wait on, has
already been added. The result: the neighbor meant to occupy SMILES position 2 (by
written position) silently lands at position 3 and vice versa, flipping the parity
computed from it. Ring-*closing* bonds are unaffected (both endpoints already known,
materialized immediately) — this is why the bug is intermittent, not universal, and why it
surfaced as "order sensitivity" rather than a flat wrong-every-time bug: which respelling
of a given molecule happens to put a stereocenter's ring bond on the opening side is itself
representation-dependent.

The fix did not need to invent new machinery: the SMILES parser already builds the correct
sequence at parse time — `Molecule::stereo_neighbor_order` (`crates/chematic-core/src/molecule.rs`),
populated via a dedicated slot-based mechanism (`StereoEntry::PendingRing` in
`crates/chematic-smiles/src/parser.rs`) that resolves each ring-closure entry to its real
partner *without* disturbing its textual position in the sequence. `chematic-smiles`'s own
canonical-SMILES writer already reads this field
(`crates/chematic-smiles/src/canonical.rs:827`) to re-derive `@`/`@@` correctly for a new
output order — which is exactly why `canonical_smiles`'s stereo output was never affected by
this bug. `crates/chematic-chem/src/cip.rs` simply never read it. Fix: extracted a shared
`stereo_neighbors()` helper that prefers `Molecule::stereo_neighbor_order` (falling back to
the old adjacency-order heuristic only for molecules with no parse-time stereo data, e.g.
built directly via `MoleculeBuilder`), and pointed both `assign_tetrahedral()` and the
public `tetrahedral_stereo_neighbors()` (previously an independent duplicate of the same
buggy logic) at it.

**Root cause 2 — unmasked by fix 1, not caused by it: incomplete CIP double-bond
duplication.** Fixing the neighbor order alone, verified against the RDKit oracle, produced
824 newly-correct stereocenters but **15 new wrong ones** that had been correct before —
concerning enough to stop and re-verify rather than declare victory (per this project's
standing "verify at full-corpus scale, not spot checks" discipline). Direct trace of the
smallest example (`C=C(C)[C@@H]1CN[C@@H](C(=O)O)[C@@H]1CC(=O)O`, atom 3) showed the *new*
(correct) neighbor order was being combined with a rank order that was itself wrong: CIP
rule requires a double bond `A=B` to duplicate its partner into **both** atoms' own
substituent spheres (a phantom "B" in A's list *and* a phantom "A" in B's list), but
`cip_branch_spheres()` (`crates/chematic-chem/src/cip.rs`) only added the "arrival" side
(B's own sphere gets a phantom-A when B is *expanded*, having been reached via the double
bond) — never the "departure" side (A's own sphere, populated while iterating A's
neighbors, never got a second phantom entry for B). A vinyl/methyl-substituted stereocenter
neighbor (`C(=CH2)(CH3)-`) therefore scored a substituent set of `(C,C)` instead of the
correct `(C,C,C)`, silently losing a length-based tie-break it should have won. This half of
the duplication rule had been missing since the sphere-expansion code was first written; it
had simply never been exercised by a case where getting it right vs. wrong changed the
final R/S — until fix 1 stopped accidentally cancelling it out. Fix: add the missing
departure-side phantom (`is_double` check inside the neighbor-iteration loop, pushing an
extra copy of the child key into the same layer) — 5 lines, `BondOrder::Double` only
(aromatic/Mancude-ring duplication is a separate, harder problem, deliberately not
attempted here — see below).

**Causal verification.** The double-bond fix, applied on top of the order fix, introduced
**zero** new regressions relative to the order-fix-alone state (`post2 - post1` = ∅) and
resolved 6/15 of the order-fix's regressions outright (the other 9 all involve an aromatic
ring directly bonded to the stereocenter — the deliberately out-of-scope Mancude case, see
below). Net, full-corpus, both fixes together: **990 → 132 mismatches against RDKit's
per-atom CIP oracle (76.22% → 96.83%, 4163 stereocenters)**; **9 residual regressions
remain** relative to the original 990-mismatch baseline, all attributable to the
now-precisely-characterized Mancude gap, none to either fix newly introduced. Cross-checked
against the metric this investigation started from: InChI's own order-only mismatch rate
dropped from 13.4% to 3.5% (`scripts/ecfp4_agreement.py` tier 6, n=1000), and its
InChI-*specific* excess over `canonical_smiles`'s own baseline order-sensitivity —
previously 74% (99/134) — is now **0%**: the residual 3.5% is fully explained by the same
shared ordering-layer churn `canonical_smiles` already has, no separate InChI-specific
mechanism detected at this sample size. Full workspace test suite (673+ tests across 18
crates), `native-inchi` feature tests, and `bash scripts/check.sh` all pass; 2 new
regression tests added (`test_tetrahedral_stable_when_ring_bond_opens_before_other_neighbors`,
`test_tetrahedral_double_bond_duplicates_into_own_sphere`, both in `crates/chematic-chem/src/cip.rs`).

**Remaining, explicitly out of scope for this fix**: 132/4163 stereocenters (3.2% vs. the
legacy RDKit oracle; ~155/4186, 3.7%, vs. the more authoritative `rdCIPLabeler` modern
oracle — see below) still disagree. **Attempted and reverted, not merely deferred**: a
follow-up round categorized the residual as (at least) four distinct mechanisms — 96
aromatic-ring-adjacent, 15 phosphorus stereocenters, 1 triple-bond-adjacent, 43
uncharacterized (counts per the frozen `validation/cip_label_corpus.jsonl` corpus,
measured against the modern oracle; an earlier same-round estimate against the legacy
oracle had said "10 triple-bond, 6 phosphorus" — superseded by the frozen, oracle-pinned
numbers here) — and extended the double-bond duplication fix above to
`BondOrder::Triple` (2 phantom duplicates per side instead of 1) as the cleanest
mechanical analog. Result: **net negative**, 16 newly-wrong stereocenters vs. 1 newly-fixed
(traced by hand; a stereocenter directly bonded to an alkyne outranked a structurally
richer ring branch purely because triple-bond duplication concentrates 3 atoms into one
BFS shell while the ring's richness is spread across deeper shells). Root cause:
`cip_branch_spheres`/`compare_branches` pools *all* atoms at each BFS depth into one
sorted multiset and compares shell-by-shell — not the true CIP hierarchical-digraph
algorithm, which recursively compares branch-by-branch, following the highest-priority
sub-branch first. Adding a locally-correct duplication rule to this approximate comparator
is whack-a-mole: the double-bond fix above netted positive by (this corpus's) luck of
distribution, the triple-bond attempt netted negative — neither outcome validates the
comparator itself. **Reverted, not committed.** Properly closing the residual requires
replacing the shell-multiset comparator with a true recursive branch-by-branch traversal —
a materially larger undertaking (RDKit's own `CIPLabeler` is thousands of lines), not a
rule addition, and a scope decision for a future round. The aromatic bucket (96, largest)
carries an *additional* risk on top of the comparator issue: `kekulize()`'s matching is
atom-order-dependent, so Kekulize-then-duplicate could reintroduce the exact order-sensitivity
the fix above just closed for heteroaromatic rings — not gated/traced this round. Also
found while re-checking the oracle: legacy `Chem.AssignStereochemistry`/`_CIPCode` and
modern `rdCIPLabeler.AssignCIPLabels` disagree on 43 atom-cases in this corpus (the legacy
oracle has its own non-trivial error floor) — of the residual, 123 are wrong under *both*
oracles (solid floor of real chematic bugs), 9 were legacy-oracle errors (chematic was
actually right), 32 are new-only-under-the-modern-oracle. The pre-existing
`aromatic_context` 41/1000 flag-assignment defect is unrelated and also still open.

**Design RFC for a proper fix, this round**: `docs/cip_accurate_rfc.md` — a design-only
(no engine code) proposal for a separate "Accurate CIP" engine (new `chematic-cip` crate,
true recursive branch-by-branch hierarchical-digraph comparator, dual-mode API alongside
the existing fast/approximate engine, 4 milestones gated at 98%/99%/99.5% modern-oracle
agreement). The 155-case residual is frozen, reproducible, and oracle-pinned at
`validation/cip_label_corpus.jsonl` (manifest line records RDKit/chematic versions;
regenerate via `scripts/cip_ground_truth.py --freeze`) so future milestone PRs have a
fixed baseline to diff their engine's output against.

Reproduce: `.venv/bin/python scripts/cip_ground_truth.py ~/Downloads/SMILES.csv --freeze`
(RDKit CIP oracle comparison + frozen corpus, new this round); `.venv/bin/python scripts/ecfp4_agreement.py
~/Downloads/SMILES.csv --layer2-sample 1000` (InChI order-only / InChI-specific-excess
numbers, tier 6).

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
| After calling `apply_aromaticity()` as documented (**pre-fix, historical**) | **130/1000 (13.0%)** *still* mismatch — now **41/1000 (4.1%)** post-fix, see below |
| Of that residual: still explained by aromaticity-perception disagreement between the two spellings | **41/130 (~32%)** — same *class* of bug as the known `aromatic_context` regression; literal same-code-path attribution unverified (see below); **still open, unaffected by this round's fix** |
| Of that residual: **not** explained by aromaticity perception at any granularity checked | **89/130 (~68%)** — confirmed via 3 independent checks; root cause found and **FIXED** a follow-up round (`implicit_h` count, not perception — see the mechanism note below) |

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

**Mechanism of the shared residual — RESOLVED in a follow-up round, see the
["Defect found and FIXED"](#defect-found-and-fixed-this-round-implicit-h-count-lost-across-the-kekuléaromatic-boundary)
section above for the full writeup.** At the time this worked example was first written,
per-atom aromatic-flag assignment, per-bond *type* assignment, and SSSR ring topology were
all ruled out for 89/130, atom-order-dependence was ruled out for ECFP4's share of all
130, and a one-molecule spot check stalled on atom-index correspondence (the Kekulé-origin
respelling has a different internal atom order, so raw per-index atom-table comparison
isn't meaningful without a real correspondence step). The two candidates floated at the
time — (a) a positional bond-structure difference below the multiset's resolution, (b) a
symmetric aromaticity-flag swap — were **both wrong**: building the missing
correspondence step (RDKit canonical-rank-based atom mapping,
`scripts/aromaticity_mechanism_probe.py`) found 0 flag/bond/ring-membership diffs for all
89, and a full 7-field `atom_table` diff under that same correspondence found the actual
answer: **`implicit_h` count**, differing in 89/89 cases — a variable neither candidate
nor the original multiset check had ever compared. Now fixed; the shared residual is
41/1000 (the pre-existing, distinct `aromatic_context`-class defect only). InChI's own
13.4% order-sensitivity defect was a **separate, untouched blank cell** at the time this
paragraph was written — since fixed (13.4%→3.5%, InChI-specific excess 74%→0%) in a later
round, not part of this mechanism; see
["`assign_cip()` was wrong on ~24% of stereocenters"](#defect-found-and-fixed-this-round-assign_cip-was-wrong-on-24-of-stereocenters-not-just-order-unstable)
below.

**Blast radius beyond ECFP4**: `initial_atom_id` (the invariant seed carrying the
aromaticity byte) is shared by ECFP6 and FCFP4, and a 300-mol spot check confirms both
inherit the same naive representation-dependence (ECFP6 94.0%, FCFP4 94.0%, vs. ECFP4's
92.2%) — this is a shared-seed-function defect, not ECFP4-specific. Not separately
measured post-`apply_aromaticity()` for ECFP6/FCFP4.

The 92.2% naive case is arguably working as documented — CLAUDE.md/README already state
Kekulé input needs `apply_aromaticity()` first, and this is the first time that contract
has been checked specifically against fingerprint output. The (pre-fix) 13.0% residual was
not excused by that contract: a caller following the documented procedure still got a
different fingerprint for the same molecule — and, per the tier 6 findings above, the same
molecule also canonicalized and InChI-generated differently. About a third of that residual
matches the same *class* of bug as an already-known, already-deferred limitation (still
open); about two-thirds was the `implicit_h`-loss bug, now fixed — see
[Defect found and FIXED this round](#defect-found-and-fixed-this-round-implicit-h-count-lost-across-the-kekuléaromatic-boundary)
above.

Reproduce: `.venv/bin/python scripts/ecfp4_agreement.py ~/Downloads/SMILES.csv --aromaticity-sample 1000 --layer2-sample 1000`
and `.venv/bin/python scripts/aromaticity_mechanism_probe.py ~/Downloads/SMILES.csv --sample-n 1000`
