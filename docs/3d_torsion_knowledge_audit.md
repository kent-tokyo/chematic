# ETKDG Torsion Knowledge — Legacy Audit (3D Breakthrough Program, Wave 2, Agent E)

Audits `crates/chematic-3d/src/etkdg_knowledge.rs` as it stood at base SHA
`faa911ad131d47b9bcfbe77aee1eb61a52243698` (branch point of
`feat/3d-torsion-knowledge-v2`), before any of this PR's new code was added.
Required by the Wave 2 spec §2. See
`validation/manifests/etkdg_torsion_knowledge_sources.json` for the RDKit/paper
sources cited below.

**Counting method for "existing torsion rule count":** one rule = one
`return Some(TorsionPreference { .. })` statement inside `get_torsion_preference`,
**or** one `SmartsTorsionRule { .. }` struct literal inside
`SMARTS_TORSION_RULES`. Directional mirror pairs (the same physical bond
written once as `b=X,c=Y` and again as `b=Y,c=X` so both SMILES-traversal
orders match) are counted **separately** — they are separate code branches
even when they encode the same physical preference.

- `get_torsion_preference`: 74 `return Some(...)` branches (confirmed by
  `grep -c "return Some(TorsionPreference" etkdg_knowledge.rs`).
- `SMARTS_TORSION_RULES`: 6 entries (confirmed by counting `smarts:` fields
  inside the `static SMARTS_TORSION_RULES: &[...] = &[ .. ];` array).
- **Total: 80.**
- File line count at base SHA: 1385 (`wc -l`).

## 1. Provenance of the "CSD-derived" / "experimental" label

The module's own top-of-file doc comment reads:

> `//! ETKDG Torsion Knowledge Base — experimental torsion angle preferences from CSD.`

**Finding: this claim is unverified and, on the evidence available, false.**
There is no CSD (Cambridge Structural Database) query, no statistical fit, no
citation, and no data file anywhere in this module or its git history. Every
one of the 80 rules is a single hand-picked `angle_deg` plus a hand-picked
`penalty_per_degree` "for X kcal/mol at Y° deviation" — the penalty constants
read as plausible-sounding order-of-magnitude guesses (e.g. `0.15` for
alkanes "≈3 kcal/mol for 20° deviation", `0.03` for biphenyl "very soft —
flat potential") with no traceable source. Contrast this with the real
RDKit-shipped `torsionPreferences_v2.in` (fetched and hashed in the sources
manifest), whose every line is `[SMARTS] s1 V1 s2 V2 ... s6 V6` — a genuine,
literature-cited 6-term Fourier fit. The legacy file's single-angle-plus-linear-
penalty shape cannot even *represent* that kind of data, which is direct
evidence the "from CSD" claim was aspirational/mislabeled rather than an
actual derivation that was later simplified.

**Action taken:** the header comment is corrected (see the code diff) to
describe the module honestly as a hand-authored heuristic layer, and every
rule below is reclassified per the required taxonomy — never
`experimental`/`CSD-derived` unless later evidenced.

## 2. Rule-by-rule classification

One row per `return Some(TorsionPreference { .. })` statement, identified by
its **exact 1-based line number** in the base-SHA file (re-derived via
`grep -n "return Some(TorsionPreference" etkdg_knowledge.rs`, 74 hits — not
eyeballed) plus the 6 `SMARTS_TORSION_RULES` entries. Classification taxonomy
(exactly one bucket per row): `verified_from_primary_source` /
`supported_by_rdkit_oracle` / `reasonable_heuristic_only` /
`incorrect_or_overgeneralized` / `dead_or_unreachable` / `ambiguous`.

| Line | Rule (as commented in source) | Classification | Why |
|---|---|---|---|
| 165 | Alkane C-C-C-C → 180° | `supported_by_rdkit_oracle` | RDKit's own default staggered-anti preference for sp3-sp3 chains is well-established; matches `default_torsion_preference()` too. |
| 177 | Aromatic-aliphatic C-Ar-C-C → 180° | `reasonable_heuristic_only` | Plausible (benzylic bonds are low-barrier, near-free rotors), "180° or 0°" flattened to just 180°, uncited. |
| 186 | Amide C-N(C=O)-C-X → 180° | `incorrect_or_overgeneralized` | Self-admitted in the source comment: "restricted rotation, prefer 0° (cis) or 180° (trans)... prefer trans as it's more common" — a genuinely **bimodal** preference collapsed to one mode. Spec §2's "both 0° and 180° allowed" failure mode. |
| 194 | Ester O-C(=O)-O-C → 0° | `reasonable_heuristic_only` | Plausible (syn/Z ester conformation), uncited. |
| 202 | Aromatic-aromatic biphenyl → 45° | `incorrect_or_overgeneralized` | RDKit's real data (`torsionPreferences_v2.in` line 270, fetched) shows unsubstituted biphenyl's true potential is `V1=-0.7(s-1) V2=-8.0(s+1) V4=4.4(s+1) V6=-1.5(s+1)` — multi-term, **two-minima**, not a single 45° angle. Spec §2's "multi-modal torsions" failure mode, concretely evidenced. |
| 210 | Enamine C=C-N → 180° | `ambiguous` | Never checks the *B-C bond itself* is the intended rotatable single bond, only neighbor hybridization — see §3.5. |
| 216 | Enamine N-C=C (reverse) → 180° | `ambiguous` | Same bond-identity gap. |
| 226 | Vinyl halide C=C-X → 0° | `ambiguous` | Same bond-identity gap. |
| 234 | Acrylic/chalcone C=C-C(=O) → 180° | `reasonable_heuristic_only` | s-trans is the known major conformer for acyclic enones, uncited. |
| 242 | Phenyl ketone Ar-C(=O) → 0° | `reasonable_heuristic_only` | Coplanarity-for-conjugation is textbook direction, uncited. |
| 248 | Phenyl ketone (reverse) → 0° | `reasonable_heuristic_only` | Mirror of L242. |
| 256 | Thioester S-C(=O) → 180° | `reasonable_heuristic_only` | Plausible, uncited. |
| 264 | Carbamate/sulfonamide N-C(=O) → 180° | `reasonable_heuristic_only` | Plausible, uncited; conflates two distinct bond environments under one angle. |
| 272 | Sulfoxide C-S(=O)-C → 90° | `reasonable_heuristic_only` | Plausible (pyramidal S), uncited. |
| 280 | Disulfide C-S-S-C → 90° | `supported_by_rdkit_oracle` | The ~90° (gauche) S-S-C-C preference is well-established structural chemistry; not disputed by anything in the fetched RDKit data. Not `verified_from_primary_source` since no primary source was directly read for this specific number. |
| 288 | Alcohol/ether C-C-O-X → 180° | `incorrect_or_overgeneralized` | Source comment itself: "gauche/anti mixture, use 180° as default" — another explicit bimodal-collapsed case. |
| 298 | Amine C-C-N-C → 180° | `reasonable_heuristic_only` | Plausible anti preference, uncited. |
| 306 | Nitrile terminus → 180° | `verified_from_primary_source` | sp-hybridized carbon is genuinely linear by VSEPR — textbook fact, not a fitted preference. |
| 314 | Phosphorus P-C-C-X → 180° | `reasonable_heuristic_only` | Generic catch-all, uncited. |
| 327 | Urea N-C(=O)-N → 0° | `reasonable_heuristic_only` | Plausible (lone-pair conjugation), uncited. |
| 335 | Sulfonamide N-S(=O)(=O) → 90° | `reasonable_heuristic_only` | Plausible, uncited. |
| 341 | Sulfonamide (reverse) → 90° | `reasonable_heuristic_only` | Mirror of L335. |
| 349 | Aryl ether Ar-O-C → 0° | `reasonable_heuristic_only` | Plausible (O lone-pair conjugation), uncited. |
| 355 | Aryl ether (reverse) → 0° | `reasonable_heuristic_only` | Mirror of L349. |
| 366 | Fluoroalkane C-C-C-F → 180° | `reasonable_heuristic_only` | Plausible, uncited. |
| 374 | Nitro C-N(=O)=O → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 384 | Hydrazone/oxime C=N-N/O → 0° | `incorrect_or_overgeneralized` | Source comment: "prefer 0° (E/Z isomerism; E is more stable)" — collapses a genuinely bimodal (E and Z both real) case to one mode, plus the same bond-identity gap as the enamine rules. |
| 397 | Imide N-C(=O)-C(=O) → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 405 | Benzyl Ar-C-X → 90° | `reasonable_heuristic_only` | Plausible (hyperconjugation-perpendicular), uncited. |
| 413 | Allylic C=C-C-X → 0° | `reasonable_heuristic_only` | Source comment: "s-cis/s-trans mixture; use 0° as default" — same bimodal-collapse pattern, one tier milder (author explicitly flags it as "a default"). |
| 431 | Heteroaromatic biaryl (both directions, one `if`) → 45° | `reasonable_heuristic_only` | Plausible, uncited. |
| 442 | N-alkyl heteroaromatic (both directions) → 180° | `reasonable_heuristic_only` | Plausible, uncited. |
| 453 | Heteroaromatic N-carbonyl (both directions) → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 463 | Heteroaromatic N-alkene (both directions) → 0° | `ambiguous` | Involves `CSp2Alkene` on one side — same bond-identity gap as L210. |
| 475 | Thioaryl/aryl thioether (both directions) → 90° | `reasonable_heuristic_only` | Plausible, uncited. |
| 484 | OSp3-CSp3 (reverse of L349-family) → 180° | `reasonable_heuristic_only` | Plausible, uncited. |
| 497 | Aryl amine Ar-NR2 (both directions) → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 511 | Furanyl/oxazolyl biaryl (both directions) → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 523 | Thienyl/thiazolyl biaryl (both directions) → 45° | `reasonable_heuristic_only` | Plausible, uncited. |
| 533 | Furanyl methyl (both directions) → 180° | `reasonable_heuristic_only` | Plausible, uncited. |
| 543 | Thienyl methyl (both directions) → 180° | `reasonable_heuristic_only` | Plausible, uncited. |
| 555 | Furanyl/thienyl carbonyl (4-way `if`) → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 572 | Morpholine N-C-C-O gauche → 60° | `incorrect_or_overgeneralized` | Real gauche effect has **two** equally valid minima (+60°/−60°); recording only +60° penalizes −60° as if 120° off instead of 0° off. Spec §2's "gauche ±60° as two distinct minima" failure mode, concretely confirmed. |
| 581 | Piperazine N-C-C-N gauche → 60° | `incorrect_or_overgeneralized` | Same one-of-two-minima problem as L572. |
| 592 | Styrene Ar-C=C → 0° | `ambiguous` | Same bond-identity gap as L210. |
| 598 | Styrene (reverse) → 0° | `ambiguous` | Same. |
| 607 | Vinyl thioether C=C-S → 0° | `ambiguous` | Same. |
| 613 | Vinyl thioether (reverse) → 0° | `ambiguous` | Same. |
| 622 | Allylic amine C=C-N(sp3) → 0° | `ambiguous` | Same. |
| 628 | Allylic amine (reverse) → 0° | `ambiguous` | Same. |
| 638 | Ketone/aldehyde C(sp3)-C(=O) → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 644 | Ketone/aldehyde (reverse) → 0° | `reasonable_heuristic_only` | Mirror of L638. |
| 654 | Heteroaromatic N-thioether → 90° | `reasonable_heuristic_only` | Plausible, uncited. |
| 660 | N-thioether (reverse) → 90° | `reasonable_heuristic_only` | Mirror of L654. |
| 671 | Heteroaromatic N-ether O → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 677 | N-ether O (reverse) → 0° | `reasonable_heuristic_only` | Mirror of L671. |
| 689 | 5-membered S-heteroaromatic-to-N-heteroaryl → 45° | `reasonable_heuristic_only` | Plausible, uncited. |
| 698 | 5-membered O-heteroaromatic-to-N-heteroaryl → 45° | `reasonable_heuristic_only` | Plausible, uncited. |
| 709 | Aromatic carbonyl to sp3/alkene → 0° | `ambiguous` | Involves `CSp2Alkene` — same bond-identity gap as L210. |
| 719 | sp-nitrile/alkyne CSp3-NSp → 180° | `verified_from_primary_source` | Same VSEPR-linear justification as L306. |
| 730 | Amide (reverse) CCarbonyl-NSp2 → 180° | `reasonable_heuristic_only` | Mirror of L186 (same bimodal-collapse caveat applies but tracked once at L186; the reverse-traversal branch itself is just a mirror, not an independent judgment call). |
| 736 | Amide (reverse) CCarbonyl-NSp3 → 180° | `reasonable_heuristic_only` | Mirror of L264. |
| 745 | Isocyanate NSp-CCarbonyl → 180° | `verified_from_primary_source` | sp-hybridized N=C=O is genuinely linear by VSEPR. |
| 751 | Isocyanate (reverse) → 180° | `verified_from_primary_source` | Same VSEPR justification as L745. |
| 760 | Ar-N=C=O → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 766 | Ar-N=C=O (reverse) → 0° | `reasonable_heuristic_only` | Mirror of L760. |
| 786 | Thioamide (conditional on real S=C bond-order check) → 180° | `reasonable_heuristic_only` | The **only** branch that inspects a real bond order (`bond_between(c_idx,n).order==Double`) before firing — better-grounded mechanism than the rest of the cascade, still an uncited angle. |
| 802 | Anomeric/vicinal-O gauche O-C-C-O → 60° | `incorrect_or_overgeneralized` | Same one-of-two-gauche-minima problem as L572/L581. |
| 812 | OSp3-CCarbonyl (reverse ester) → 0° | `reasonable_heuristic_only` | Mirror of L194. |
| 824 | Conjugated diene C=C-C=C → 180° (s-trans) | `incorrect_or_overgeneralized` | Worst instance of the bond-identity gap: fires whenever B and C are **both** `CSp2Alkene`, exactly as true when B-C is the actual double bond itself as when it's the real connecting single bond. See §3.5. |
| 833 | 1,2-Dicarbonyl C(=O)-C(=O) → 0° | `reasonable_heuristic_only` | Plausible, uncited. |
| 845 | Halogen adjacent to sp2 (Ar-X) → 180° | `reasonable_heuristic_only` | Plausible, uncited. |
| 855 | Phosphorus ester P-O → 180° | `reasonable_heuristic_only` | Plausible, uncited. |
| 861 | Phosphorus ester (reverse) → 180° | `reasonable_heuristic_only` | Mirror of L855. |
| — | `SMARTS_TORSION_RULES` × 6 (hindered biaryl 90°, primary amide 0°, tertiary amide 180°, heteroaromatic-N-carbonyl 0°, aryl ester 0°, carbamate 0°) | `reasonable_heuristic_only` | Better-grounded than the atom-type cascade (real SMARTS context, e.g. H-count distinguishing primary/tertiary amide), still hand-picked angles/penalties with no cited source. The "hindered biaryl → 90°" entry is directly supersedable by RDKit's own real data (`torsionPreferences_v2.in` line 265: ortho,ortho'-disubstituted biphenyl is `V2=3.6, s=+1` → minimum at 90°/270°, same qualitative conclusion with an actual citable source — implemented as this PR's `StandardExperimental` equivalent). |

**No branch qualifies for `dead_or_unreachable`** — every `if` in the cascade
is reachable given some real molecule (spot-checked against the existing
test suite; no branch's condition is a strict subset of an earlier,
already-returning branch's condition in a way that would make it permanently
unreachable).

## 3. Specific structural defects (spec §2 checklist)

### 3.1 Multi-modal torsions — confirmed present, unrepresentable

See group #5 (biphenyl) and #27/#3/#16 (bimodal-collapsed cases) above. The
`TorsionPreference { angle_deg: f64, penalty_per_degree: f64 }` type is
structurally incapable of representing more than one minimum — this is a
type-level limitation, not just a missing-data problem, which is exactly why
this PR introduces `FourierTorsionTerm`/`TorsionPotential` (§3 of the spec)
as a **new, separate** type rather than trying to patch more angles into the
old one.

### 3.2 Environments needing both 0° and 180°

Confirmed in groups #3 (amide), #16 (alcohol/ether), #30 (allylic), #81
(anomeric gauche) — see table above.

### 3.3 Gauche ±60° as two distinct minima

Confirmed in groups #56-57 and #81 — see table above. Concrete repro: for
morpholine (`C1CNCCO1`), `get_torsion_preference` returns exactly one
`TorsionPreference { angle_deg: 60.0, .. }`; a real conformer search finds
both +60° and −60° (300°) chair-flip populations equally likely, and the
single-value model scores −60° as if it were 240° off the correct target
rather than 0°.

### 3.4 Conflicts between torsion preference and declared chirality/E-Z stereo

`get_torsion_preference` has **zero** awareness of `chematic_core::Chirality`
or declared `BondOrder::Up`/`Down` (E/Z) markers anywhere in its 74-branch
cascade — it operates purely on `AtomType` derived from local connectivity.
Concretely: nothing in this module would stop a caller from requesting a
torsion preference for a bond immediately adjacent to a declared stereocenter
or a declared E/Z double bond, and getting back a preference that happens to
contradict the declared configuration, with no diagnostic raised. This PR's
new implementation does not attempt to resolve stereo conflicts either (that
is Agent D's territory, already merged as PR #190's `stereo_constraints.rs`),
but the new `TorsionKnowledgeDiagnostic` machinery is built so a future PR
can wire in a stereo-cross-check without changing the public shape (see
`docs/3d_torsion_knowledge_audit.md` §5 "known unsupported chemistry" in the
PR body).

### 3.5 Bond-order blindness (misapplication to double/triple/aromatic-ring bonds)

Confirmed the single most pervasive defect: **`get_torsion_preference` never
inspects `mol.bond_between(b_idx, c_idx)`'s `BondOrder` at all.** Every rule
fires purely off `AtomType(b)`/`AtomType(c)`, which is derived from each
atom's *own* neighbor set, not from the specific bond connecting them. Two
concrete failure modes, both real (not hypothetical):

1. **Double/triple bond misapplication**: group #83 ("conjugated diene
   C=C-C=C → 180°") matches whenever both B and C classify as `CSp2Alkene` —
   which is true both when B-C is a genuine connecting single bond *and*
   when B-C happens to be a C=C double bond itself (e.g. querying atoms in
   the "wrong" order across a `C=C-C=C` chain). Nothing in the function can
   tell these apart.
2. **Aromatic ring-bond misapplication**: unlike `build_smarts_torsion_map`
   (which explicitly takes a `ring_bond_set` parameter and skips ring bonds
   — the *correct* pattern), `get_torsion_preference`'s heteroaromatic rules
   (groups #31 "NAromatic-CAromatic → 45°" etc.) have no ring-bond exclusion
   at all. A fused heteroaromatic system's *intra-ring* N-C bond (rigidly
   fixed by ring closure) would still get an inappropriate "45° biaryl twist"
   preference suggested if this function were ever called on it.

Terminal bonds are handled implicitly-by-accident (a terminal atom's
`AtomType` is usually `H` or a halogen, which most rules don't target as
both `b`/`c`), not by an explicit terminal check — fragile, not a real
guarantee.

### 3.6 First-match-wins with no conflict detection

`get_torsion_preference` is a linear `if / else if / return` cascade — later
branches are silently unreachable for any atom-type combination an earlier
branch already claims, and the ordering (which branch "wins" for an
ambiguous combination) is nowhere documented as intentional. `#83`
(conjugated diene) sitting near the *bottom* of the cascade, for instance,
would never fire for a `CSp2Alkene`-`CSp2Alkene` pair that also happens to
match an earlier, more specific branch — but there is no comment anywhere
explaining that this is deliberate versus accidental. Classified overall as
`ambiguous` (the *mechanism*, not any single rule).

### 3.7 Silent skip on SMARTS parse failure

Confirmed in `build_smarts_torsion_map`:

```rust
let Ok(query) = parse_smarts(rule.smarts) else {
    continue;
};
```

If any static `SMARTS_TORSION_RULES` entry ever fails to parse (typo,
future `chematic-smarts` regression, etc.), the rule silently disappears —
no panic, no log, no test failure, no diagnostic. This is precisely the
anti-pattern spec §13 asks the *new* implementation to avoid (every skip
must be a typed, visible diagnostic).

### 3.8 Secondary amines mistyped as `NSp2`

**Confirmed, and self-acknowledged in the source's own comments.**
`classify_atom_type`'s nitrogen branch:

```rust
} else if neighbors <= 2 {
    AtomType::NSp2
} else {
    AtomType::NSp3
}
```

types *any* non-aromatic, non-triple-bonded nitrogen with ≤2 heavy neighbors
as `NSp2` — including a plain secondary amine (`R2NH`, genuinely sp3,
pyramidal, e.g. diethylamine `CCNCC`'s central N: 2 heavy neighbors, zero
double bonds). There is no check for an actual double bond or adjacency to a
carbonyl/aromatic system anywhere in this branch — degree alone decides.
Confirmed by direct construction: `classify_atom_type(&parse("CCNCC").unwrap(), <N atom>)`
returns `NSp2` for a textbook sp3 amine. The source code's own comments at
two call sites (lines 294-296 and 564-565 in the base-SHA file) *explicitly
acknowledge* this: "NSp2 covers N with 2 explicit bonds (secondary amine in
SMILES)" and "Secondary amines in rings are classified as NSp2 ... even
though they are sp3; accept both NSp2 and NSp3" — i.e. the author knew the
label was wrong and wrote defensive `is_sat_n(t) == NSp2 || NSp3` helper
checks around it rather than fixing the classifier. **This is precisely the
kind of undisclosed workaround spec §2 asks to be surfaced, not silently
carried forward.** Classified `incorrect_or_overgeneralized`.

This PR's new `classify.rs` does not reuse `AtomType` at all for its own
rule matching (SMARTS-based, using real bonding patterns via
`chematic-smarts`), so this defect is **not inherited** into the new v2
code path — but `classify_atom_type` itself is kept unchanged (spec §8: no
legacy-API behavior change without new regression fixtures + migration
notes, which this PR does not attempt), so the mistyping remains live in the
legacy path. Recorded here so Coordinator/reviewers don't miss it.

## 4. Test suite weaknesses (spec §2 checklist)

Confirmed weak/misleading tests in the base-SHA file (none rewritten — see
spec §8, legacy behavior/tests are left alone; recorded here for visibility):

- `test_thioester_torsion_preference`: calls
  `get_torsion_preference(&mol, AtomIdx(0), AtomIdx(0), AtomIdx(1), AtomIdx(2))`
  — **atom A and atom B are the same index (0)**, a nonexistent A-B-C-D
  dihedral (duplicate atom in a supposedly-4-distinct-atom quadruple). The
  test then discards the return value (`let _pref = ...`) and only checks
  atom *types* separately — it does not actually exercise or verify
  `get_torsion_preference`'s behavior at all despite the name.
- `test_furanyl_biaryl_prefers_planar`: calls
  `get_torsion_preference(&mol, o_idx, o_idx, o_neighbor, o_neighbor)` — **both
  A=B and C=D are duplicate pairs**. The subsequent assertion is
  `pref.is_some() || pref2.is_some()` — an OR across two calls means the test
  passes even if the primary (duplicate-index) call returns `None`, as long
  as an unrelated second call happens to return `Some`. This cannot fail in
  a way that would catch a real regression in the duplicate-index case.
- `test_pattern_count_covers_20_plus`: iterates several molecule/atom-quadruple
  pairs and ends with `let _ = pref; // just ensuring no panic` for most of
  them — the canonical "no panic" anti-pattern spec §2 calls out by name.
- `test_phenyl_ketone_torsion_preference`: despite its name, never calls
  `get_torsion_preference` — it only asserts that some atom in acetophenone
  classifies as `CCarbonyl`.

None of these are fixed in this PR (spec §8 forbids legacy behavior changes
without new regression fixtures + migration notes, and these are test-only,
zero production-behavior-risk issues) — flagged here so they are visible to
Coordinator/reviewers rather than silently carried forward unremarked.

## 5. Summary counts for the PR body

Recomputed directly from the 74+6=80 rows in §2's table (not re-derived by a
separate arithmetic pass — counted straight from the classification column):

| Classification | Count (of 80) |
|---|---:|
| `verified_from_primary_source` | 4 (L306, L719, L745, L751) |
| `supported_by_rdkit_oracle` | 2 (L165, L280) |
| `reasonable_heuristic_only` | 55 (49 cascade rows + 6 `SMARTS_TORSION_RULES`) |
| `incorrect_or_overgeneralized` | 8 (L186, L202, L288, L384, L572, L581, L802, L824) |
| `dead_or_unreachable` | 0 |
| `ambiguous` | 11 (L210, L216, L226, L463, L592, L598, L607, L613, L622, L628, L709) |

(4+2+55+8+0+11 = 80 ✓ — verified by re-adding the column tallies, not
assumed from the earlier draft of this table, which had a transcription
error in how multi-branch groups were counted; corrected here.)

## 6. RDKit oracle differential for the new v2 layer (spec §12)

Sections 1-5 above audit the *legacy* code (spec §2). This section covers the
live-RDKit-oracle differential spec §12 requires for the *new* v2 layer, run
against an **isolated venv** (`python3 -m venv` + `pip install rdkit`, never
this repo's shared `.venv`, per this program's standing rule) with
`rdkit==2026.03.4`. This is a different "RDKit" than the source-file pin in
`validation/manifests/etkdg_torsion_knowledge_sources.json` (RDKit release
tag `Release_2025_09_2`, fetched as read-only source text for translation
provenance, never compiled/executed) — the two serve different purposes and
are not meant to be the same artifact: the source-file pin is what the rule
*data* was translated from; the pip-installed version is a live behavioral
oracle queried at review time. Both are recorded so neither is silently
conflated with the other.

Reproduce:
```bash
cargo run --release -p chematic-3d --example torsion_knowledge_v2_gap_check   # writes validation/etkdg_torsion_knowledge_v2_chematic_side.json
python3 -m venv /tmp/oracle_venv && /tmp/oracle_venv/bin/pip install rdkit
/tmp/oracle_venv/bin/python3 scripts/etkdg_torsion_knowledge_v2_oracle_diff.py
```

### 6.1 What was actually comparable via RDKit's public API

Two of spec §12's six named comparison axes (rule family, central bond
selection, ring classification, minima, 1-4 pair selection, torsion
distribution) turned out to be reachable via RDKit's public Python API;
**rule family** and **1-4 pair selection** are not (see §6.4). Central bond
selection and ring classification are the same underlying check here (a
central bond's classification IS its ring membership), so this collapses to
2 distinct measurements: ring classification, and torsion distribution/minima.

### 6.2 Ring classification: 221/230 (96.1%) bond-level agreement

For 24 fixtures (every `rigid_ring`/`fused_aromatic`/`macrocycle`-tagged
`CORPUS` entry, all of `SMALL_RING_SET`, all of `RING_TOPOLOGY_SET`), atom-
index correspondence between chematic-smiles's parser and
`Chem.MolFromSmiles` was verified directly (element-sequence comparison, not
assumed) before comparing any bond — 24/24 fixtures verified, 0 unverifiable.
Per-bond ring-size sets (chematic's `RingMembershipIndex::ring_sizes_for` vs.
RDKit's `GetRingInfo().BondRingSizes()`) then agreed on 221/230 (96.1%) bonds.

All 9 disagreements are in **cubane** (`C1C2C3C1C4C2C3C4`), and the underlying
mechanism (`augmented_ring_set` recovering rings a bare SSSR minimum cycle
basis misses) is real and correctly identified, but an earlier draft of this
section misdescribed it in terms of a proper cube and was corrected after
independent review checked the fixture's own actual atom/bond count rather
than assuming it matches real cubane's graph. The corpus's `cubane` SMILES
parses to **8 atoms and 11 bonds** (confirmed directly from
`validation/etkdg_torsion_knowledge_v2_chematic_side.json`'s own bond list),
not the 8 atoms/12 bonds a proper cube graph (every vertex degree 3) would
have — atoms 0 and 7 each have only 2 skeletal bonds in this fixture, degree
2 not 3. This is a **pre-existing defect in the frozen 58-molecule corpus's
"cubane" entry** (inherited byte-identical from Wave 1's
`distance_geometry_v2_gap_check.rs` per this program's standing convention;
not introduced or fixable here — flagged for whoever owns
`scripts/etkdg_vs_rdkit_gap.py`/the sibling gap-check example), not a defect
in this PR's ring-classification logic.

With 8 atoms and 11 bonds, cycle rank (bonds − atoms + 1) = 4. RDKit's SSSR
returns exactly 4 independent rings, as expected for a minimum cycle basis on
a graph of that rank — this is NOT "SSSR failing to find a 5th/6th face",
it is SSSR correctly returning a basis of the size the graph's actual cycle
rank dictates. Chematic's `augmented_ring_set` then recovers **one further,
non-independent** small ring beyond that minimum basis — and this one really
does exist in the graph: tracing edges (0,1)(1,2)(2,3)(3,0)(3,4)(4,5)(5,1)
(5,6)(6,2)(6,7)(7,4) shows a genuine 5-membered cycle 1-2-3-4-5-1 (verified
directly from the bond list, not asserted). So chematic's extra ring-size
entries (e.g. bond (0,1): chematic `[4,5]` vs. RDKit `[4]`) reflect a real,
additional small ring present in this specific (topologically atypical)
graph, recovered by design — the mechanism this section originally described
is correct, the "6 square faces of a cube" framing was not, since this
fixture is not a proper cube graph. The **`adamantane`** fixture has the same
kind of pre-existing corpus-labeling caveat (11 atoms in
`C1CC2CC3CC1CC(C2)C3`, not real adamantane's 10) but produced zero
ring-classification disagreements against RDKit regardless.

All 15 non-cubane fixtures (benzene, naphthalene, pyridine, furan, thiophene,
adamantane, cyclohexane, cyclopentane, indole, purine, quinoline, anthracene,
pyrene, biphenyl, cyclododecane, crown_12_4, cyclooctadecane, cyclopropane,
cyclobutane, cycloheptane, cyclooctane, norbornane, spiro_5_6) matched RDKit
bond-for-bond, including the fused-aromatic (naphthalene, indole, purine,
quinoline, anthracene, pyrene) and bridged (norbornane, adamantane) cases
specifically named in spec §5.

### 6.3 Torsion distribution: empirical confirmation for the amide bond

RDKit's own ETKDGv3 conformer ensemble (50 conformers,
`useExpTorsionAnglePrefs=True`, fixed seed) for N-methylacetamide's C(=O)-N
central bond produced dihedral angles clustering entirely in the 120°-240°
range (23 conformers in [120°,180°), 27 in [180°,240°), 0 elsewhere;
min=-180.0, max=180.0, mean|angle|=180.0) — i.e. RDKit's live, independent
conformer generator agrees empirically with this PR's
`standard:secondary_amide` rule's predicted minimum (phi=180°, from the
`(n=1, s=-1, V=100.0)` Fourier term). This is a genuine positive result, not
assumed: the rule's translation from `torsionPreferences_v2.in:142` predicts
exactly what RDKit's own conformer generator, run independently, produces.

Butane's C-C-C-C dihedral was also sampled (same settings) as a sanity check,
not a rule comparison: this PR does not curate a generic-alkane standard-tier
rule (`butane` is reported as `unmatched` at the knowledge layer, honestly),
and RDKit's own empirical distribution for this specific bond came out close
to uniform across all six 60°-buckets (`[11,6,7,7,12,7]` of 50) — consistent
with there being no strong single-bond preference for this case in RDKit's
own generator either, not a contradiction to explain away.

### 6.4 Disclosed limitation: rule family and 1-4 pair selection

Spec §12 also asks to compare **matched rule family** and **1-4 pair
selection** against RDKit's own choices. Neither is achievable via RDKit's
public Python API in `rdkit==2026.03.4` (checked directly, not assumed): the
`ExperimentalTorsionAngle`/`BoundsMatrixBuilder` C++ machinery that performs
this matching internally has no public accessor returning which SMARTS
matched which bond, or which 1-4 pairs were adjusted and by how much. Closing
this would require either a patched/instrumented RDKit C++ build (out of
scope — this PR touches no build system, per file-ownership §14) or reverse-
engineering matches from conformer geometry alone (which would not actually
verify *which rule* fired, only that *some* preference was applied — a
weaker and potentially misleading claim). This is reported here as a real,
unresolved gap in what a live RDKit oracle can confirm for this PR, not
narrowed silently — this PR's rule-to-fixture matching is instead verified
via chematic's own SMARTS-parse tests (every rule's SMARTS parses; the
generic macrocycle/small-ring/standard rules were confirmed, in this same
review pass, to actually fire on the fixtures they were written for — see
`rules_macrocycle.rs`'s regression tests, added after independent review
found the macrocycle tier matched zero potentials on several corpus
fixtures) and translation provenance via direct source citation (§1 above
and the sources manifest), not a live RDKit differential.

Raw results: `validation/etkdg_torsion_knowledge_v2_chematic_side.json`
(chematic side) and `validation/etkdg_torsion_knowledge_v2_rdkit_oracle_diff.json`
(RDKit side + diff), both written by the reproduce commands above.

### 6.5 Second-round correction: two of the new generic rules cited a shadowed line

A second independent review pass on the two generic macrocycle rules added
in §6.4's fix (`ring_ch2_ch2_chain` and `ring_ch2_ether_chain` as they
existed at that point) found that citing "which RDKit line matches this
SMARTS pattern" is not the same claim as "which RDKit line actually governs
this bond in RDKit's real behavior" — because RDKit's matching loop is
first-match-wins (`TorsionPreferences.cpp`'s `doneBonds` bitset: a bond is
scored by whichever pattern matches it FIRST in file order; every later
pattern is skipped for that bond once one has matched — confirmed by
directly reading the matching loop, not assumed from the `.in` file's
comments alone).

Both original picks turned out to be shadowed by an earlier, broader pattern
in the same file that also matches the same fixture:

- The alkane-chain rule cited line 245 (`[CX4H2;r{9-}][CX4H2;r{9-}]`,
  requires H2 on both atoms, term `(3,+1,4.0)`). Line 244, immediately
  before it, is a strict superset (`[CX4&r][CX4&r]`, no H-count
  restriction, term `(3,+1,2.0)`) and therefore always wins for
  cyclododecane's backbone bonds. Checked lines 11-243 for any earlier
  pattern that could match first on a plain saturated hydrocarbon (no
  O/N/S/aromatic anywhere in the fixture): none exists. Renamed
  `macrocycle:ring_ch2_ch2_chain` → `macrocycle:ring_generic_cx4_chain`,
  re-cited to line 244, coefficient corrected 4.0 → 2.0.
- The ether-chain rule cited line 65 (`[C][CX4H2][OX2][!#1]`, atom4 = any
  non-H, term `(1,+1,2.0)`, single well). Line 60, earlier in the file,
  requires atom4 be aliphatic carbon *specifically* (`[C]`, a strict subset
  of `[!#1]`) with term `(3,+1,4.0)` (three wells) — and crown ether's
  `-O-CH2-CH2-O-` repeat unit genuinely has a plain aliphatic carbon on
  that position, so line 60 wins. Checked lines 11-59: line 27 (an
  O-C-O-C acetal pattern) requires the central carbon bonded to two
  *different* ring oxygens, which this repeat unit does not have (each
  carbon has exactly one O neighbor); lines 29-55 all require an aromatic
  atom; lines 56-59 require H-counts/atom-types this fixture's plain CH2
  doesn't have. `macrocycle:ring_ch2_ether_chain` kept its name, re-cited
  to line 60, coefficients corrected from a single 1-fold term to a 3-fold
  term (2.0 → 4.0, periodicity 1 → 3).

An attempted empirical cross-check (RDKit's own ETKDG conformer ensemble for
the crown-ether C-O bond, `useMacrocycleTorsions=True`) came back
inconclusive — a 60-bucket histogram of `[3,10,16,20,4,7]` out of 60
conformers is not a clean signature of either a single well or a symmetric
three-well distribution, most likely because the whole-molecule bounds
matrix couples this bond to its ring neighbors' own torsion terms, which a
single-bond histogram cannot isolate. The source-code control-flow argument
above is a structural fact about the cascade (not a statistical inference),
so it was treated as decisive over the noisier empirical proxy; both are
recorded here rather than only the one that "worked", per this program's
practice of reporting a negative/inconclusive check rather than omitting it.

Separately (not a citation error, an added rule): `macrocycle:lactam_n_calpha`
was added citing line 113 — the N-Calpha bond adjacent to the C(=O)-N bond
that has no real coverage for the NX3H1+CX4H2 combination (§ above). This
sharpens that gap claim: the unbranched `macrocyclic_amide` fixture's
C(=O)-N bond genuinely matches nothing, but its immediately adjacent N-CH2
bond does. The known-gaps claim in the PR body was corrected accordingly to
name the specific bond, not the whole fixture.
