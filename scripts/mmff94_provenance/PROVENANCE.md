# MMFF94 source provenance (Phase 1B-0, issue #227)

Pinned, not moving-`main`. All RDKit citations below resolve to one fixed
commit; re-fetch commands are given so this can be regenerated if the pin is
ever bumped, but no code in this repository reads RDKit's git history at
build or run time.

## Pinned RDKit revision

- **pip package installed for this audit**: `rdkit==2026.3.3` (`python -c
  "import rdkit; print(rdkit.__version__)"` → `2026.03.3`)
- **Matching RDKit git tag**: `Release_2026_03_3`
- **Commit SHA**: `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`
- **Resolved via**: `gh api repos/rdkit/rdkit/tags --paginate` (GitHub's tags
  API, not a local moving clone) — the tag→SHA mapping is immutable once
  published, so this is a fixed reference, not a moving target.
- **License**: BSD-3-Clause (RDKit `license.txt` at this commit) — reuse of
  parameter table data with attribution is permitted; this repo already has
  one precedent (`mmff94_numeric.rs`'s `defaultMMFFPBCI`-derived PBCI table,
  cited the same way before this PR).

## Source files and tables, by MMFF94 subsystem

All paths are relative to the RDKit repo root at commit
`e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`. Raw URL pattern:
`https://raw.githubusercontent.com/rdkit/rdkit/e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f/<path>`.

| Subsystem | File | Table / function | Notes |
|---|---|---|---|
| Atom typing rules (the actual classification logic — aromatic 5-ring/6-ring alpha/beta detection, N-oxide, imidazolium, etc.) | `Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp` | `MMFFMolProperties::setMMFFHeavyAtomType`, `RingMembershipSize`, `isAtomNOxide`, `isRingAromatic` | 3726 lines total; this repo ports only the aromatic C/N/O/S 5-ring/6-ring block (lines ~503–800 at the pinned commit), not the full aliphatic/metal/halogen switch — see the PR body for exact scope. |
| Numeric type definitions + equivalence-class levels | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFDef` string table (parsed by `MMFFDefCollection`) | Format: `SYMBOL\tTYPE\tEQ_LEVEL2\tEQ_LEVEL3\tEQ_LEVEL4\tEQ_LEVEL5\tDESCRIPTION...`; lines starting with `*` are secondary/alias symbols, skipped by RDKit's own parser (`inLine[0] != '*'`). Frozen copy: `scripts/mmff94_provenance/rdkit_defaultMMFFDef.txt` (extracted verbatim, C-string-escapes decoded, from the pinned commit — not hand-transcribed). |
| Atom-type intrinsic properties (atomic number, coordination, valence, lone-pair/multiple-bond/aromaticity/linear/sbmb flags) | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFProp` string table (parsed by `MMFFPropCollection`) | Format: `atype\taspec(atomic#)\tcrd\tval\tpilp\tmltb\tarom\tlin\tsbmb`. Frozen copy: `scripts/mmff94_provenance/rdkit_defaultMMFFProp.txt`. This is the field this PR's semantic-compatibility gate is built on — `aspec` is ground truth for "which element is this numeric type allowed to represent." |
| Bond types / bond-stretch parameters | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFBond` (params), `defaultMMFFBndk` (empirical-rule Badger's-rule-like constants), `defaultMMFFHerschbachLaurie`, `defaultMMFFCovRadPauEle` (covalent radii + Pauling electronegativities feeding the empirical bond rule) | Confirms MMFF94 *does* define a real empirical bond-stretch fallback (Halgren Part V) — chematic has none. Not implemented in this PR (see bond-fallback classification below); cited for the future-PR decision only. |
| Angle types / angle-bend parameters | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFAngleData` | Matches chematic's `MMFF94_ANGLE_ENERGY` table structurally (angle_type, i, j, k, ka, theta0). |
| Stretch-bend | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFStbn` (specific rows), `defaultMMFFDfsb` (default/generic stretch-bend constants by periodic-table-row pair — RDKit's own empirical fallback for stretch-bend) | chematic's `mmff94_stbn` has a partial fallback chain already (angle-type→0, then generic `(0,0,type_j,0)`); RDKit's `defaultMMFFDfsb` is periodic-row-keyed, not identical in shape. **Priority 2 update (issue #227, 2026-08-05)**: `defaultMMFFDfsb` (29 rows, frozen copy `scripts/mmff94_provenance/rdkit_defaultMMFFDfsb.txt`, extracted programmatically not hand-transcribed) is now ported as a **diagnostic-only** classifier in `mmff94_term_coverage_audit.rs` (`dfsb_default_resolvable` field) — NOT wired into chematic-ff's production energy/gate path. Algorithm verified against `Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp` at the pinned commit: `MMFFMolProperties::getMMFFStretchBendParams` (lines ~3566-3612) tries `MMFFStbnCollection::getMMFFStbnParams` (specific/generic by MMFF *type*, single exact lookup after I/K canonicalization — no equivalence-class step) first, and **only on failure** falls back to `MMFFDfsbCollection::getMMFFDfsbParams(periodicRow(atom1), periodicRow(atom2), periodicRow(atom3))`, keyed by periodic-table row via `getPeriodicTableRow` (`AtomTyper.cpp:251-264`: atomic number 1-2→row 0, 3-10→row 1, 11-18→row 2, 19-36→row 3, 37-54→row 4, else 0), canonicalized so `min(row1,row3) <= max(row1,row3)`. Confirmed **no equivalence-class (`eqLevel`) step exists anywhere in the stretch-bend resolution path** — `eqLevel` is used only by RDKit's angle/torsion/OOP fallback functions (`AtomTyper.cpp:527,552,743,768,862`), not stretch-bend — so a `equivalence_fallback_resolvable` bucket does not apply to stretch-bend under RDKit's real algorithm, and building one for it would misrepresent RDKit's own behavior, not just be incomplete. **Priority 2B update (issue #227, 2026-08-06)**: the same table is now ported into PRODUCTION `chematic_ff::mmff94_stbn` (`crates/chematic-ff/src/mmff94_energy/oop_stbn.rs`), unconditional, not gated — the diagnostic-only `dfsb_default_resolvable` field above was replaced by a `dfsb_resolved: bool` field (`mmff94_term_coverage_audit.rs`) that keeps the type-only diagnostic (`present_at_different_classification`, via the newly-split-out `mmff94_stbn_type_only`) and the final production-resolution question separate on purpose: of the 2,107 instances the type-only lookup misses, 1,680 are genuine table gaps that Dfsb closing matches RDKit's real behavior for, but 427 are routing-bug candidates (a real, correctly-typed parameter exists at a *different* classification code) that Dfsb *also* happens to resolve — coverage achieved, but chematic uses RDKit's generic default instead of the correctly-routed specific parameter for those 427, a parameter-selection-parity gap, not fixed by this port (see `validation/results/mmff94_coverage_227_term_audit_summary.json`'s `stretch_bend_dfsb_resolution` for the exact split). **Priority 2C update (issue #227, 2026-08-09)**: root-caused the 427, diagnostic-only, in `mmff94_stbn_equivalence_diagnostic_227.rs` — chematic's production stretch-bend code uses the **angle type** (0-8, `angle_type_for`) directly as the `MMFF94_STBN` key, but RDKit's real `getMMFFStretchBendParams` (`AtomTyper.cpp:3566-3612`) computes a *distinct*, finer-grained **stretch-bend type** (0-11) via `getMMFFStretchBendType(angleType, bondType1, bondType2)` (`AtomTyper.cpp:2480-2508`, ported verbatim) and `MMFFStbnCollection::getMMFFStbnParams`'s I/K canonicalization (`Params.h:601-663`: swap iff `iAtomType > kAtomType`, or tie-break on raw `bondType1 < bondType2` when equal) — not an `eqLevel` gap, confirming the Priority 2 finding above still holds. Also ported and cross-checked: `getMMFFAngleType`'s real ring-offset formula (`AtomTyper.cpp:2412-2447`: `angleType = size; if (bondTypeSum) angleType += bondTypeSum + size - 2`) disagrees with `angle_type_for`'s own ring-offset table for `bt_sum=2` (3-ring) and `bt_sum∈{1,2}` (4-ring) — a real, independent, second bug, measured LATENT on the 265-molecule corpus (0/113 reachable ring-embedded angle triples exercise the diverging branches); and `isAngleInRingOfSize3or4` (`AtomTyper.cpp:357-395`), which is local bond-adjacency, NOT SSSR-based (0/10,107 triples disagreed with chematic's SSSR-based ring check on this corpus, but the two are not the same algorithm). Live-oracle cross-check via `MMFFMolProperties.GetMMFFStretchBendParams`/`GetMMFFBondStretchParams` (`scripts/mmff94_stbn_oracle_validate_227.py`, same pinned RDKit build) confirms the ported formula matches RDKit exactly on 255/427 candidates (100%, zero exceptions) where chematic's own bond-order/aromaticity perception agrees with RDKit's; the remaining 172/427 are confounded by a separate, pre-existing aromaticity-perception gap (chematic trusts certain lowercase-aromatic input SMILES RDKit's sanitizer kekulizes instead — same mechanism also explains all 172 (molecule, triple) instances shared between `Angle`'s 277 and `StretchBend`'s 427 populations, set-identity verified). See the PR for issue #227 for the full breakdown. **Priority 2C production fix (issue #227, 2026-08-10)**: the classification-key bug diagnosed above is now fixed in PRODUCTION. `mmff94_stbn`/`mmff94_stbn_type_only` (`crates/chematic-ff/src/mmff94_energy/oop_stbn.rs`) take the `MMFF94_STBN` key as `stretch_bend_type: u8` (0-11), no longer `angle_type: u8` (0-8) — computed by a new `pub` `stretch_bend_type_for` (`crates/chematic-ff/src/mmff94_minimizer.rs`, ported verbatim from the diagnostic's `rdkit_stretch_bend_type`/`resolve_rdkit` arg-canonicalization logic). `angle_type_for`'s independently-diagnosed `bt_sum=2`-for-3-ring / `bt_sum∈{1,2}`-for-4-ring formula bug is also fixed in the same PR (root-cause dependency: a wrong `angle_type` feeds directly into `getMMFFStretchBendType`'s first argument), still measured LATENT on this 265-molecule corpus (0/113 reachable) — fixed for correctness, not because it moved any corpus number itself. Measured effect on the 427 `routing_bug_candidate` population, cross-validated two independent ways (a live RDKit oracle re-run via `scripts/mmff94_stbn_oracle_validate_227.py`, and a direct production-code cross-check against the frozen diagnostic's per-row predictions, zero mismatches either way): **220** moved from RDKit's generic Dfsb default to the correct, specific `MMFF94_STBN` parameter (parameter-selection parity achieved); **27** now correctly contribute ZERO stretch-bend energy (`mmff94_stbn_type_only` now finds the real `(0.0, 0.0)` row directly, `Some((0.0, 0.0))` — energy-equivalent to RDKit's own `isDoubleZero`-gated `None` drop, though not identical at the coverage-reporting layer: RDKit's API returns `None` here so `Mmff94CoverageReport` now counts these as resolved hits where RDKit would count them absent, a reporting nuance not an energy one; chematic previously injected a nonzero generic Dfsb value here, which was the real bug); **8** correctly remain on the generic Dfsb fallback (RDKit's own real algorithm also falls through here — chematic's output was already numerically right, now for the right reason too); **172** are the aromaticity-perception-confounded instances from the paragraph above, unchanged and still out of scope for this fix. All four counts verified by a direct per-row join (not just aggregate-count arithmetic): zero of the 172 confounded rows land in the 220 or 27 buckets. `mmff94_term_coverage_audit.rs`'s StretchBend `routing_bug_candidate` count itself moves 427 → 180 (the 220+27 that fully resolve no longer appear as "missing" at all; the 8+172 still do, now correctly classified as `stretch_bend_type` rather than the old wrong `angle_type`); `table_gap` (1,680) is untouched, as expected (this fix only corrects routing, not genuine data gaps). Coverage caveat: `stretch_bend_type_for`'s `ta == tc` branch is exercised 35/427 times in this corpus but never with `bond_type_ij != bond_type_jk` among those 35, so its asymmetric-bond-type corner is carried over unvalidated from the already-merged diagnostic, not independently re-verified here. See the PR for issue #227 for the full before/after corpus benchmark comparison. |
| Out-of-plane (Wilson angle) | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFOop` (MMFF94), `defaultMMFFsOop` (MMFF94s variant, not used by chematic) | |
| Torsion | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFTor` (MMFF94), `defaultMMFFsTor` (MMFF94s variant) | **Diagnostic update (issue #227, 2026-08-10, terminology corrected 2026-08-10 per reviewer request on PR #275)**: unlike stretch-bend (above), torsion's real resolution path (`Code/ForceField/MMFF/Params.h`'s `MMFFTorCollection::getMMFFTorParams`, lines 822-937 at the pinned commit) DOES run a genuine 4-stage `eqLevel` canonical-type-substitution ladder (the torsion/angle/OOP-only mechanism `PROVENANCE.md` already cited at `AtomTyper.cpp:743,768,862`), PLUS a second axis chematic has no equivalent of: `getMMFFTorsionType` (`AtomTyper.cpp:2528-2571`) computes a primary AND a "secondary" (pre-ring-override) torsion-type code, and the whole ladder is retried with the secondary code whenever the primary one is exhausted. Ported as a **diagnostic-only** tool, `crates/chematic-3d/examples/mmff94_torsion_equivalence_diagnostic_227.rs`, cross-checked against a live RDKit oracle for all 1,107 `routing_bug_candidate` Torsion instances (`validation/results/mmff94_coverage_227_term_audit_summary.json`). Below, every number is labeled **(oracle-validated)** -- a live `GetMMFFTorsionParams` call on that exact row, all 1,107 checked -- or **(self-port estimate)** -- derived only from comparing this diagnostic's own RDKit-formula port against chematic's code or against itself, no live oracle call for that specific number -- immediately after the number, not as a single end-of-paragraph disclaimer (an earlier version of this entry put one disclaimer at the paragraph's end, which read as if the whole paragraph, including the 76.3% figure, were oracle-measured; it is not).

Five-number breakdown, all 1,107 candidates, all **(oracle-validated)**: raw table/ladder row found = 853/1,107; valid non-zero table resolution = 851/1,107 (423 exact + 428 equivalence_level_3); explicit-zero row dropped = 2/1,107 (a real, found `MMFF94_TORSION_ENERGY` row with V1=V2=V3=0.0 that RDKit's own `isDoubleZero` gate drops to "no term", not a port bug); empirical-rule resolution = 254/1,107 (this diagnostic's table+ladder port finds nothing, but the live oracle still returns a nonzero term -- direct evidence of RDKit's separate Halgren empirical rule, `getMMFFTorsionEmpiricalRuleParams`, `AtomTyper.cpp:2874-3080`, chematic has no equivalent at all); final unresolved/no-term = 2/1,107 (the identical 2 rows as explicit-zero row dropped, by construction). Sum check: 851 + 254 + 2 = 1,107; 853 = 851 + 2 -- "853" and "851" are deliberately different numbers, not interchangeable synonyms for "resolved": 853 is every row where the table/ladder found *some* row (including the 2 RDKit itself does not count as a real term), 851 is the oracle-confirmed genuine non-zero resolution subset. Self-predicted `(V1,V2,V3)` values match the oracle exactly on 853/853 non-empirical rows (oracle-validated, 0 unexplained discrepancies).

Two findings, both checked directly rather than assumed from the stretch-bend PR's result or from each other: (1) `torsion_type_for`'s base-case classification FORMULA itself (`crates/chematic-ff/src/mmff94_minimizer.rs`) disagrees with RDKit's real one on 1,107/1,107 (100%) of the routing candidates (**self-port comparison** between chematic's code and this diagnostic's independently-ported `getMMFFTorsionType` -- RDKit's Python API does not expose the internal classification code directly, so this specific 100% figure is not itself an oracle call, though it is indirectly corroborated by the 853/853 oracle-validated value-match above) AND 10,325/13,530 (76.3%) of ALL torsion instances in the corpus (**self-port estimate, not independently oracle-validated at this scale** -- a corpus-wide sweep, not just the 1,107 candidates) -- chematic classifies purely from atom-type membership in `MLTB_TYPES` (`(MLTB(tj),MLTB(tk))` -> 0/1/2), never consulting the j-k bond's own order/MMFF bond type, whereas RDKit's real formula is `torsionType = bondTypeJK` (`getMMFFBondType`, 0 unless the bond is SINGLE and sbmb/aromatic-flagged on both ends) with a narrow empirical override to 2; of the 9,216 torsions (**self-port estimate, not oracle-validated**) that still resolve to SOME value today via chematic's crude fallback despite the classification mismatch, this diagnostic's own table+ladder port confirms 1,792 (**self-port estimate, not oracle-validated for this full sweep** -- only the 1,107-candidate population was oracle-checked) carry a numerically different `(V1,V2,V3)` than RDKit's real answer, a silent-wrong-parameter population an order of magnitude larger than the 1,107 instances `mmff94_term_coverage_audit.rs` can see (it only logs misses). (2) The eqLevel ladder itself **measures as contributing ZERO cases beyond chematic's EXISTING, UNMODIFIED `mmff94_torsion_energy` fallback chain** -- decisive check: `existing_fallback_resolves_with_corrected_code` = 853/1,107 and `existing_fallback_value_matches_ladder` = 853/853 (**self-port comparison between chematic's existing code and this diagnostic's ladder port -- not itself a live oracle call, but TRANSITIVELY oracle-confirmed**: the ladder port's own 853 values it matches are independently oracle-validated above, so chematic's existing fallback reproducing them exactly also confirms chematic's existing-fallback-plus-corrected-classification values against the oracle for those 853 rows), i.e. chematic's current lookup function, given ONLY the corrected classification code and no other change, already resolves every candidate this diagnostic's custom ladder resolves, to the identical value (all 428 "equivalence_level_3" hits are exactly chematic's own existing `(tors_type,0,tj,tk,0)` double-wildcard probe; ladder stages 1-2, the only genuinely NEW substitutions, never fire once). Revised recommendation (narrower than the ladder-focused hypothesis this diagnostic set out to test): production fix is `getMMFFTorsionType` ported into `torsion_type_for` -- a classification-only change, no lookup/fallback-chain modification needed -- which alone closes 851/1,107 (76.9%, oracle-validated); do NOT additionally build an eqLevel ladder, it is real in RDKit's source but measured as latent (0 incremental effect) on this corpus, same "real mechanism, unexercised" verdict the stretch-bend PR reached for a different sub-bug; Halgren's empirical rule is a separate, larger follow-up (the remaining 254/1,107 = 22.9%, oracle-validated, no existing chematic-ff equivalent). Not implemented by this diagnostic, which is read-only and does not modify `crates/chematic-ff/src` or `crates/chematic-3d/src`.

**Production fix (issue #227, 2026-08-10)**: the classification-formula bug diagnosed above is now fixed in PRODUCTION, exactly per the diagnostic's own narrower recommendation (no eqLevel ladder, no Halgren empirical rule). `torsion_type_for` (`crates/chematic-ff/src/mmff94_minimizer.rs`) now takes `mol: &Molecule` plus `ti`/`tl` (breaking signature change) and computes the base code as `bond_type_for(tj, tk, order_jk)` (reusing the already-correct, unmodified `bond_type_for`) with the empirically-required override to type 2 (`bond_type_jk==0 && order_jk==Single && (bond_type_ij==1 || bond_type_kl==1)`), instead of the old `(MLTB(tj),MLTB(tk))->0/1/2` atom-type-membership rule. The ring-4/5 override is also replaced end to end: a new private `ring_size_4_or_5` ports `isTorsionInRingOfSize4or5` (`AtomTyper.cpp:403-447`) faithfully -- local bond-adjacency, NOT SSSR-based -- and the 5-ring branch now additionally requires `ti==1 || tj==1 || tk==1 || tl==1`, a condition the old SSSR-based check had no equivalent of at all. Verified fresh against the actual production code (not restated from the diagnostic's self-port estimates): re-running `mmff94_term_coverage_audit` post-fix finds Torsion `routing_bug_candidate` 1,107 -> 254 (`table_gap` unchanged at 14, `torsions_missing` 1,121 -> 268); of the original 1,107, 853 now resolve to a raw table row via chematic's existing, unmodified `mmff94_torsion_energy` fallback chain alone (851 valid non-zero + 2 explicit-zero rows RDKit's own `isDoubleZero` gate would also drop, not counted as resolved), matching the diagnostic's own 851/1,107 (76.9%) prediction exactly; the remaining 254 need RDKit's out-of-scope Halgren empirical rule, also matching exactly. Beyond the 1,107 candidates, a corpus-wide before/after sweep of ALL 13,530 torsion instances (frozen old-formula copy vs. the new production code) found: 10,617 unchanged, 1,792 changed to a numerically different `(V1,V2,V3)` (matching the diagnostic's own self-port estimate of 1,792 exactly), 853 newly resolved (== the 853 above), and **0 newly lost** (no torsion that used to resolve to a wrong value now resolves to nothing). Oracle-validated via live RDKit (`rdkit==2026.3.3`): all 1,792 changed-value rows checked (not a sample) against `GetMMFFTorsionParams` -- **1,776/1,792 (99.1%)** match the new, post-fix value exactly, **0** match the OLD (pre-fix) value instead (i.e. zero cases where this fix made a previously-correct value wrong), and the remaining 16/1,792 (0.9%) are a pre-existing, out-of-scope MMFF94 aromaticity-perception gap for charged aromatic (pyridinium-type) rings unrelated to torsion classification itself (chematic types the ring carbons aromatic type 37, RDKit types them non-aromatic type 2 for these specific rings; both sides agree the torsion classification is type 0 either way). The 853 newly-resolved rows were already oracle-validated at 853/853 (100%) by the diagnostic itself. Also empirically characterized (not just reasoned about): chematic's torsion-enumeration loops have no `i==l` guard (a 3-membered-ring degenerate torsion where the two outer atoms coincide), which theoretically could make `ring_size_4_or_5`'s local-adjacency check misfire for a substituted 3-ring; measured 33 such instances in this corpus, 0 of which trigger the ring override either before or after this fix, and a constructed methylcyclopropane check against the live RDKit oracle confirms the port is faithful (RDKit's own real algorithm computes the identical ring-size-5 local adjacency for this case, but its own type-1 gate -- using RDKit's correct CR3R=22 ring-carbon typing, which chematic also assigns correctly here -- prevents the override from firing, matching chematic's result exactly). See the PR for issue #227 for the full before/after 265-molecule energy-pipeline benchmark comparison. |
| van der Waals | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFVdW` | |
| Charges (partial bond charge increments) | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFPBCI`, `defaultMMFFChg` | `defaultMMFFPBCI` is already the cited source for chematic's existing `pbci_for` table (pre-dates this PR). |
| Aromaticity perception feeding MMFF typing | `Code/GraphMol/Aromaticity.cpp` | `setMMFFAromaticity` (module-level function, not a `MolOps` member despite earlier notes in this file placing it there) | Priority 1A (issue #227): ported as `compute_mmff94_aromatic_view` in `mmff94_numeric.rs` — a **partial, behaviorally-calibrated** port, not a full one: every rule is a direct, line-cited port (ring-by-ring pi-electron counting at lines ~955-1035, the exocyclic-double-bond/NOS lone-pair-bonus rules, the multi-pass resolution loop) except the hybridization gate at line 1023 (`atom->getHybridization() != Atom::SP2`), approximated as `total_degree(atom) > 3` since chematic has no general hybridization-inference engine to port this faithfully. Measured gap on the 265-molecule Wave 1 corpus (`scripts/mmff94_hybridization_gate_gap_227_report.py`): 4,128/4,172 (98.9%) ring C/N atoms same decision as RDKit, 44 where the approximation under-triggers (misses a real pyramidal-SP3 ring N), 0 where it over-triggers, 0 unclassified — see `validation/results/mmff94_hybridization_gate_gap_227_report.txt`. RDKit's own general aromaticity model (distinct from both this MMFF-specific one and from chematic's `chematic_perception::apply_aromaticity`) is not relevant to MMFF typing and is out of scope here. |

**Production fix (issue #227, 2026-08-13): nitrile/sulfonamide/nitro/azide-N
and charged-sulfoxide-S typing.** Root-caused via a new diagnostic,
`scripts/mmff94_angle_bond_gap_classify.py`, that classifies each unique
`mmff94_strict`-blocking Angle/Bond `table_gap` tuple (issue #227's 2026-08-10
status comment: 97 Angle + 5 Bond instances, 27 + 3 unique atom-type tuples)
against RDKit's real resolution path (table lookup incl. the `eqLevel`
equivalence ladder, vs. the eq.18-20 empirical rule, vs. a chematic-side
atom-type mismatch) by parsing RDKit's own source tables (`defaultMMFFDef`,
`defaultMMFFAngleData`, `defaultMMFFBond`, `defaultMMFFBndk`,
`defaultMMFFHerschbachLaurie`, `defaultMMFFCovRadPauEle`) directly and
cross-checking against a live oracle. Result (instance-weighted): Angle 46%
(45/97) `type_mismatch`, 44% (43/97) genuinely need the empirical rule
(eq.18-20, tracked separately — see Torsion's own empirical-rule entry
above for the same mechanism class), 9% (9/97) resolvable via the `eqLevel`
ladder (also not yet ported); Bond 100% (5/5) `type_mismatch`. All 15 Angle
+ 3 Bond unique `type_mismatch` tuples traced to exactly the "5 small
pre-existing gaps" already named in this file's Priority-1A-2 history
(nitrile-N approximation, NSO2 sulfonamide/cyano-N, azide/diazo typing,
charge-shortcut masking nitro-N, charged-sulfoxide-S) — zero new typing
bugs found. Fixed in `assign_n_type`/`assign_s_type`
(`crates/chematic-ff/src/mmff94_numeric.rs`), each condition a direct,
line-cited port of RDKit's real `case 7`/`case 16` cascade
(`AtomTyper.cpp` lines ~971-1481 for N, ~1815-1917 for S at the pinned
commit): nitrile/isocyanide N (degree-1, triple-bonded) → type 42 (NSP);
terminal azide/diazo N → type 47 (NAZT); central charged cumulated
azide/diazo N → type 53 (`=N=`); nitro N → type 45 (NO2/NO3, was
incorrectly returning 46 "N=O" nitroso before this fix, unreachable in
practice because the pre-existing `charge > 0 -> 34` shortcut fired first
for every real charge-separated nitro group); sulfonamide/sulfonate N
(attached to a P/S bonded to ≥2 terminal oxygens) → type 43 (NSO2/NSO3,
reusing `classify_n_c3_carbon_context`'s pre-existing but previously-unwired
`is_cyano_like` field for the cyanamide half of the same RDKit flag);
charged sulfoxide S (e.g. `[S+]([O-])`, MMFF94's only valid charge-separated
spelling) → type 17 (S=O) — the old `assign_s_type` only counted *explicit
double bonds* to O, missing every charge-separated form entirely.
Re-measured on the full 265-molecule corpus: `mmff94_strict`'s own
bond+angle coverage gate (`Mmff94CoverageReport::bond_angle_fully_covered`)
moves from 107 to 103 failing molecules (4 flip to pass: `chembl_tier_b_0050`
sulfonamide, `_0080` azide, `_0159` nitro, `_0192` nitrile; **zero**
molecules regress from pass to fail), `bonds_missing` 84→80, `angles_missing`
374→358. Full `cargo test --workspace` (all crates) green throughout. Next:
the `eqLevel` equivalence ladder (9% of the Angle residual) and the
eq.18-20 Bond-stretch/Angle-bend empirical rule (44%) remain unimplemented,
tracked as this issue's next two stages.

**Production fix (issue #227, 2026-08-13, Stage B): RDKit's `eqLevel`
atom-type-equivalence ladder for Angle lookup.** Ported
`MMFFAngleCollection::operator()`'s real 4-stage canonical-type-
substitution ladder (`Code/ForceField/MMFF/Params.h` at the pinned commit)
as `eq_level`/`MMFF94_EQ_LEVEL` in
`crates/chematic-ff/src/mmff94_energy/angle.rs`: `angle_type` and the
central atom type stay fixed; `type_i`/`type_k` are substituted through
RDKit's real Level 2/3/4/5 equivalence classes (55-entry table, extracted
from `defaultMMFFDef`) before falling through to chematic's pre-existing
(non-RDKit, kept as a safety net) `angle_type=0` fallback. This is a
genuinely different, independent axis from the `type_mismatch` fixes in
Stage A: it also incidentally rescues many `routing_bug_candidate`
instances (which the "97 table_gap" framing didn't originally count), not
just the pure `table_gap` population. Effect was much larger than Stage A's
own diagnostic predicted (~9% of the *table_gap-only* residual) precisely
because of this: full 265-molecule corpus re-measurement moved
`mmff94_strict`'s bond+angle gate from **103 → 84 failing molecules** (19
more flip to pass, zero regressions), `angles_missing` 358→191. One
pre-existing test (`angle_type_for_butadiene_sp2_single_bond_...`) had a
stale assertion updated: it originally demonstrated "the old
hardcoded-angle_type=0 bug was a clean miss, not just less-specific" — that
premise no longer holds now that `eqLevel` makes angle_type=0 resolve too
(to a genuinely *different*, less-specific value, 118.043° vs the correct
121.55°) — updated to demonstrate the more precise point instead: getting
`angle_type` wrong now silently returns a *different wrong* parameter, the
same "silent wrong parameter" failure class as the #236 furan collision.
Re-running `scripts/mmff94_angle_bond_gap_classify.py` post-fix: the
Angle/Bond `table_gap` residual collapses to 8 + 1 unique tuples (from 27 +
3), **100% `empirical_rule`** (0 remaining `type_mismatch`/
`equivalence_table`) — Stage C (eq.18-20) is now the only remaining
mechanism gap, and 7/8 of its unique tuples are already oracle-confirmed by
this same diagnostic script.

## Halgren primary literature (secondary/theoretical cross-reference, not the implementation source)

- T. A. Halgren, "Merck Molecular Force Field. I. Basis, Form, Scope,
  Parameterization, and Performance of MMFF94," *J. Comput. Chem.* **17**,
  490–519 (1996).
- T. A. Halgren, "MMFF VI. MMFF94s option for energy minimization studies,"
  *J. Comput. Chem.* **20**, 720–729 (1999) — MMFF94s variant (RDKit's
  `defaultMMFFsOop`/`defaultMMFFsTor`), not used by chematic; noted for
  completeness only.
- T. A. Halgren, "MMFF VII. Characterization of MMFF94, MMFF94s...," *J.
  Comput. Chem.* **20**, 730–748 (1999).

This repo implements against RDKit's transcription of the Halgren tables
(RDKit's own `Params.cpp` comments cite "Copyright (c) Merck and Co., Inc.,
1994, 1995, 1996" directly on the table data), not against a fresh reading
of the papers — RDKit is the practical interoperability target (this whole
program's benchmark is "close the gap to RDKit"), and its numbering is
independently verified against the Halgren type names it carries in the
same table (`CB`, `C5A`, `NPYD`, etc., matching Halgren's own Table I
symbols).

## Regeneration

```bash
SHA=e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f
curl -sL "https://raw.githubusercontent.com/rdkit/rdkit/$SHA/Code/ForceField/MMFF/Params.cpp" -o /tmp/Params.cpp
# then re-run scripts/gen_mmff94_numeric_type_registry.py --extract /tmp/Params.cpp
# to refresh scripts/mmff94_provenance/rdkit_default{MMFFDef,MMFFProp}.txt
python3 scripts/gen_mmff94_numeric_type_registry.py
```

`rdkit_defaultMMFFDfsb.txt` was extracted from the same pinned `Params.cpp` (the
`defaultMMFFDfsb` string literal, lines ~4894-4932 at this commit) via a small
one-off Python script that regex-extracts every quoted string segment,
concatenates them, and unescapes `\t`/`\n` — not hand-transcribed. Re-run the
same extraction against a fresh `/tmp/Params.cpp` if the pin is ever bumped;
`getPeriodicTableRow`'s 5-way atomic-number bucketing (ported into
`mmff94_term_coverage_audit.rs` as `rdkit_periodic_table_row`) is from
`Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp` at the same pinned
commit, lines 251-264.
