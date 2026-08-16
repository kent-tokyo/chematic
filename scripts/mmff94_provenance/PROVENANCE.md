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
| Charges (partial bond charge increments) | `Code/ForceField/MMFF/Params.cpp` | `defaultMMFFPBCI`, `defaultMMFFChg` | `defaultMMFFPBCI` is already the cited source for chematic's existing `pbci_for` table (pre-dates this PR). **Phase 2 update (issue #227, 2026-08-16)**: `mmff94_charges_numeric`'s BCI bond-type source, flagged as a known follow-up in the Torsion entry below, investigated and fixed — see the dedicated Charges/BCI entry after the Torsion entry. |
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

**Production fix (issue #227, 2026-08-13, Stage C1: wildcard theta0 table
restoration).** Discovered while investigating why 7 of the 8 residual
Angle `empirical_rule` tuples were classified `zero_ka_table_row` rather
than `no_table_row`: RDKit's real `defaultMMFFAngleData` has 2,342 rows;
chematic's first port (Stage A/B state) had only 2,245 — missing exactly
the 97 rows where `type_i == type_k == 0`. These carry `ka == 0.0`
(RDKit's `isDoubleZero` placeholder) and exist purely to supply a
central-atom-type-only default `theta0`; RDKit's real
`getMMFFAngleBendParams` (`AtomTyper.cpp` lines ~3538-3554) reuses that
row's `theta0` verbatim when found (skipping the theta0 sub-rule
entirely) and only derives `ka` empirically — the mechanism 7/8 of the
residual tuples actually use. This PR restores only the missing table
DATA and the minimal guard needed to keep using it safe: `MMFF94_ANGLE_ENERGY`
2245→2342 rows (`mmff94_energy::tests::table_sizes` updated), and
`mmff94_angle_energy`'s `search` closure changed from `.map` to `.and_then`
so a `ka == 0.0` hit is never surfaced as a real parameter (a physically-
invalid zero-force-constant `Some` would otherwise regress the function).
This is **provably a no-op for every pre-existing input**: the old table
had zero `ka == 0.0` rows, so the new filter has nothing to filter among
them; only the 97 newly-visible wildcard rows are affected, and the
filter is exactly what stops them from being misused. Verified two ways:
(1) `cargo test --workspace` green, unchanged, before and after; (2) a
full 265-molecule corpus re-measurement via
`mmff94_strict_gate_remeasure_227` shows **zero status differences**
across all 265 molecules (`Ok`/`MissingParameters`/`MinimizationFailed`/
`UnsupportedAtomType` counts identical, per-molecule join confirmed, not
just aggregate counts) — this PR deliberately does not yet make any of
these 97 rows resolve to a real value; that is Stage C2's job. Split out
from a combined Stage C implementation specifically so this data-
completeness fix (small, mechanically verifiable, zero risk) could be
reviewed and land independently of C2's substantive new empirical-rule
logic.

**Production fix (issue #227, 2026-08-13, Stage C2: genuine empirical
Bond/Angle rule).** Stacked on Stage C1. Ports RDKit's real
`getMMFFBondStretchEmpiricalRuleParams`/`getMMFFAngleBendEmpiricalRuleParams`
(`Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp` at the pinned
commit) as `bond_empirical`/`mmff94_bond_energy_resolved`
(`crates/chematic-ff/src/mmff94_energy/bond.rs`) and
`angle_empirical_theta0`/`angle_empirical_ka`/`mmff94_angle_energy_resolved`
(`.../angle.rs`), strictly *after* the existing exact-table/eqLevel-ladder
chain in both cases — the empirical rule never overrides a real table
hit. New source tables, transcribed verbatim from the pinned `Params.cpp`
(cross-checked against `scripts/mmff94_angle_bond_gap_classify.py`'s own
independently-fetched copies): `MMFF94_COV_RAD_PAU_ELE` (18 rows),
`MMFF94_BNDK` (58 rows), `MMFF94_HERSCHBACH_LAURIE` (25 rows — the raw
table's third numeric column, `dp_ij`, is confirmed genuinely unused by
RDKit's own real `kb` formula by direct source read, so not stored), and
`MMFF94_ANGLE_Z`/`MMFF94_ANGLE_C` (11/7 rows). Atomic numbers and the
central atom's `crd`/`val`/`mltb`/`lin` come from the existing
`mmff94_numeric_type_registry.rs`, not re-transcribed.

**Resolution provenance.** Added `Mmff94Resolution`
(`.../mmff94_energy/mod.rs`): `DirectTable` / `EquivalentType { level }` /
`GenericAngleTypeFallback` (chematic's pre-existing, non-RDKit
`angle_type=0` net) / `EmpiricalBond` / `EmpiricalAngle`. The plain
`mmff94_bond_energy`/`mmff94_angle_energy` (unchanged since Stage A/B/C1)
keep their bare `Option<Params>` signature; `mmff94_bond_energy_resolved`/
`mmff94_angle_energy_resolved` are new, additive functions returning
`Option<(Params, Mmff94Resolution)>` — needed because eqLevel substitution
(Stage B) can land on a real row that is nonetheless the wrong parameter
for the original triple's chemistry (the #236 furan-collision failure
class), so `Some(...)` alone never proved correctness.

**The one unconfirmable tuple, left fail-closed:**
`(angle_type=0, type_i=43, type_j=18, type_k=63)`, `chembl_tier_b_0022`.
Numeric type 63 (C5A) has no row in the 55-entry eqLevel table (types
1-55 only). RDKit's real `MMFFDefCollection::operator()` returns
`nullptr` for such a type under both build variants, and
`MMFFAngleCollection::operator()`'s real eqLevel loop dereferences that
unchecked — undefined behavior in RDKit's own C++, confirmed by direct
source read. The live oracle nonetheless returns a finite, plausible
value for this triple (`ka=1.281584861919745, theta0=104.6`), confirmed
*stable* across 20 atom-index renumberings of the source molecule and
across a second, structurally unrelated molecule with the same local
(N43)-(S18)-(C63) pattern — a value deterministic in the atom types
alone, not per-call noise. This is consistent with either a real,
unidentified RDKit resolution path or RDKit deterministically reading
the same out-of-bounds heap memory adjacent to its own 55-row static
table on every call in the same process; the two are observationally
identical from the Python binding, and no source-level mechanism was
found that would make this a *defined* resolution. No direct/eqLevel
table row exists for this tuple at any ordering (confirmed by a fresh
scan of the parsed `defaultMMFFAngleData`). Per instruction, excluded
from this PR's empirical-parity claim and left fail-closed:
`mmff94_angle_energy_resolved` gates the empirical path on
`has_eq_level_row(type_i) && has_eq_level_row(type_k)`
(`TableSearch::UndefinedSubstitution`, distinct from a genuine
`TableSearch::Exhausted`) — narrowly scoped to this specific condition,
not disabling empirical resolution for any of the other 7 confirmed
tuples, nor for `type_i`/`type_k` ≤ 55, nor for a type > 55 that hits the
exact table directly. Regression-tested
(`angle_empirical_fails_closed_for_undefined_eq_level_substitution`).

**Ring-size-3-or-4 detection.** Ported `isAngleInRingOfSize3or4`
(`AtomTyper.cpp` lines ~357-398) verbatim as
`is_angle_in_ring_of_size_3_or_4` (`.../mmff94_minimizer.rs`, `pub` for
reuse from chematic-3d) — local bond adjacency, NOT SSSR, deliberately
distinct from the pre-existing SSSR-based `atoms_share_ring_of_size`
used for angle-*type* classification. Direct source read
(`AtomTyper.cpp` lines ~2746-2785) confirms the ring-size 60°/90°
override applies only in the from-scratch theta0 branch, never in the
`ka==0.0`-row theta0-reuse branch (Stage C1) — chematic's port preserves
this asymmetry exactly (`angle_empirical_theta0` only runs for
`TableSearch::Exhausted`, never `TableSearch::ZeroKa`).

**Production wiring.** `mmff94_bond_energy_resolved`/
`mmff94_angle_energy_resolved` now back chematic-ff's own
`bond_energy`/`angle_energy`/`stretch_bend_energy`
(`.../mmff94_minimizer.rs`) — the actual production energy/gradient
functions, not just a coverage check — and chematic-3d's independent
`compute_mmff94_coverage`, so the strict-gate's coverage decision and
the minimizer's actual energy computation stay consistent. The two
flanking bonds' `r0` are resolved once per angle triple and passed in by
the caller, matching RDKit's real `getMMFFAngleBendParams` (which
requires both flanking `getMMFFBondStretchParams` calls to succeed
before even attempting the empirical path). `mmff94_term_coverage_audit.rs`
gained `empirical_resolved` per-row fields and
`bonds_final_unresolved`/`angles_final_unresolved` molecule-level counts,
mirroring the pre-existing `stbn_final_unresolved`/`dfsb_resolved`
type-only-vs-final-resolution split (Priority 2B) — the type-only
`bonds_missing`/`angles_missing` axis is unaffected by this PR (80/191,
unchanged), while `bonds_final_unresolved`/`angles_final_unresolved`
drop to 0/25.

**End-to-end confirmation via the real policy path.** Two pre-existing
`chematic-3d` tests exercising `[C@H](F)(Cl)Br` under
`minimize_with_policy` — `mmff94_strict_refuses_chfclbr_with_missing_angle_params`
and `mmff94_with_uff_fallback_falls_back_and_reports_why_on_chfclbr` —
asserted the *old* "no MMFF94 table entry, must fall back to UFF"
behavior for this molecule's 3 halogen-C-halogen angles (exactly 3 of
the 7 oracle-confirmed empirical tuples). Both now fail their old
premise and are rewritten (`mmff94_strict_now_resolves_chfclbr_via_empirical_angle_rule`,
`mmff94_with_uff_fallback_no_longer_needs_to_fall_back_on_chfclbr`) to
assert the new, correct behavior: `Mmff94BondAngleStrict` now succeeds
directly, with a real nonzero empirical angle energy before minimization
and zero UFF fallback needed — direct, independent confirmation (a
second code path, not just the unit tests above) that the empirical
rule is genuinely wired into production, not just measured by a
custom-built harness.

**Corpus re-measurement (full 265-molecule Wave 1 corpus, production
`minimize_with_policy` via `mmff94_strict_gate_remeasure_227`).**
From the post-#314/#315-merge baseline (measured fresh on that exact
commit, not assumed from an earlier session's number): 265 total,
`Ok`=178, `MissingParameters`=83, `UnsupportedAtomType`=1,
`MinimizationFailed`=3 (87 failing). Stage C1 alone: **zero status
differences** across all 265 molecules (see C1's own entry above) — the
data restoration has no behavioral effect by itself. Stage C2 (this
PR): `Ok`=248, `MissingParameters`=13, `UnsupportedAtomType`=1,
`MinimizationFailed`=3 (**17 failing, 87→17**). Verified by a full
per-molecule join, not aggregate counts: **zero regressions** (no
previously-`Ok` molecule became non-`Ok`); the `MinimizationFailed` set
(`chembl_tier_b_0028`/`0029`/`0030`) and `UnsupportedAtomType` set
(`force_field_unsupported_probe`) are byte-identical before/after — not
newly caused, not masked. All 70 `Ok`-transitions are
`MissingParameters` → `Ok` (83→13 exactly accounts for the 70). Of the
13 molecules still `MissingParameters`, 12 have real, unrelated genuine
`Exhausted`-path Angle gaps (caffeine + 11 `chembl_tier_b` molecules,
25 angle instances total under the coverage audit's `angles_final_unresolved`)
and 1 is `chembl_tier_b_0022`, the fail-closed tuple above — none of
these were part of this PR's empirical-parity claim (the claim covers
exactly the 7 oracle-confirmed tuples plus the 1 Bond tuple). This
result was reproduced twice: once on a combined (pre-split) Stage C
implementation, and again independently after the C1/C2 split — the
split version's per-molecule results are byte-identical to the combined
version's (0 mismatches), confirming the split preserved exact
functional equivalence.

**Empirical-rule application counts** (issue #227 acceptance-gate
requirement: report term/molecule counts where empirical fired,
separately from the strict-gate pass-count): the empirical Bond rule
resolves 1 unique `(bond_type, ti, tj)` tuple; the empirical Angle rule
resolves 7 of 8 unique table_gap tuples (1 left fail-closed). Applied
across the 265-molecule corpus, this newly covers all Bond gaps
(`bonds_final_unresolved` 1→0 at the type-diagnostic level) and reduces
`angles_final_unresolved` from the type-only-diagnostic's 191 down to
25 — the remaining 25 are governed by different, unrelated causes (12
molecules, listed above) than what this PR's 7+1 confirmed tuples
address.

**Quality gates**: `cargo test --workspace` green (verified on both the
C1-only and C1+C2 states), `cargo clippy --workspace --all-targets -- -D
warnings` clean, `cargo fmt --all -- --check` clean.

**Production fix (issue #227 Phase 1, 2026-08-15): torsion parameter gap
root-caused to a bond-order-source bug, NOT a missing Halgren empirical
rule.** Fresh T0 audit on this exact commit (`mmff94_term_coverage_audit.rs`,
before any Phase 1 change): `torsions_missing` **257** instances across
**62/265** molecules (Tier A=1, Tier B=61) — every number below names its
producing tool and denominator per the acceptance gate's own requirement, so
none of these are interchangeable with the Bond/Angle counts elsewhere in
this file. Of the 257: **254** have `present_at_different_classification =
Some` (a real row exists in chematic's own, unmodified 926-row
`MMFF94_TORSION_ENERGY` table at a *different* classification code); **3**
are absent at every code 0..=8 (`table_gap`).

**One molecule excluded from every count above and below**:
`force_field_unsupported_probe` (SMILES `[P](C)(C)(C)=C`, `primary_category:
"force_field_unsupported"` in its own manifest entry — a deliberate
fail-closed probe fixture, confirmed by name/category, not investigated
further) throws a typing error (`atom 0 (Element(15)) was assigned MMFF94
numeric type 20 (CR4R), whose registry element is Element(6)`) before torsion
enumeration ever runs. All corpus-wide counts in this section are therefore
effectively over a 264/265 typing-succeeded population; this is stated
explicitly here rather than silently narrowing the denominator. Left
unfixed, per the Phase 1 directive's guidance not to touch atom typing
mid-PR for a molecule whose whole point is to be unsupported.

**Two hypotheses investigated and falsified, both against the live oracle,
before the real cause was found** (see
`crates/chematic-3d/examples/mmff94_torsion_empirical_diagnostic_227.rs` for
the tooling; `rdkit==2026.03.4`, this session's actual installed version —
not the 2026.03.3 pin earlier entries in this file used, see
`validation/results/pipeline_v2_vs_rdkit_environment_record.json` for the
prior, already-documented precedent that this exact version bump does not
confound status/coverage transitions):

1. *A from-scratch Halgren empirical torsion rule.* Structure derived from
   OpenBabel's `forcefieldmmff94.cpp` (GPL-2.0, a DIFFERENT project from
   RDKit — not copied from, hypothesized from its comments citing Halgren
   Part V pp. 631-632 Table X by page/table number, since the primary paper
   is paywalled and no other public source reproduces the formula; see the
   session's own research notes for the full source survey). Rule (b)
   ("aromatic b-c central bond", the branch that structurally matches
   254/254 of these instances, all of which have `bond_order_jk ==
   Aromatic`) predicts a UNIFORM `V1=0, V2=6.0, V3=0` for every case (β=6.0,
   π_bc=0.5, U_C=2.0 for all-carbon central atoms — the formula has no
   dependence on the terminal i/l atoms at all). **(oracle-validated, all
   254, not a sample): 0/254 match.** Real oracle values cluster into 7
   distinct `(V1,V2,V3)` tuples that vary with the terminal atoms — evidence
   this is not a central-bond-only empirical formula at work, ruling out
   rule (b) (and, by the same central-bond-only-input argument, every other
   Halgren empirical-rule branch (c)/(d)/(g)/(h) considered).
2. *RDKit's real eqLevel canonical-type-substitution ladder applied to the
   torsion's terminal atoms.* Ported the SAME `MMFF94_EQ_LEVEL` table
   Angle's Stage B already uses in production (`mmff94_energy/angle.rs`),
   substituting `ti`/`tl` through stages 3/4/5 with `tj`/`tk`/`tors_type`
   held fixed (mirroring Angle's real mechanism). **(oracle-validated
   indirectly via chematic's own 926-row table, all 254): 0/254 additional
   hits** beyond what the existing exact/wildcard chain already finds.

**The real cause**: for all 254, `mmff94_torsion_energy` at a DIFFERENT
classification code than the one `torsion_type_for` computed returns a value
that matches the live oracle exactly — **(oracle-validated, all 254): 254/254
row-level matches, zero exceptions**. `torsion_type_for`'s formula itself
(fixed 2026-08-10, see the Torsion table row above) is correct given its
input; the input was wrong. `assign_mmff94_numeric_types` already computes
an MMFF-specific re-perceived molecule
(`compute_mmff94_aromatic_view` — Kekulized, RDKit-`setMMFFAromaticity`-
matching bond orders, ported in Priority 1A) to derive atom TYPES correctly,
then discarded that molecule; every bond-order-dependent classification call
(`bond_type_for`/`angle_type_for`/`torsion_type_for`/`stretch_bend_type_for`)
kept reading `BondOrder` from the CALLER's original, un-reperceived molecule.
For caffeine's own pyrimidinedione ring (already the worked example in this
file's aromaticity-parity row above), chematic's general perception marks
bond 5-6 `BondOrder::Aromatic`; RDKit's real sanitizer Kekulizes the same
bond to `Single` — **(oracle-validated): confirmed 254/254 via
`MolFromSmiles(...).GetBondBetweenAtoms(j,k).GetIsAromatic() == False`**, and
independently corroborated by chematic's own pre-existing,
already-oracle-validated
`validation/results/mmff94_aromaticity_bond_parity_227_oracle.json` dump,
which already recorded `bond_aromatic["5-6"]: false` for caffeine before
this fix existed. `bond_type_for`'s own "`BondOrder::Aromatic` forces
`bond_type=0`" rule is itself correct (re-confirmed here, not just assumed
from the earlier benzene check) — it was simply being fed the wrong bond
order.

**Fix**: `assign_mmff94_numeric_types_with_view`
(`crates/chematic-ff/src/mmff94_numeric.rs`) returns `(types, mmff_mol)`
instead of discarding the re-perceived molecule; `assign_mmff94_numeric_types`
is now a thin wrapper. Threaded through chematic-ff's 5 production
energy/gradient entry points (`mmff94_total_energy`, `mmff94_torsion_scan`,
`mmff94_energy_breakdown`, `minimize_mmff94_full`, `minimize_mmff94_lbfgs`)
and chematic-3d's `compute_mmff94_coverage` (via `run_mmff94_bridge`, the
production `Mmff94BondAngleStrict`/`Mmff94WithUffFallback` gate) — a
classification/bond-order-source fix, not a new resolution tier. This is a
BROADER fix than Torsion alone: `bond_type_for`/`angle_type_for`/
`stretch_bend_type_for` share the same root cause, and all three improved as
a direct, unavoidable consequence (reported as its own line item, not folded
into the torsion numbers, same 265-molecule corpus, same audit tool):
`bonds_missing` (type-only) 80→1, `angles_missing` (type-only) 191→46,
`stbn_type_only_missing` (never gated by `Mmff94BondAngleStrict` at any
policy, reported for completeness only) 1865→1694.

**Known remaining inconsistency, explicitly not fixed here (follow-up, not
silent)**: `mmff94_charges_numeric` (`mmff94_numeric.rs`) independently reads
`bond.order` from the caller's original molecule for MMFF bond-charge-
increment (BCI) contributions, the same root-cause shape as the fix above,
unaddressed — electrostatics are out of scope for a torsion-parameter-gap
PR. Not measured to move any number in this PR (charges are not part of the
`Mmff94BondAngleStrict` coverage gate), flagged for whoever next touches
MMFF94 electrostatics.

**Renumbering-invariance, checked directly rather than assumed**:
`compute_mmff94_aromatic_view`'s Kekulization
(`chematic_core::kekulize`, the same blossom-matching solver canonical
SMILES already depends on, not a new per-bond heuristic written for this
fix) can in principle have genuine ties for a symmetric ring (e.g.
benzene's own two equally-valid alternating patterns) — a real risk that
the fix's re-perceived bond order could silently depend on atom traversal
order. Checked empirically, not assumed: `caffeine_reperceived_bond_order_is_invariant_under_atom_renumbering`
renumbers caffeine 32 ways (`deterministic_permutation`, the same
xorshift64-based generator issue #227 Priority 1A-2's own order-independence
tests already use) and re-identifies the C5A(63)-C=O(3) ring-fusion bond by
its unique atom-TYPE-pair content signature (not by index) in each variant
— **32/32 renumberings give the identical bond order**, no dependence
found. `benzene_reperceived_ring_bond_orders_are_invariant_under_atom_renumbering`
does the same for the textbook genuine-Kekule-tie case (two valid
alternating patterns) and confirms all 6 ring bonds resolve uniformly to
`Aromatic` (RDKit's real ring-level promotion for an accepted aromatic
ring, not a residual single/double pattern) across all 32 renumberings too
— whichever choice the Kekulizer's internal tie-break makes never leaks
into the final classification input. No renumbering-dependence bug found;
reported as a real, checked negative result, not assumed from the
pre-existing atom-TYPE permutation-invariance tests (which test a different
question — type assignment, not the bond ORDER this fix's classification
path additionally depends on).

**Second, contained addition**: `torsion_no_term_by_design`
(`mmff94_minimizer.rs`) and `Mmff94Resolution::NoTermByDesign`
(`mmff94_energy/mod.rs`). Halgren's real empirical-rule cascade omits the
torsion term entirely whenever either central (j-k) atom has MMFF's `lin`
flag (type 4 CSP, type 53 `=N=`, type 61) — rotating around a bond whose
other end is a linear 180° center changes no real geometry. This is the
complete, exact explanation for all 3 `table_gap` instances (2 in
`chembl_tier_b_0001`, an aryl-nitrile Ar-Ar-C#N shape with central k=type 4;
1 in `chembl_tier_b_0080`, an Ar-N=[N+]=[N-] cumulated-azide shape with
central k=type 53) — **(oracle-validated, all 3): `GetMMFFTorsionParams`
returns `None` for all 3, and each central atom's registered `linear` flag
matches RDKit's own MMFF atom type at that atom exactly** (types also
independently confirmed correct: chematic's assigned types match the
oracle's `GetMMFFAtomType` exactly on all 4 atoms of all 3 instances).
Wired into `compute_mmff94_coverage`'s torsion loop (a new
`torsions_no_term_by_design: usize` field on `Mmff94CoverageReport`,
included in `torsions_total`, excluded from `torsions_missing`) and into
`mmff94_term_coverage_audit.rs` identically, so these 3 no longer trip
`include_torsion_oop_in_gate` for their 2 molecules.

**Post-fix T-final audit** (`mmff94_term_coverage_audit.rs`, same 265-molecule
corpus, same denominators as T0 above): `torsions_missing` 257→**0**,
`torsions_no_term_by_design`=**3** (a denominator correction, not a
resolution — these never needed a parameter). `bonds_missing` (type-only)
80→**1**; `angles_missing` (type-only) 191→**46**; bond+angle
gate-would-fail 14→**13**.

**Corpus re-measurement, two named, distinct entry points, per-molecule join
(not aggregate arithmetic), zero regressions on both**:

- `minimize_with_policy(..., Mmff94BondAngleStrict, ...)` via
  `mmff94_strict_gate_remeasure_227.rs` (`generate_coords` starting geometry
  — the same tool/policy shape as `pipeline_v2_mmff94_strict`'s force-field
  arm, `include_torsion_oop_in_gate=false`, i.e. torsion coverage is NOT
  gated here by design): `Ok` 248→249/265 (**+1**, `MissingParameters`
  13→12, `UnsupportedAtomType`=1 unchanged, `MinimizationFailed`=3
  unchanged). Small delta expected: this gate never checked torsion
  coverage at all, so only the bond/angle side-effect of the fix (not the
  torsion fix itself) can move it.
- `minimize_with_policy_gated(..., Mmff94BondAngleStrict, ..., true, true)`
  via `mmff94_strict_gate_remeasure_227.rs --complete-bonded-term-gate`
  (matches the `chematic_pipeline_v2_mmff94_strict_complete_bonded_term_gated`
  arm's gate shape): `Ok` 187→**249**/265 (**+62**, `MissingParameters`
  74→**12**, `UnsupportedAtomType`=1 unchanged, `MinimizationFailed`=3
  unchanged). This is the gate the torsion fix (both the bond-order-source
  fix and the `NoTermByDesign` correction) actually targets. An intermediate
  measurement with the bond-order-source fix alone (before
  `torsion_no_term_by_design` was wired into `compute_mmff94_coverage`)
  showed `Ok`=247, `MissingParameters`=14 — the remaining 2
  (`chembl_tier_b_0001`, `chembl_tier_b_0080`) are exactly the 2 molecules
  `NoTermByDesign` fixes, confirmed by name via the same per-molecule join.

Both final re-measurements independently verified via a full per-molecule
join against their own pre-fix baseline on the same commit (not assumed from
the audit tool's aggregate numbers): **zero success→failure regressions on
either gate**.

**Recommendation and scope decision**: do NOT implement Halgren's empirical
torsion rule. On this corpus, as of this exact codebase state, zero real
instances need it — every apparent candidate was either a classification/
bond-order-source bug (254) or a genuinely term-free case RDKit also has no
parameter for (3). Shipping an empirical formula whose one falsifiable
prediction already failed against the live oracle would be strictly worse
than shipping none (per the Phase 1 directive's own "if you don't get an
improvement, report that fact plainly"). If a future corpus expansion
surfaces a genuine Torsion `table_gap` residual that is NOT explained by a
linear central atom, re-open this investigation with a fresh oracle
differential — do not assume this finding transfers without re-checking, the
same discipline this entry itself applied to the pre-2026-08-10 diagnostic's
now-superseded 254-instance empirical-rule hypothesis.

**Production fix (issue #227 Phase 2, 2026-08-16): BCI bond-type source —
oracle-confirmed real bug, fixed, NOT the same root cause as Phase 1's.**
Phase 1's Torsion entry above explicitly flagged `mmff94_charges_numeric`
(`crates/chematic-ff/src/mmff94_numeric.rs`) as reading `bond.order` from
the caller's original, un-reperceived molecule — "the same root-cause shape"
as the bug it had just fixed for bond/angle/torsion/stretch-bend
classification — and left it unaddressed as electrostatics-out-of-scope.
Phase 2's brief was to check this empirically before assuming it, per the
directive's explicit falsify-before-fix instruction. **The hypothesis as
literally stated was only half right**: the BCI code path had a compound
bug, not a single view-source bug.

*Source-level evidence first* (`Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp`,
pinned commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`): RDKit's real
`computeMMFFCharges` (lines 3071-3488) calls `unsigned int bondType =
this->getMMFFBondType(bond)` at line 3472 — textually the SAME
`getMMFFBondType` method `getMMFFBondStretchParams` calls at line 3500, on
the identical `bond` object from the identical (sanitized/Kekulized) `mol`
built once per `MMFFMolProperties` construction. There is no separate
"charge bond type" concept in RDKit's own algorithm; `getMMFFBondType`
itself (lines 2457-2475) is a *narrow* function returning 0 unless the bond
is formally `Bond::SINGLE` **and** both atom types are flagged `sbmb`/`arom`
(RDKit's "single bond between two conjugation-capable atoms" special case)
— it is never a function of bond multiplicity (Double/Triple/Aromatic all
collapse to 0 unconditionally, confirmed by a direct grep: `getMMFFChgParams`
is called from exactly this one call site in the whole pinned RDKit source,
always fed `getMMFFBondType`'s 0/1 result — no other bond-type value ever
reaches it).

*What chematic actually had*: a **second, independent, private**
`bond_type_for(order: BondOrder) -> u8` lived in `mmff94_numeric.rs`
(distinct from `crate::mmff94_minimizer::bond_type_for(ti, tj, order)`, the
function Phase 1 already fixed and oracle-validated for bond/angle/
torsion/stretch-bend), mapping `Single/Up/Down -> 0, Double -> 1, Triple ->
2, Aromatic -> 4` — a bond-*multiplicity* encoding with no resemblance to
RDKit's real `getMMFFBondType` formula at all, and no atom-type/sbmb
dependence whatsoever. This is why the population framing in the original
"same root-cause shape" hypothesis was incomplete: a plain "does the
Kekulized view differ from the original view" sample would systematically
undersample this bug, because it also fires on molecules where the two
views **agree** (e.g. any real, unambiguous C=O double bond, correctly
`BondOrder::Double` on both the original and reperceived molecule, still
got the wrong `bond_type` under the old formula).

*Table-level check* (zero-compile, pure text parsing, done before any oracle
call): chematic's `MMFF94_CHG` table itself is a faithful, byte-identical
port of RDKit's real `defaultMMFFChg` (both 498 rows, both with bond-type
column values `{0, 1, 4}` — confirmed by parsing a fresh
`/tmp/Params.cpp` fetch against the same pinned commit and diffing the
parsed row sets). The 3 `bond_type=4` rows (atom-type pairs 58/36, 58/37,
58/57) are demonstrably **unreachable** under RDKit's real algorithm too
(the single `getMMFFChgParams` call site never produces anything but 0 or
1) — vestigial data RDKit's own C++ never queries either, not evidence the
table needed a bond_type=4 concept. Of the table's 50 `bond_type=1` rows,
15 atom-type pairs also have a `bond_type=0` row for the same pair; 6 of
those 15 have a **materially different** `bci` value between the two rows
(up to 0.24 e-, e.g. `(4,9)`: −0.300 at bt=1 vs −0.106 at bt=0) — bounding
the maximum achievable per-contribution effect before running any molecule,
and confirming the Phase 1 directive's Step-1.4 off-ramp ("keyed by atom
TYPE pairs, not bond order, so the view doesn't matter") does NOT apply
here: bond type is a real, value-changing key in this table.

*Empirical falsification, full corpus (not a sample, same discipline as
Phase 1's 257/257 and 1,107/1,107 full-population checks)*: dumped
chematic's then-current (pre-fix) `mmff94_charges_numeric` output for every
heavy atom in all 264 typing-succeeded molecules
(`crates/chematic-3d/examples/mmff94_bci_charges_dump_227.rs`) and RDKit's
real `GetMMFFPartialCharge` for the same atoms, same molecules, no
embedding/conformer needed (topology-only, same precedent as the Torsion
investigation) via a live oracle
(`scripts/mmff94_bci_charges_oracle_227.py`, `rdkit==2026.03.4`). **Before
any fix**: 1,687/6,693 heavy atoms (25.2%) across 206/264 molecules (78.0%)
differed from the oracle by more than 1e-6 e⁻ — mean |Δ| 0.0189 e⁻, p90
0.076 e⁻, p99 0.239 e⁻, max 1.0 e⁻ (`chembl_tier_b_0080`, a separate,
unrelated cause — see below). This single number already falsifies any
"BCI already correct" or "bond-type doesn't matter here" hypothesis outright
— the pre-existing gap was large and corpus-wide, not a narrow edge case.

**Fix**: `mmff94_charges_numeric` now (1) calls
`assign_mmff94_numeric_types_with_view(mol)` and reads bond order from its
returned `mmff_mol` (the same reperceived view Phase 1 already threads
through the other four term kinds) instead of the caller's original `mol`,
closing the view-source half of the bug; and (2) the private, wrong
`bond_type_for(order)` is deleted entirely and replaced with a call to
`crate::mmff94_minimizer::bond_type_for(ti, tj, order)` — the SAME
already-fixed, already-oracle-validated function bond-stretch/angle/torsion/
stretch-bend classification uses, not a new third implementation — closing
the formula half of the bug. No new resolution tier, no table change: the
existing `lookup_chg_contribution`'s `unwrap_or_else(|| pbci_for(ti).0 -
pbci_for(tj).0)` fallback is itself RDKit's own real behavior
(`AtomTyper.cpp:3478-3479`, `mmffChgParams.second ? ... : ((*mmffPBCI)(atomType)->pbci
- (*mmffPBCI)(nbrAtomType)->pbci)`) and is unchanged.

**Post-fix re-measurement, same tool, same 264 molecules, same oracle
(literally the same dump script re-run, so before/after is a true diff, not
two different measurement methods)**: 67/6,693 atoms (1.0%) across 11/264
molecules (4.2%) still differ from the oracle — mean |Δ| 0.00116 e⁻, p90
0.0 e⁻ (exact match at the 90th percentile), p99 0.0144 e⁻, max 1.0 e⁻
(same outlier molecule, unchanged — see below). **Zero regressions,
verified by a genuine per-atom join (not aggregate-count arithmetic)**: 0
atoms that matched the oracle exactly before the fix now mismatch; 1,620
atoms moved from mismatched to exact-match; 5,006 were already exact-match
and remain so; 67 remain mismatched both before and after (the residual
below, unmoved by this fix either direction — direct evidence the fix
neither caused nor masked this separate gap).

**Regression-pinned unit test** (`acetone_carbonyl_charges_match_rdkit_oracle_after_bond_type_fix`,
`crates/chematic-ff/src/mmff94_numeric.rs`): acetone's C=O bond (type 3 —
type 7, `BondOrder::Double`) has a `MMFF94_CHG` row only at `bond_type=0`
(`(0, 3, 7, −0.57)`), no `bond_type=1` row for that pair. Under the old
formula, `Double -> 1` made this bond MISS the table and silently fall back
to the generic `pbci_for(3) − pbci_for(7)` PBCI difference; under the fixed
formula every Double/Triple/Aromatic bond maps to 0 unconditionally
(matching RDKit's real `getMMFFBondType`), landing on the real row directly.
Expected values (O = −0.57, carbonyl C = +0.448, both methyl C = +0.061)
are copied verbatim from a live RDKit oracle query, not derived from this
fix's own output.

**Residual, 67/6,693 atoms / 11/264 molecules, explicitly NOT investigated
further (out of scope for a BCI bond-type fix, flagged as follow-up)**: the
3 largest-magnitude residual molecules
(`chembl_tier_b_0080`/`_0159`/`_0161`) all show chematic's and RDKit's MMFF
atom TYPES agreeing exactly at every mismatched atom, with a charge
difference of almost exactly 0.5 e⁻ (one pair, the azide central/terminal N
of `chembl_tier_b_0080`, differs by 1.0 e⁻ — additive across the two
adjacent atoms of the same 0.5 e⁻-shaped effect) — a pattern consistent
with `mmff94_charges_numeric`'s formal-charge/`fcadj` redistribution step
(equation 15's `v·ΣformalCharge` term and/or the anionic-neighbor-leak step)
for charge-separated species (nitro N type 45, azide N types 47/53, charged
sulfoxide S type 17, O2CM type 32), NOT the bond-type BCI step this fix
addresses. Confirmed independent of this fix by construction (these 67
atoms are the "unchanged mismatch both before and after" set from the
per-atom join above — literally unmoved by the fix in either direction).
Not root-caused further here; flagged for whoever next touches MMFF94
formal-charge handling.

**Renumbering invariance, checked directly** (this fix depends on
`mmff_mol`'s Kekulized bond order, the same genuine-Kekule-tie-sensitive
input Phase 1's own reviewer follow-up required a test for):
`mmff94_charges_numeric_is_invariant_under_atom_renumbering`
(`crates/chematic-ff/src/mmff94_numeric.rs`) renumbers caffeine 32 ways
(`deterministic_permutation`/`rebuild_with_order`, the same helpers Phase
1's own renumbering tests use) and confirms the molecule's sorted per-atom
charge multiset — the only renumbering-stable representation of "which
charges this molecule produces", since individual atom indices are not
stable under relabeling — is identical (within 1e-9, float-summation-order
noise) across all 32 variants. No renumbering-dependence found.

**Quality gates**: `cargo test -p chematic-ff` 178 -> 182 passed (4 new
tests), 0 failed; `cargo clippy -p chematic-ff --all-targets --all-features
-- -D warnings` clean; `cargo fmt --all -- --check` clean.

**Follow-up investigation (issue #227 Phase 2, 2026-08-16): the BCI fix's one
new stereo violation, `chembl_tier_b_0082`.** The 3-state
`pipeline_v2_mmff94_strict` re-measurement (this file's Torsion/Charges
entries; full numbers in
`validation/results/mmff94_bci_gap_227_phase2_report.md`) found exactly one
molecule whose declared E/Z stereo flips from satisfied to violated between
State 2 (post-torsion-fix, pre-BCI-fix) and State 3 (post-BCI-fix): a
per-atom regression the earlier `per_molecule_join_regressions: 0` headline
number does NOT capture, since that field only tracks pipeline STATUS
transitions (success/typed_failure/timeout), not stereo-quality changes
within an otherwise-successful call — stated explicitly here and corrected
everywhere this file and the Phase 2 report use "0 regressions" language.

*Characterization* (`crates/chematic-3d/examples/mmff94_bci_stereo_drift_diagnostic_227.rs`,
using the new, purely additive
`chematic_3d::stereo_constraints::{debug_double_bond,debug_all_double_bonds}`
diagnostic — production `verify_double_bond`/`verify_stereo` untouched):
`chembl_tier_b_0082`'s single declared bond (`sub1=11 end1=13=end2=14
sub2=15`, declared trans/`same_side=false`) is Satisfied identically at
pre-minimization (`ForceFieldPolicy::None`, angle -140.7°) and State 2 final
(angle 166.9°, still trans-side) — both states share the same embedding
seed, and BCI charges do not affect embedding, so this identity is expected,
not a coincidence. State 3 final: angle 0.286° — `same_side=true`, Violated,
89.7° past the boundary (not a marginal case; a genuine ~167° rotation of
the alkene's own substituent dihedral occurred during minimization).

*Reproducibility*: `stereo_before` (post-embedding) is Satisfied in both
State 2's and State 3's own committed dump rows
(`mmff94_bci_gap_227_state{2,3}_*_chematic_rows.jsonl`, same
`embed_seed=20260801`/`max_attempts=8` `pipeline_v2_vs_rdkit_dump.rs` always
uses for `chematic_pipeline_v2_mmff94_strict` — confirmed this arm goes
through the seeded `EmbedParameters` path, not the legacy
`generate_coords_etkdg` entry point some other arms use) — the divergence
is real and reproducible at this fixed seed, not attributable to a
different embedding. A multi-seed sweep was not run (the single-seed
reproduction plus the RDKit oracle comparison below were judged sufficient
to determine the fix tier; not claimed as an exhaustive seed-independence
proof).

*RDKit oracle comparison (the key discriminator)*: chematic's own
`verify_stereo` judge, already applied to RDKit's saved geometries by
`pipeline_v2_vs_rdkit_common_scorer.rs`'s existing `score_rows` step (not a
new check), shows this bond **Satisfied on all 4 RDKit arms**
(`rdkit_etkdgv3_raw/uff/mmff94/best_of_n`) for `chembl_tier_b_0082` —
RDKit's own real MMFF94 minimizer, which always had correct BCI charges
(RDKit never had this bug), does not reproduce the flip. This rules out "the
corrected electrostatics legitimately crossed a real, RDKit-shared torsion
barrier" — it is a **chematic-specific gap**, not an expected physical
consequence of the charge fix.

*Energy-term isolation*: not run as a separate charge-swap experiment (the
oracle comparison above already answers the load-bearing question —
whether this is chematic-specific — more directly and with less risk of a
confounded methodology than re-scoring the same converged geometry with a
different charge set would have). `force_field_converged: false` at the
200-iteration cap in BOTH State 2 and State 3 (from the committed dump
rows) is the one other concrete data point: this molecule's minimization
does not converge either way, consistent with (not proof of) a
minimizer-robustness explanation rather than a clean two-basin energy
comparison.

*Architectural context, confirmed by direct source read, not assumed*:
`crates/chematic-ff/src/mmff94_minimizer.rs` (the MMFF94 minimizer itself)
has zero references to stereo anywhere in the file — the force field has no
notion of declared chirality/E-Z at all. `pipeline_v2.rs`'s own module docs
already document exactly this failure class for two earlier molecules
(`chembl_tier_b_0076`/`chembl_tier_b_0083`, found during the v0.14.0 release
gate: "a force field has no notion of declared chirality/E-Z and can walk a
geometry back across whichever stereo boundary a naive post-hoc repair
would have fixed") — `chembl_tier_b_0082` is a third instance of the SAME
already-known, already-documented architectural gap, not a new failure
class the BCI fix introduced from scratch; the BCI fix only changed which
molecule's already-marginal minimization trajectory happened to cross it.

*Fix tier chosen, in the directive's own priority order*: (1) a root-cause
fix (real stereo-awareness inside MMFF94 minimization) would be a large,
general architecture change — the project's own prior work on this exact
failure class explicitly deferred a structurally similar composition
question ("deliberately deferred rather than decided by omission") rather
than rushing one; out of scope here. (2) Converting this to a typed failure
under `StereoPolicy::Ignore` (the policy `chematic_pipeline_v2_mmff94_strict`
actually uses) would change default-arm behavior broadly, explicitly
flagged as needing its own separate authorization, not decided
unilaterally in this PR. (3) `StereoPolicy::RepairAndVerify` gets a new
post-minimization repair-and-reverify step (`crates/chematic-3d/src/pipeline_v2.rs`)
— empirically verified SAFE and EFFECTIVE for this exact case before
implementing: `repair_stereo` on State 3's violated geometry succeeds,
producing angle -179.7° (89.7° past the boundary, a robust result, not
marginal), with `worst_bond_length_ratio`/`gross_clash_count` both
IDENTICAL before and after the repair (0.1397/0, unchanged) — the
9-atom/13.8Å-displacement reflection is large in absolute terms but
introduces no bond-length or clash degradation, confirmed directly rather
than assumed. Implemented as an additive, fail-closed step: accepted only
if repair succeeds AND the reverified result has zero violations AND the
repaired geometry stays within `MAX_SANE_BOND_LENGTH` — any rejection falls
through to the original, unmodified `FinalStereoViolation` failure.
`StereoPolicy::Ignore`/`VerifyOnly` are completely unaffected by
construction (the new code path is nested inside the existing
`stereo_policy != Ignore` block and additionally gated on
`== RepairAndVerify`). (4) Not needed: tier 3 succeeded.

*What remains unrecovered*: `StereoPolicy::Ignore` — the policy this PR's
own 3-state measurement arm (`chematic_pipeline_v2_mmff94_strict`) actually
uses — never repairs, by design (its whole point is to never gate on
stereo), so `chembl_tier_b_0082`'s violation is real and NOT recovered
under the arm this Phase's headline numbers are measured on. This is a
known, named, reproducible, now-tested residual (see
`chembl_tier_b_0082_ez_bond_survives_bci_fix_under_repair_and_verify_not_under_ignore`,
`crates/chematic-3d/src/pipeline_v2.rs`), not silently absorbed into any
"0 regressions" claim.

**Test coverage gap, stated plainly rather than hidden**: no test in this
PR exercises the new post-minimization repair step's OWN fail-closed path
(repair genuinely fails, or succeeds but leaves a violation, or produces an
unsound geometry) with a real, from-scratch "unrepairable" molecule — a
grep of this file's existing test suite found zero pre-existing tests that
exercise `repair_stereo` returning `Err` at all (stage 8's own,
already-shipped repair-failure path has the same gap already, predating
this PR). The new code's fail-closed structure (three ANDed guard
conditions, falling through to the unchanged, pre-existing
`FinalStereoViolation` return otherwise) is directly auditable in the diff
and structurally identical in shape to stage 8's own `match repair_stereo
{...}` gate, but is not independently exercised by a dedicated
failing-repair integration test here. Flagged as a known gap, not silently
omitted.

**Quality gates (this addition)**: **under `--release`** (used for
faster iteration while developing this fix), `cargo test -p chematic-3d
--lib` showed 3 failures — 2 timing-race timeout tests
(`timeout_zero_fails_closed_with_typed_timeout`/
`timeout_failure_still_carries_evidence_computed_before_it_tripped`) and 1
already-known-jitter-molecule `atorvastatin_fragment` regression test in
`distance_geometry_v2.rs` — verified reproducing identically on the commit
immediately before this addition (none touch this PR's changed files, so
not caused by it). **Root-caused, not just deferred**: all 3 are
`--release`-build-specific — the timeout tests rely on a sub-millisecond
budget (`total_timeout_ms: Some(0)`) being exceeded by *any* nonzero
elapsed time, which release-mode optimization can defeat entirely for a
tiny molecule (pipeline completes in <1ms, `elapsed > 0` never trips); the
`atorvastatin_fragment` case is the same "borderline numerical case,
build-profile-sensitive" pattern the `chembl_tier_b_0082` investigation
above independently found for `stereo_before`'s violation count
(distance-geometry embedding's eigendecomposition + retry loop can pick a
different one of `max_attempts` attempts across optimization levels, even
at a fixed seed). **Confirmed via the actual gate command**: `cargo test
--workspace --all-features` (no `--release` — the command
`scripts/`/CI/the PR's own test plan actually runs) is fully green, 0
failures, workspace-wide, including all 3 of these — the debug-profile
build simply runs slower/differently enough that none of these
build-profile-sensitivity edge cases manifest. `cargo test -p chematic-3d
--lib` (same, no `--release`): 534 -> 537 passed (3 new/updated tests:
`repair_and_verify_recovers_post_minimization_stereo_violation`,
`post_minimization_stereo_repair_is_a_no_op_when_nothing_needs_recovering`,
`chembl_tier_b_0082_ez_bond_survives_bci_fix_under_repair_and_verify_not_under_ignore`),
0 failed. The `chembl_tier_b_0082` golden regression test's own assertions
were written to avoid pinning the exact violation count for this reason —
it pins the ROBUST property (`Ignore` never repairs;
`RepairAndVerify` always ends fully satisfied) rather than the
build-profile-sensitive exact stereo-before/after counts.
`cargo clippy -p chematic-3d --all-targets --all-features -- -D warnings`
clean; `cargo fmt --all -- --check` clean.

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
