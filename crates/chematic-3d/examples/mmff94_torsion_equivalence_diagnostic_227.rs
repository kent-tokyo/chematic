//! DIAGNOSTIC-ONLY (issue #227). Does NOT touch any production code path.
//!
//! Investigates the 1,107 Torsion `routing_bug_candidate` instances found by
//! `mmff94_term_coverage_audit.rs` (re-verified fresh before writing this
//! file: `torsions_missing=1121`, `routing_bug_candidate=1107`,
//! `table_gap=14`, unchanged from the frozen
//! `validation/results/mmff94_coverage_227_term_audit_summary.json`) -- a
//! real row exists in chematic-ff's `MMFF94_TORSION_ENERGY` table for the
//! exact atom-type quadruple, just not at the classification code
//! chematic-ff's production `torsion_type_for` computed for it.
//!
//! **Starting hypothesis, stated explicitly per the task**: a sibling
//! diagnostic (PR #273, `diag/mmff94-stretch-bend-routing-candidates-227`)
//! found stretch-bend's 427 routing candidates were NOT caused by a missing
//! `eqLevel` equivalence-fallback ladder -- RDKit's real stretch-bend
//! resolution path (`MMFFStbnCollection::getMMFFStbnParams`, `Params.h:601-`
//! `663`) has no `eqLevel` step at all, only one exact lookup then the
//! periodic-row Dfsb default. That PR's root cause was a distinct
//! classification-key bug (chematic used angle-type as the stretch-bend key
//! instead of computing RDKit's own finer-grained stretch-bend type).
//!
//! **Verified independently for torsion, from the same pinned RDKit source**
//! (`Code/ForceField/MMFF/Params.h`, `MMFFTorCollection::getMMFFTorParams`,
//! lines 822-937 at commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f` --
//! same pin as `scripts/mmff94_provenance/PROVENANCE.md`): torsion's real
//! resolution path is DIFFERENT in kind from stretch-bend's. It genuinely
//! does run a 4-stage `eqLevel` canonical-type-substitution ladder (`Params.h`
//! `AtomTyper.cpp:743,768,862` already cited in `PROVENANCE.md` as torsion/
//! angle/OOP-only, confirmed here by direct reading of the loop body), PLUS
//! a second, independent fallback axis chematic has no equivalent of at all:
//! RDKit's `getMMFFTorsionType` (`AtomTyper.cpp:2528-2571`) computes a
//! **pair** of torsion-type codes -- a primary code and a "secondary" code
//! (the classification the torsion WOULD have gotten before a ring override
//! bumped it to type 4 or 5, or literal `0` for non-ring torsions) -- and
//! retries the *entire* 4-stage ladder a second time with the secondary code
//! whenever the primary-code ladder is exhausted. chematic's
//! `torsion_type_for` returns a single `u8`, has no secondary-code concept,
//! and chematic's own `mmff94_torsion_energy` fallback chain (exact, reverse,
//! 2 single-end wildcards, double-wildcard, hardcoded-type-0 double-wildcard
//! -- 7 tiers total, see its own doc comment) is a much cruder approximation
//! of RDKit's real ladder: it wildcards straight to atom-type `0` rather than
//! substituting through the 4 real canonical equivalence levels, and its
//! final "type-0 generic" tier is a single double-wildcard probe, not a full
//! second 4-stage ladder. (Spoiler, quantified in Results below: this
//! structural gap is real but MEASURES as contributing zero cases beyond
//! chematic's existing crude fallback on this corpus -- don't stop reading
//! at "the ladder exists," the attribution is the whole point.)
//!
//! **A second, independent question this file also checks (learned from the
//! stretch-bend PR: don't assume the *classification formula itself* is
//! right just because a fallback ladder exists for it)**: does
//! `torsion_type_for`'s own primary-code FORMULA match RDKit's real
//! `getMMFFTorsionType`? It does not, structurally. chematic classifies the
//! non-ring base code purely from atom-type membership in the static
//! `MLTB_TYPES` set (`(MLTB(tj), MLTB(tk))` -> 0/1/2, ignoring the actual j-k
//! bond order entirely). RDKit classifies it from the **j-k bond's own MMFF
//! bond type** (`bondTypeJK = getMMFFBondType(bondJK)`, which is 0 unless the
//! bond is SINGLE and sbmb/aromatic-flagged on both ends -- a double/triple/
//! aromatic j-k bond always gets `bondTypeJK=0`, regardless of `tj`/`tk`'s
//! own MLTB membership), with an empirically-required override to type 2
//! (`bondTypeJK==0 && bond is SINGLE && (bondTypeIJ==1 || bondTypeKL==1)`).
//! These are genuinely different formulas that can and do disagree (see
//! `classification_mismatch_primary` in the stderr summary below) -- this
//! file reports both axes (ladder-gap vs. classification-formula bug)
//! separately, not conflated, exactly as the stretch-bend PR's review
//! insisted on for `angle_type_for` vs `getMMFFAngleType`.
//!
//! **Methodology note on `empirical_rule` vs `unresolved`**: RDKit's torsion
//! path has a THIRD fallback beyond the table+ladder -- Halgren's empirical
//! rule (`getMMFFTorsionEmpiricalRuleParams`, `AtomTyper.cpp:2874-3080`,
//! ~200 lines of per-element-pair V1/V2/V3 formulas keyed off atom
//! coordination/lone-pair/multiple-bond properties). This file does NOT
//! hand-port that formula a third time (chematic has no equivalent of it
//! either, a separate, larger gap out of scope for this diagnostic). Instead,
//! since every candidate is cross-checked against a live RDKit oracle anyway
//! (see below), the empirical-rule/unresolved split is determined
//! *empirically from the oracle's own returned value*, not self-computed:
//! when this file's OWN from-scratch table+ladder port finds nothing
//! (`selected_parameter_kind = "table_unresolved"` in the raw JSONL below)
//! but the live oracle still returns a nonzero torsion term, that is direct
//! evidence RDKit's empirical rule produced it; if the oracle ALSO returns
//! `None`, that is a genuine residual gap even under RDKit's complete
//! algorithm. `scripts/mmff94_torsion_oracle_validate_227.py` performs this
//! relabeling and writes the fully-enriched final JSONL (with
//! `used_exact`/`used_equivalence`/`used_empirical`/`used_unresolved`) as a
//! separate file, `..._oracle_enriched.jsonl` -- kept distinct from this
//! file's own raw, oracle-independent prediction so the two are never
//! silently merged.
//!
//! Everything RDKit-derived below is a *fresh, independent* port, not a call
//! into chematic-ff's own fallback logic: `MMFF94_TORSION_ENERGY` (chematic's
//! already-ported copy of RDKit's `defaultMMFFTor` table, exported `pub` from
//! `chematic_ff::mmff94_energy`) is reused as raw DATA only, probed here with
//! an from-scratch exact-lookup helper (`raw_lookup_torsion`) that does none
//! of chematic's own wildcarding -- so this diagnostic cannot be "correct by
//! construction" from sharing the bug it measures. `bond_type_for` is reused
//! as-is (same precedent as the stretch-bend PR: it is a direct, previously-
//! verified port of RDKit's `getMMFFBondType`, not part of the classification
//! logic under test). The `eqLevel` canonical-type table is parsed at run
//! time from the ALREADY-frozen, provenance-cited
//! `scripts/mmff94_provenance/rdkit_defaultMMFFDef.txt` (extracted verbatim
//! from the pinned commit; format documented in `PROVENANCE.md`), not
//! hand-transcribed here.
//!
//! **Results (live RDKit oracle `rdkit==2026.3.3`, matching the pinned commit's
//! release tag, ALL 1,107 candidates checked, not a sample)**:
//! - `classification_mismatch_primary = 1107/1107 (100%)`: not a single one
//!   of the 1,107 candidates has chematic's `torsion_type_for` output
//!   matching RDKit's real primary classification code. This is the
//!   dominant effect, confirmed structurally above, not a coincidence of
//!   this particular slice: a corpus-wide sweep over ALL 13,530 torsion
//!   instances (not just the 1,107 "missing" ones) finds
//!   `all_torsions_mismatch = 10,325/13,530 (76.3%)` -- the classification
//!   formula disagrees on the large majority of ALL torsions in this corpus,
//!   most of which currently still resolve to SOME value via chematic's own
//!   crude wildcard-to-0 fallback (`all_torsions_hit_with_mismatched_code =
//!   9,216`), of which this file's own (non-oracle-validated for this
//!   full-population sweep -- only the 1,107-candidate population below was
//!   oracle-checked) table+ladder port CONFIRMS `1,792` also carry a
//!   numerically different `(V1,V2,V3)` than what RDKit's real
//!   classification+ladder selects (0 "undetermined" -- see
//!   `all_torsions_hit_with_mismatched_code_undetermined` in the stderr
//!   summary) -- a silent wrong-parameter population an order of magnitude
//!   larger than the 1,107 "missing" instances this diagnostic was scoped
//!   to, invisible to `mmff94_term_coverage_audit.rs` entirely (it only logs
//!   misses).
//! - `ladder_resolves_same_code = 0`: feeding chematic's OWN (wrong)
//!   classification code through a fully-correct eqLevel ladder resolves
//!   **none** of the 1,107 candidates -- confirming a ladder alone, bolted
//!   onto the existing (buggy) classification, would fix nothing.
//! - Using RDKit's REAL classification (primary + secondary) through the
//!   same ladder: **851/1,107 (76.9%) resolve via the table + eqLevel
//!   ladder** (423 at the exact/own-type level, 428 at ladder stage 3),
//!   **254/1,107 (22.9%)** require RDKit's separate Halgren empirical-rule
//!   fallback (which chematic has no equivalent of at all -- confirmed, not
//!   guessed, by checking the oracle still returns a nonzero term where this
//!   file's own table+ladder port finds nothing), and **2/1,107 (0.2%)**
//!   land on a real, explicit all-zero `MMFF94_TORSION_ENERGY` row that
//!   RDKit's own `isDoubleZero` gate drops to "no term" (matching the
//!   sibling PR's `found_but_zero_dropped` pattern) -- genuinely zero
//!   contribution under RDKit's real algorithm either way, not a residual gap.
//! - **The eqLevel ladder measures as CONTRIBUTING NOTHING INCREMENTAL on
//!   this corpus** (checked directly, not assumed from the 76.9% figure
//!   above, which could equally be explained by a real ladder effect): ALL
//!   428 stage-3 hits resolve to key `(tors_type, 0, tj, tk, 0)` --
//!   `EQ_LEVEL5` is 0 for essentially every organic atom type (verified: 0
//!   rows in `rdkit_defaultMMFFDef.txt` have `EQ_LEVEL2 != TYPE`, so stage 0
//!   is always a literal exact match; separately, every stage-3 hit in this
//!   population happens to be a full double-wildcard), which is EXACTLY
//!   chematic's own EXISTING tier-4 `search(tors_type, 0, tj, tk, 0)` probe.
//!   Stages 1 and 2 (genuine equivalence-class substitutions, not
//!   wildcarding) never fire once in this population (`0` hits at either
//!   level). Decisive check: `existing_fallback_resolves_with_corrected_code
//!   = 853/1107` and `existing_fallback_value_matches_ladder = 853` --
//!   chematic's `mmff94_torsion_energy`, completely UNMODIFIED, fed ONLY the
//!   corrected classification code (no eqLevel ladder port involved at all),
//!   already resolves every single one of the 853 candidates this file's
//!   custom ladder resolves, to the IDENTICAL value. The eqLevel ladder is
//!   real in RDKit's source (as `PROVENANCE.md` already documented) but
//!   **measures as latent on this corpus** -- the same "real mechanism, 0
//!   measured incremental effect" verdict the sibling stretch-bend PR
//!   reached for `angle_type_for`'s ring-offset bug (0/113 reachable), just
//!   for a different mechanism.
//! - **853/853 (100%) of the self-predicted `(V1,V2,V3)` values for
//!   non-empirical rows match the oracle exactly** (0 unexplained
//!   discrepancies) once the 2 zero-dropped rows are accounted for --
//!   validating both this file's `getMMFFTorsionType` port and its eqLevel
//!   ladder port bit-for-bit against the real library, not just directionally.
//!
//! **Bottom line**: the task's starting hypothesis ("a missing eqLevel
//! ladder is torsion's real mechanism, unlike stretch-bend") is only
//! half right, in a way that matters for production scope. The ladder
//! mechanism DOES genuinely exist in RDKit's source (unlike stretch-bend) --
//! but it measures as contributing ZERO cases beyond what chematic's
//! existing wildcard-to-0 fallback already reaches, once fed a correct
//! classification code. The classification-formula bug (chematic's
//! atom-type-membership-based `(MLTB(tj),MLTB(tk))` rule vs. RDKit's real
//! j-k-bond-type-based `bondTypeJK` rule) is the ENTIRE fixable story for
//! 851/1,107 (76.9%) of this population, and a distinct, much LARGER,
//! corpus-wide bug on its own (76.3% of ALL 13,530 torsions, not just these
//! 1,107). Recommendation for the next (production) step, narrower than
//! originally hypothesized: (1) port `getMMFFTorsionType` faithfully into
//! `torsion_type_for` -- a classification-only fix, no lookup/fallback-chain
//! changes needed at all, closes 851/1,107 (76.9%) of this diagnostic's
//! population using chematic's EXISTING `mmff94_torsion_energy` unmodified,
//! plus an unmeasured (self-port-estimated at up to 1,792, not yet
//! oracle-validated) share of the larger 9,216-instance silent-wrong-value
//! population found by the full-corpus sweep; (2) do NOT additionally build
//! an eqLevel ladder as part of this fix -- it is real in RDKit but measured
//! as contributing nothing incremental here; revisit only if a future,
//! larger, or differently-shaped corpus exercises stages 1/2 (never
//! observed on this one); (3) separately scope Halgren's empirical rule as
//! its own, larger follow-up (closes the remaining 254/1,107 = 22.9%, a
//! real, non-trivial mechanism with no existing chematic-ff equivalent at
//! all -- do not fold it into (1)'s estimate).
//!
//! Run:
//! ```text
//! cargo run --release -p chematic-3d --example mmff94_torsion_equivalence_diagnostic_227 \
//!   > validation/results/mmff94_torsion_equivalence_diagnostic_227.jsonl \
//!   2> validation/results/mmff94_torsion_equivalence_diagnostic_227_stderr.log
//! ```
//! Oracle validation (writes the enriched final JSONL + prints the headline
//! resolution counts):
//! ```text
//! .venv/bin/python scripts/mmff94_torsion_oracle_validate_227.py \
//!   validation/results/mmff94_torsion_equivalence_diagnostic_227.jsonl
//! ```

use std::collections::BTreeMap;

use chematic_core::{AtomIdx, Molecule};
use chematic_ff::mmff94_energy::MMFF94_TORSION_ENERGY;
use chematic_ff::{assign_mmff94_numeric_types, bond_type_for, mmff94_torsion_energy};
use chematic_perception::find_sssr;
use serde_json::{Value, json};

// ── Corpus loading (same manifests as mmff94_term_coverage_audit.rs) ────────

struct CorpusMol {
    tier: String,
    name: String,
    smiles: String,
}

fn load_corpus() -> Vec<CorpusMol> {
    let mut out = Vec::new();
    for (tier, path) in [
        (
            "A",
            "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
        ),
        (
            "B",
            "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
        ),
    ] {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
        let manifest: Value =
            serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {path}: {e}"));
        for m in manifest["molecules"].as_array().expect("molecules array") {
            out.push(CorpusMol {
                tier: tier.to_string(),
                name: m["name"].as_str().unwrap().to_string(),
                smiles: m["smiles"].as_str().unwrap().to_string(),
            });
        }
    }
    out
}

// ── eqLevel canonical-type table (parsed from the frozen provenance copy) ──

/// `type -> [eq_level2, eq_level3, eq_level4, eq_level5]`, parsed from
/// `scripts/mmff94_provenance/rdkit_defaultMMFFDef.txt` (RDKit's
/// `defaultMMFFDef`, format per `PROVENANCE.md`:
/// `SYMBOL\tTYPE\tEQ_LEVEL2\tEQ_LEVEL3\tEQ_LEVEL4\tEQ_LEVEL5\t...`). Lines
/// whose first field is `*` are secondary/alias symbol rows, skipped exactly
/// as RDKit's own parser skips them (`inLine[0] != '*'`,
/// `MMFFDefCollection`'s constructor).
fn load_eq_level_table(path: &str) -> BTreeMap<u8, [u8; 4]> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    let mut table = BTreeMap::new();
    for line in text.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() || fields[0] == "*" || fields.len() < 6 {
            continue;
        }
        let Ok(atype) = fields[1].parse::<u8>() else {
            continue;
        };
        let Ok(e2) = fields[2].parse::<u8>() else {
            continue;
        };
        let Ok(e3) = fields[3].parse::<u8>() else {
            continue;
        };
        let Ok(e4) = fields[4].parse::<u8>() else {
            continue;
        };
        let Ok(e5) = fields[5].parse::<u8>() else {
            continue;
        };
        table.insert(atype, [e2, e3, e4, e5]);
    }
    table
}

/// `eqLevel[idx]` for a given atom type, falling back to the type itself if
/// absent from the table (should not happen for any real MMFF94 numeric
/// type; counted separately if it does, see `eq_level_fallback_count`).
fn eq_level(table: &BTreeMap<u8, [u8; 4]>, t: u8, idx: usize, fallback_count: &mut usize) -> u8 {
    match table.get(&t) {
        Some(levels) => levels[idx],
        None => {
            *fallback_count += 1;
            t
        }
    }
}

// ── RDKit-faithful ports (independent of chematic-ff's production code) ────

/// Port of RDKit's `isTorsionInRingOfSize4or5` (`AtomTyper.cpp:403-447`).
/// Purely local bond-adjacency, NOT SSSR-based: 4-ring iff i-l are directly
/// bonded; 5-ring iff i and l (excluding their ring neighbours j and k
/// respectively) share a common neighbour. Returns 0 if neither.
fn rdkit_ring_size_4_or_5(mol: &Molecule, i: AtomIdx, j: AtomIdx, k: AtomIdx, l: AtomIdx) -> u8 {
    if mol.bond_between(i, l).is_some() {
        return 4;
    }
    let nbrs_i: std::collections::BTreeSet<AtomIdx> = mol
        .neighbors(i)
        .map(|(n, _)| n)
        .filter(|&n| n != j)
        .collect();
    let has_common = mol
        .neighbors(l)
        .map(|(n, _)| n)
        .filter(|&n| n != k)
        .any(|n| nbrs_i.contains(&n));
    if has_common { 5 } else { 0 }
}

/// `(primary_torsion_type, secondary_torsion_type)`. Port of RDKit's
/// `getMMFFTorsionType` (`AtomTyper.cpp:2528-2571`). `secondary` is 0 unless
/// a ring override (type 4 or 5) fired, in which case it is the type the
/// torsion would have gotten WITHOUT the ring override -- RDKit's real
/// last-resort fallback code, not reproducible from chematic's
/// `torsion_type_for` (which has no secondary-code concept at all).
#[allow(clippy::too_many_arguments)]
fn rdkit_torsion_type(
    mol: &Molecule,
    i: AtomIdx,
    j: AtomIdx,
    k: AtomIdx,
    l: AtomIdx,
    ti: u8,
    tj: u8,
    tk: u8,
    tl: u8,
) -> (u8, u8) {
    let order_ij = mol.bond_between(i, j).expect("i-j bond").1.order;
    let order_jk = mol.bond_between(j, k).expect("j-k bond").1.order;
    let order_kl = mol.bond_between(k, l).expect("k-l bond").1.order;
    let bond_type_ij = bond_type_for(ti, tj, order_ij);
    let bond_type_jk = bond_type_for(tj, tk, order_jk);
    let bond_type_kl = bond_type_for(tk, tl, order_kl);

    let mut torsion_type = bond_type_jk;
    // MMFF.IV page 609's simple condition fails CYGUAN01 in RDKit's own test
    // suite; this empirically-corrected condition is RDKit's real code.
    if bond_type_jk == 0
        && order_jk == chematic_core::BondOrder::Single
        && (bond_type_ij == 1 || bond_type_kl == 1)
    {
        torsion_type = 2;
    }

    let mut secondary = 0u8;
    let ring_size = rdkit_ring_size_4_or_5(mol, i, j, k, l);
    if ring_size == 4 && mol.bond_between(i, k).is_none() && mol.bond_between(j, l).is_none() {
        secondary = torsion_type;
        torsion_type = 4;
    } else if ring_size == 5 && (ti == 1 || tj == 1 || tk == 1 || tl == 1) {
        secondary = torsion_type;
        torsion_type = 5;
    }

    (torsion_type, secondary)
}

/// Exact-or-full-reverse lookup into chematic's OWN raw `MMFF94_TORSION_ENERGY`
/// table, with NO wildcarding of any kind -- deliberately more restrictive
/// than chematic's own `mmff94_torsion_energy`, so this file's own
/// eqLevel-driven substitution (below) is the only source of wildcarding in
/// this diagnostic's ladder. "Full reverse" (all four positions, not just
/// the outer two) matches chematic's own `mmff94_torsion_energy`'s reverse
/// tier (`(ri,rj,rk,rl) = (tl,tk,tj,ti)`), which is itself how a
/// single-direction-canonicalized table (RDKit's own storage convention) is
/// probed from an uncanonicalized query -- trying both directions against a
/// canonical-single-direction table is equivalent to canonicalizing the
/// query first, without needing to replicate RDKit's exact swap-rule code.
fn raw_lookup_torsion(tt: u8, a: u8, b: u8, c: u8, d: u8) -> Option<(f64, f64, f64)> {
    let find = |t0: u8, t1: u8, t2: u8, t3: u8, t4: u8| {
        MMFF94_TORSION_ENERGY
            .binary_search_by_key(&(t0, t1, t2, t3, t4), |&(u0, u1, u2, u3, u4, ..)| {
                (u0, u1, u2, u3, u4)
            })
            .ok()
            .map(|idx| {
                let (.., v1, v2, v3) = MMFF94_TORSION_ENERGY[idx];
                (v1, v2, v3)
            })
    };
    find(tt, a, b, c, d).or_else(|| find(tt, d, c, b, a))
}

/// One eqLevel-ladder stage's `(iWildCard, lWildCard)` index pair into the
/// `[eq_level2, eq_level3, eq_level4, eq_level5]` array, per RDKit's real
/// loop (`Params.h:853-861`): stage 0 = own type (level 2, `EQ_LEVEL2` is
/// always identical to the atom's own type in the source table, confirmed by
/// inspecting every row of `rdkit_defaultMMFFDef.txt` -- matching the
/// angle-ladder code's own comment "we skip 1-1-1 since Level 2 === Level
/// 1"), stage 1 = (level3, level5), stage 2 = (level5, level3), stage 3 =
/// (level5, level5). `j`/`k` are NEVER substituted at any stage.
const LADDER_STAGES: [(usize, usize); 4] = [(0, 0), (1, 3), (3, 1), (3, 3)];

struct LadderHit {
    stage: u8,
    key: (u8, u8, u8, u8, u8),
    value: (f64, f64, f64),
}

/// Runs the 4-stage eqLevel ladder for one torsion-type code, first hit wins
/// (matches RDKit's own `while` loop, which stops advancing `iter` once a
/// match is found in the un-forced case -- see the full two-pass driver
/// below for the one case where a later pass can still override an early
/// hit).
fn ladder_pass(
    eq_table: &BTreeMap<u8, [u8; 4]>,
    tor_type: u8,
    ti: u8,
    tj: u8,
    tk: u8,
    tl: u8,
    fallback_count: &mut usize,
) -> Option<LadderHit> {
    for (stage, (iw, lw)) in LADDER_STAGES.iter().enumerate() {
        let ci = eq_level(eq_table, ti, *iw, fallback_count);
        let cl = eq_level(eq_table, tl, *lw, fallback_count);
        if let Some(v) = raw_lookup_torsion(tor_type, ci, tj, tk, cl) {
            return Some(LadderHit {
                stage: stage as u8,
                key: (tor_type, ci, tj, tk, cl),
                value: v,
            });
        }
    }
    None
}

struct Resolution {
    kind: &'static str, // "exact" | "equivalence_level_1" | "equivalence_level_2" | "equivalence_level_3" | "table_unresolved"
    key: Option<(u8, u8, u8, u8, u8)>,
    value: Option<(f64, f64, f64)>,
    forced_recheck_applied: bool,
}

fn kind_for_stage(stage: u8) -> &'static str {
    match stage {
        0 => "exact",
        1 => "equivalence_level_1",
        2 => "equivalence_level_2",
        3 => "equivalence_level_3",
        _ => unreachable!(),
    }
}

/// Full two-pass RDKit-real torsion-table resolution: primary-code ladder,
/// then (always, if the primary ladder failed entirely, OR -- the one
/// documented empirical quirk in RDKit's own code, `Params.h:841-852` -- if
/// the primary ladder ONLY succeeded at its last, most-wildcarded stage AND
/// `tor_primary == 5` AND a nonzero secondary code exists) a full second
/// 4-stage ladder using the secondary code, which OVERRIDES the primary hit
/// if it also succeeds. Faithful restatement of the C++ `while` loop's
/// control flow (verified by hand-tracing against the pinned source), not a
/// literal transliteration of its unusual `iter`-reset structure.
#[allow(clippy::too_many_arguments)]
fn resolve_torsion(
    eq_table: &BTreeMap<u8, [u8; 4]>,
    tor_primary: u8,
    tor_secondary: u8,
    ti: u8,
    tj: u8,
    tk: u8,
    tl: u8,
    fallback_count: &mut usize,
) -> Resolution {
    let primary_hit = ladder_pass(eq_table, tor_primary, ti, tj, tk, tl, fallback_count);

    let need_secondary = match &primary_hit {
        None => true,
        Some(h) => h.stage == 3 && tor_primary == 5 && tor_secondary != 0,
    };

    if !need_secondary {
        let h = primary_hit.expect("need_secondary is false only when primary_hit is Some");
        return Resolution {
            kind: kind_for_stage(h.stage),
            key: Some(h.key),
            value: Some(h.value),
            forced_recheck_applied: false,
        };
    }

    let secondary_hit = ladder_pass(eq_table, tor_secondary, ti, tj, tk, tl, fallback_count);
    let forced = primary_hit.is_some(); // secondary pass ran despite an existing primary hit
    match secondary_hit {
        Some(h) => Resolution {
            kind: kind_for_stage(h.stage),
            key: Some(h.key),
            value: Some(h.value),
            forced_recheck_applied: forced,
        },
        None => match primary_hit {
            Some(h) => Resolution {
                // Forced recheck found nothing; the earlier primary hit stands.
                kind: kind_for_stage(h.stage),
                key: Some(h.key),
                value: Some(h.value),
                forced_recheck_applied: forced,
            },
            None => Resolution {
                kind: "table_unresolved",
                key: None,
                value: None,
                forced_recheck_applied: false,
            },
        },
    }
}

// ── Main ─────────────────────────────────────────────────────────────────

fn main() {
    let eq_table = load_eq_level_table("scripts/mmff94_provenance/rdkit_defaultMMFFDef.txt");
    eprintln!("eqLevel table: {} atom types loaded", eq_table.len());

    let corpus = load_corpus();
    eprintln!("corpus size: {}", corpus.len());

    let mut candidate_rows: Vec<Value> = Vec::new();
    let mut kind_counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut eq_level_fallback_count = 0usize;
    let mut classification_mismatch_primary = 0usize; // computed_classification != rdkit_classification (primary)
    let mut ladder_resolves_same_code = 0usize; // ladder finds a hit using chematic's OWN classification code (pure ladder win)
    let mut ladder_resolves_only_via_rdkit_code = 0usize; // only resolves because rdkit_classification != computed_classification
    let mut forced_recheck_count = 0usize;
    // Attribution check: chematic's EXISTING mmff94_torsion_energy fallback
    // chain (exact/reverse/2 single-wildcards/double-wildcard/type0-generic
    // -- no eqLevel ladder port involved), fed ONLY the corrected
    // (RDKit-real) classification code, no other change.
    let mut existing_fallback_resolves_with_corrected_code = 0usize;
    let mut existing_fallback_value_matches_ladder = 0usize;
    let mut torsions_total = 0usize;
    let mut torsions_missing = 0usize;
    // Full-population (not just the 1,107 candidates) reachability sweep,
    // same diligence as the sibling stretch-bend PR's false-hit measurement:
    // is the classification-formula disagreement specific to the "missing"
    // population, or does it also silently affect torsions chematic
    // currently finds A value for (possibly the WRONG value, invisible to
    // mmff94_term_coverage_audit.rs, which only logs misses)?
    let mut all_torsions_mismatch = 0usize;
    let mut all_torsions_hit_with_mismatched_code = 0usize;
    let mut all_torsions_hit_with_mismatched_code_different_value = 0usize;
    let mut all_torsions_hit_with_mismatched_code_undetermined = 0usize;

    for cm in &corpus {
        let mol = match chematic_smiles::parse(&cm.smiles) {
            Ok(m) => m,
            Err(_) => continue,
        };
        let types = match assign_mmff94_numeric_types(&mol) {
            Ok(t) => t,
            Err(_) => continue,
        };
        let ring_set = find_sssr(&mol);
        let rings: Vec<Vec<AtomIdx>> = ring_set.rings().to_vec();

        for (_, bond) in mol.bonds() {
            let (j, k) = (bond.atom1, bond.atom2);
            let nbrs_j: Vec<AtomIdx> = mol.neighbors(j).map(|(nb, _)| nb).collect();
            let nbrs_k: Vec<AtomIdx> = mol.neighbors(k).map(|(nb, _)| nb).collect();
            for &i in &nbrs_j {
                if i == k {
                    continue;
                }
                for &l in &nbrs_k {
                    if l == j {
                        continue;
                    }
                    torsions_total += 1;
                    let (ti, tj, tk, tl) = (
                        types[i.0 as usize],
                        types[j.0 as usize],
                        types[k.0 as usize],
                        types[l.0 as usize],
                    );

                    // Same classification chematic-ff's production code
                    // actually computes and feeds into the torsion lookup.
                    let tt_chematic = chematic_ff::torsion_type_for(
                        &rings,
                        i.0 as usize,
                        j.0 as usize,
                        k.0 as usize,
                        l.0 as usize,
                        tj,
                        tk,
                    );

                    // Independently-derived RDKit real classification,
                    // computed for EVERY torsion (not just missing ones) so
                    // the full-population sweep below is possible.
                    let (rdkit_primary, rdkit_secondary) =
                        rdkit_torsion_type(&mol, i, j, k, l, ti, tj, tk, tl);
                    let mismatch = rdkit_primary != tt_chematic;
                    if mismatch {
                        all_torsions_mismatch += 1;
                    }

                    // chematic's OWN fallback-inclusive lookup -- same
                    // criterion mmff94_term_coverage_audit.rs uses to define
                    // "missing".
                    if let Some(chematic_value) = mmff94_torsion_energy(tt_chematic, ti, tj, tk, tl)
                    {
                        if mismatch {
                            all_torsions_hit_with_mismatched_code += 1;
                            // Does the wrong-code hit even land on the SAME
                            // parameter RDKit's real classification's ladder
                            // would pick? Three-way, not binary: our own
                            // table+ladder port might ALSO come up empty
                            // here (table_unresolved -> would need RDKit's
                            // empirical rule, self-port cannot determine the
                            // "correct" value at all) -- that case is
                            // "undetermined" by this diagnostic, NOT
                            // "confirmed different" (advisor-flagged: an
                            // earlier version of this counter conflated the
                            // two, silently inflating the mismatch count
                            // with never-checked rows).
                            let rdkit_res = resolve_torsion(
                                &eq_table,
                                rdkit_primary,
                                rdkit_secondary,
                                ti,
                                tj,
                                tk,
                                tl,
                                &mut eq_level_fallback_count,
                            );
                            let chematic_v =
                                (chematic_value.v1, chematic_value.v2, chematic_value.v3);
                            match rdkit_res.value {
                                Some(v) if v != chematic_v => {
                                    all_torsions_hit_with_mismatched_code_different_value += 1;
                                }
                                Some(_) => {} // ladder confirms chematic's current value is already numerically right, coincidentally
                                None => {
                                    all_torsions_hit_with_mismatched_code_undetermined += 1;
                                }
                            }
                        }
                        continue; // not missing, not the routing-candidate population
                    }
                    torsions_missing += 1;

                    // routing_bug_candidate vs table_gap, same discriminator
                    // as mmff94_term_coverage_audit.rs.
                    let present_at = (0..=8u8)
                        .find(|&code| mmff94_torsion_energy(code, ti, tj, tk, tl).is_some());
                    if present_at.is_none() {
                        continue; // genuine table_gap (14 of these), not this file's subject
                    }

                    if mismatch {
                        classification_mismatch_primary += 1;
                    }

                    // Resolve using chematic's OWN classification code fed
                    // through the real eqLevel ladder (tests the ladder-gap
                    // hypothesis in isolation, holding classification fixed
                    // at whatever chematic already computes).
                    let res_same_code = resolve_torsion(
                        &eq_table,
                        tt_chematic,
                        0, // chematic's torsion_type_for has no secondary-code concept
                        ti,
                        tj,
                        tk,
                        tl,
                        &mut eq_level_fallback_count,
                    );

                    // Resolve using RDKit's OWN classification code (both
                    // primary and secondary) -- the real, full algorithm.
                    let res_rdkit_code = resolve_torsion(
                        &eq_table,
                        rdkit_primary,
                        rdkit_secondary,
                        ti,
                        tj,
                        tk,
                        tl,
                        &mut eq_level_fallback_count,
                    );
                    if res_rdkit_code.forced_recheck_applied {
                        forced_recheck_count += 1;
                    }

                    if res_same_code.kind != "table_unresolved" {
                        ladder_resolves_same_code += 1;
                    } else if res_rdkit_code.kind != "table_unresolved" {
                        ladder_resolves_only_via_rdkit_code += 1;
                    }

                    // Decisive attribution check (advisor-flagged): does
                    // chematic's EXISTING, UNMODIFIED `mmff94_torsion_energy`
                    // fallback chain -- no eqLevel ladder port involved at
                    // all -- already resolve this candidate once fed only
                    // the CORRECTED (RDKit-real) classification code? If so,
                    // the fix is classification-only; the eqLevel ladder
                    // this file ports would be contributing nothing
                    // incremental on this corpus.
                    let existing_fallback_with_corrected_code =
                        mmff94_torsion_energy(rdkit_primary, ti, tj, tk, tl);
                    if let Some(v) = existing_fallback_with_corrected_code {
                        existing_fallback_resolves_with_corrected_code += 1;
                        if res_rdkit_code.value == Some((v.v1, v.v2, v.v3)) {
                            existing_fallback_value_matches_ladder += 1;
                        }
                    }

                    *kind_counts.entry(res_rdkit_code.kind).or_insert(0) += 1;

                    let (used_exact, used_equivalence, used_unresolved) = match res_rdkit_code.kind
                    {
                        "exact" => (true, false, false),
                        "table_unresolved" => (false, false, true),
                        _ => (false, true, false),
                    };

                    candidate_rows.push(json!({
                        "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                        "atoms": [i.0, j.0, k.0, l.0],
                        "atom_types": [ti, tj, tk, tl],
                        "atomic_numbers": [
                            mol.atom(i).element.atomic_number(),
                            mol.atom(j).element.atomic_number(),
                            mol.atom(k).element.atomic_number(),
                            mol.atom(l).element.atomic_number(),
                        ],
                        "computed_classification": tt_chematic,
                        "rdkit_classification": rdkit_primary,
                        "rdkit_classification_secondary": rdkit_secondary,
                        "present_at_different_classification": present_at,
                        "ladder_result_same_code_as_chematic": {
                            "kind": res_same_code.kind,
                            "key": res_same_code.key,
                            "value": res_same_code.value,
                        },
                        "selected_parameter_kind": res_rdkit_code.kind,
                        "selected_parameter_key": res_rdkit_code.key,
                        "selected_parameter_value": res_rdkit_code.value,
                        "forced_recheck_applied": res_rdkit_code.forced_recheck_applied,
                        "used_exact": used_exact,
                        "used_equivalence": used_equivalence,
                        "used_empirical": Value::Null, // determined only by scripts/mmff94_torsion_oracle_validate_227.py (see file doc comment)
                        "used_unresolved": used_unresolved,
                        "note": match res_rdkit_code.kind {
                            "exact" => "RDKit's real algorithm hits the atom-type-exact row at its own (possibly different-from-chematic's) classification code.",
                            "table_unresolved" => "Neither RDKit's real classification's primary nor secondary eqLevel ladder finds a table row -- falls to Halgren's empirical rule under RDKit's real algorithm (chematic has no equivalent of that rule at all); see the oracle-enriched JSONL for whether the empirical rule produces a nonzero term here.",
                            _ => "RDKit's real algorithm resolves this via a genuine equivalence-class (eqLevel) substitution -- a real, fixable classification/routing gap in chematic's simplified wildcard-to-0 fallback.",
                        },
                    }));
                }
            }
        }
    }

    for row in &candidate_rows {
        println!("{row}");
    }

    let total = candidate_rows.len();
    eprintln!("=== mmff94_torsion_equivalence_diagnostic_227 summary ===");
    eprintln!(
        "torsions_total={torsions_total} torsions_missing={torsions_missing} routing_bug_candidate rows examined: {total} (audit's frozen/re-verified count: 1107)"
    );
    eprintln!(
        "eq_level_fallback_count (atom type not found in def table, self-fallback used) = {eq_level_fallback_count}"
    );
    for (k, v) in &kind_counts {
        eprintln!("  selected_parameter_kind[{k}] = {v}");
    }
    eprintln!(
        "classification_mismatch_primary (rdkit_classification != computed_classification) = {classification_mismatch_primary} / {total}"
    );
    eprintln!(
        "ladder_resolves_same_code (eqLevel ladder resolves using CHEMATIC's OWN classification code -- pure ladder-gap fix, no classification fix needed) = {ladder_resolves_same_code}"
    );
    eprintln!(
        "ladder_resolves_only_via_rdkit_code (unresolved at chematic's own code, but RDKit's real classification code's ladder resolves it -- needs a classification fix, not just a ladder) = {ladder_resolves_only_via_rdkit_code}"
    );
    eprintln!(
        "table_unresolved (neither classification's ladder finds anything under RDKit's real algorithm -- falls to empirical rule / genuinely unresolved) = {}",
        kind_counts.get("table_unresolved").copied().unwrap_or(0)
    );
    eprintln!(
        "forced_recheck_count (RDKit's type-5-with-secondary forced-override quirk actually fired) = {forced_recheck_count}"
    );
    eprintln!(
        "existing_fallback_resolves_with_corrected_code (chematic's EXISTING, UNMODIFIED mmff94_torsion_energy fallback chain -- no eqLevel ladder port involved -- fed ONLY the corrected RDKit-real classification code) = {existing_fallback_resolves_with_corrected_code} / {total}"
    );
    eprintln!(
        "existing_fallback_value_matches_ladder (of those, the value matches this file's full eqLevel-ladder port exactly -- confirms whether the ladder is contributing ANYTHING incremental beyond a classification-only fix on this corpus) = {existing_fallback_value_matches_ladder}"
    );
    eprintln!(
        "--- full-population sweep (ALL {torsions_total} torsion instances, not just the 1,107 candidates) ---"
    );
    eprintln!(
        "all_torsions_mismatch (rdkit_classification != computed_classification, any hit/miss status) = {all_torsions_mismatch}"
    );
    eprintln!(
        "all_torsions_hit_with_mismatched_code (chematic's fallback-inclusive lookup currently returns SOME value despite a classification mismatch -- invisible to mmff94_term_coverage_audit.rs, which only logs misses) = {all_torsions_hit_with_mismatched_code}"
    );
    eprintln!(
        "all_torsions_hit_with_mismatched_code_different_value (of those, this file's OWN table+ladder port CONFIRMS a different value than what RDKit's real classification would select -- self-port estimate, NOT oracle-validated for this full-population sweep, only the 1,107-candidate population above was oracle-checked) = {all_torsions_hit_with_mismatched_code_different_value}"
    );
    eprintln!(
        "all_torsions_hit_with_mismatched_code_undetermined (of those, this file's own ladder ALSO comes up empty -- RDKit's real answer would need the empirical rule, this diagnostic cannot determine whether chematic's current value happens to agree or not) = {all_torsions_hit_with_mismatched_code_undetermined}"
    );
}
