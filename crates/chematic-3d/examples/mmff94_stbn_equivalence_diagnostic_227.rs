//! DIAGNOSTIC-ONLY (issue #227). Does NOT touch any production code path.
//!
//! Investigates the 427 StretchBend `routing_bug_candidate` instances found
//! by `mmff94_term_coverage_audit.rs` (a real row exists in chematic-ff's
//! `MMFF94_STBN` table for the exact atom-type triple, just not at the
//! classification code chematic-ff's production code computed for it).
//!
//! The task that produced this file started from a hypothesis: MMFF94 has a
//! "canonical type equivalence" fallback ladder (RDKit's `eqLevel`
//! mechanism, used by angle/torsion/OOP) that chematic-ff is missing, and
//! that a correctly-implemented equivalence ladder would resolve most of
//! the 427 stretch-bend routing candidates. **That hypothesis is false for
//! stretch-bend, and this file demonstrates why from RDKit's real source**
//! (`Code/ForceField/MMFF/Params.h`'s `MMFFStbnCollection::getMMFFStbnParams`,
//! pinned commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`, same pin as
//! `scripts/mmff94_provenance/PROVENANCE.md`): it does exactly ONE exact
//! lookup (after I/K canonicalization), then RDKit falls straight to the
//! periodic-row `Dfsb` default. No `eqLevel` step exists anywhere in this
//! path -- `PROVENANCE.md`'s stretch-bend row already said so before this
//! file was written, based on the same source read.
//!
//! **The real, demonstrated root cause** is a distinct, much more concrete
//! bug: chematic-ff's production stretch-bend code (`mmff94_stbn_type_only`
//! / `mmff94_stbn` in `crates/chematic-ff/src/mmff94_energy/oop_stbn.rs`,
//! called from `mmff94_minimizer.rs`) uses the **angle type** (0-8, from
//! `angle_type_for`) directly as the `MMFF94_STBN` table key. RDKit's real
//! algorithm (`AtomTyper.cpp`'s `MMFFMolProperties::getMMFFStretchBendParams`,
//! lines ~3566-3612 at the pinned commit) computes a *different*,
//! finer-grained "stretch-bend type" (0-11) via a dedicated function,
//! `getMMFFStretchBendType(angleType, bondType1, bondType2)`
//! (`AtomTyper.cpp:2480-2508`), which further splits angle types 1, 5, 7
//! into two stretch-bend types apiece depending on which of the angle's two
//! flanking bonds individually has MMFF bond-type 1 (`angle_type_for`'s
//! `bt_sum` only records the *sum*, discarding exactly the information this
//! split needs). chematic-ff never computes this second code at all --
//! it queries `MMFF94_STBN` with the angle type standing in for it.
//!
//! **Self-consistency proof this is real, with no oracle needed**: `MMFF94_STBN`'s
//! own frozen data (the ported `defaultMMFFStbn` table, verbatim from the
//! pinned commit) contains exactly one row at key 5: `(5, 22, 22, 22, 0.0,
//! 0.0)`, all type 22 (`CR3R`, cyclopropane ring carbon). Cyclopropane's
//! ring C-C-C angle has both flanking bonds single/non-aromatic/all-CR3R,
//! so both individual MMFF bond types are 0 -> `bt_sum=0` -> under
//! chematic's *own* documented ring-offset table
//! (`mmff94_minimizer.rs::angle_type_for`'s doc comment), a 3-ring with
//! `bt_sum=0` is angle type **3**, not 5. If the `MMFF94_STBN` key column
//! really were "angle type", a row keyed `5` for an all-CR3R triple would
//! be unreachable garbage -- no all-CR3R angle can ever compute
//! `angle_type_for(..) == 5` (chematic's own table: 3-ring/bt_sum=0 -> 3,
//! not 5; 5 is reserved for bt_sum=1). Under RDKit's real
//! `getMMFFStretchBendType`, angle type 3 (3-ring, bt_sum=0) maps to
//! stretch-bend type **5** -- exactly this row. The 11 key-4 rows are all
//! `CR4R`-centred (type 20) triples for the same reason (angle type 4,
//! 4-ring bt_sum=0, maps to stretch-bend type 4 -- the one case where the
//! numeral is unchanged). This is decisive: the column is `stretchBendType`,
//! not `angleType`, provable from the shipped data alone.
//!
//! **A second, independent, real bug found while tracing this**:
//! `angle_type_for`'s own ring-offset arithmetic for the `bt_sum=2` case
//! disagrees with RDKit's real formula (`AtomTyper.cpp`'s
//! `getMMFFAngleType`: `angleType = size; if (bondTypeSum) angleType +=
//! (bondTypeSum + size - 2);`). For a 3-ring: RDKit gives bt_sum=0->3,
//! bt_sum=1->5, bt_sum=2->**6**; chematic's own doc table gives bt_sum=2->
//! **8**. For a 4-ring: RDKit gives bt_sum=0->4, bt_sum=1->**7**,
//! bt_sum=2->**8**; chematic's doc table gives bt_sum=1->6, bt_sum=2->7.
//! This affects `Angle`'s own 277 routing candidates directly (same table
//! key space, 0-8) and *compounds* into `StretchBend`'s bug (a wrong angle
//! type feeds a wrong `getMMFFStretchBendType` call). This file's
//! `rdkit_angle_type` reproduces RDKit's real formula, kept deliberately
//! separate from `rdkit_stretch_bend_type` so the two bugs' contributions
//! stay distinguishable in the output (`angle_type_chematic` vs
//! `angle_type_rdkit` vs `stretch_bend_type_rdkit`).
//!
//! Also ported here, faithfully separate from chematic-ff's SSSR-based
//! ring perception: RDKit's `isAngleInRingOfSize3or4` (`AtomTyper.cpp:357`)
//! does *not* use SSSR at all -- it's a purely local bond-adjacency check
//! (3-ring: i-k directly bonded; 4-ring: i and k, excluding j, share a
//! common neighbour). `rdkit_ring_size_3_or_4` below ports that exact
//! check, independent of chematic-ff's `find_sssr`-based
//! `atoms_share_ring_of_size`, so any divergence between the two mechanisms
//! is visible in the output (`ring_size_chematic` vs `ring_size_rdkit`)
//! rather than silently assumed away.
//!
//! Everything RDKit-derived in this file is a *fresh, independent* port
//! (including its own 29-row copy of `defaultMMFFDfsb`, not a call into
//! chematic-ff's private `mmff94_dfsb_stbn`) -- deliberately not reusing
//! chematic-ff's stretch-bend classification code at all, so this
//! diagnostic cannot be silently "correct by construction" from sharing the
//! bug it's trying to measure. The one piece of chematic-ff logic reused
//! as-is is `bond_type_for` (individual bond type index, i-j and j-k
//! separately); see the live-oracle results below for exactly how its
//! *output*, given chematic's own bond-order/aromaticity *input*, checks
//! out.
//!
//! **Live RDKit oracle validation (this environment's `rdkit==2026.3.3`,
//! matching the pinned source commit's release tag)**: contrary to an
//! earlier draft of this comment, RDKit's Python binding DOES expose
//! per-term classification codes directly --
//! `MMFFMolProperties.GetMMFFStretchBendParams(mol, i, j, k)` returns
//! `(stretchBendType, kbaIJK, kbaKJI)` or `None`, and
//! `GetMMFFBondStretchParams`/`GetMMFFAngleBendParams` do the analogous
//! thing for bond/angle terms -- found by enumerating `dir()` on the
//! `MMFFMolProperties` object, not assumed. `scripts/mmff94_stbn_oracle_validate_227.py`
//! (same pinned RDKit build) reads this file's own JSONL output and calls
//! `GetMMFFStretchBendParams` for every one of the 427 candidates, using
//! the SAME molecule/atom-index triple. Result, run against the actual
//! 265-molecule corpus:
//!
//! - **255/427 (59.7%) have chematic's OWN bond-order/aromaticity
//!   perception (feeding `bond_type_for`) agreeing with RDKit's real,
//!   post-sanitization bond typing for both flanking bonds.** On this
//!   "clean" subset, this file's `rdkit_classification` /
//!   `selected_parameter_kind` / `selected_parameter_value` match the live
//!   oracle EXACTLY -- 228/228 stretch-bend-type matches (100%, zero
//!   exceptions), including all 27 `found_but_zero_dropped` rows (oracle
//!   returns `None` for every one, confirming RDKit really does drop a
//!   found-but-(0.0,0.0) row rather than falling through to Dfsb). Of the
//!   255: 220 are `exact` (chematic's current output is the generic Dfsb
//!   default; RDKit's real answer is a specific row), 27 are
//!   `found_but_zero_dropped` (RDKit has NO stretch-bend term here at all;
//!   chematic currently injects a nonzero Dfsb value), 8 are
//!   `dfsb_fallback` (RDKit ALSO falls to Dfsb -- chematic's current output
//!   is already numerically correct today, coincidentally).
//! - **172/427 (40.3%, 60 distinct molecules) have chematic's bond typing
//!   DISAGREEING with RDKit's for at least one of the two flanking bonds.**
//!   Root cause (traced by hand on the first offending triple,
//!   `chembl_tier_b_0000` atoms 6-7-8): the input SMILES writes a ring as
//!   lowercase (`n2ncc(=O)[nH]c2=O`, a pyridazine-3,6-dione-like tautomer);
//!   chematic's SMILES parser trusts lowercase input directly (see
//!   `CLAUDE.md`: "Aromatic SMILES atoms ... set `atom.aromatic = true`
//!   directly during parsing"), but RDKit's sanitizer rejects it (fails
//!   Hückel 4n+2 with two exocyclic carbonyls pulling ring electron
//!   density) and kekulizes it to a genuine alternating single/double
//!   system, so the SAME i-j bond is `BondOrder::Aromatic` in chematic
//!   (`bond_type_for` returns 0 unconditionally for `Aromatic`, per that
//!   function's own doc comment) but a real, individually-typed `SINGLE`
//!   bond between two `sbmb`-flagged atoms in RDKit (bond type 1). This is
//!   an upstream, pre-existing aromaticity-perception gap between the two
//!   engines (matching the already-documented, partial MMFF aromaticity
//!   port in `PROVENANCE.md`'s "Priority 1A" row) -- NOT a bug in this
//!   file's `getMMFFStretchBendType`/`getMMFFAngleType` port, and not
//!   fixable by correcting the stretch-bend classification alone: even a
//!   perfect classifier fed the wrong bond order still computes the wrong
//!   answer. For all 172, the live oracle confirms RDKit resolves to a
//!   SPECIFIC row every time (`stretchBendType` in `{1, 2, 3}`, never
//!   `None`, never a Dfsb-shaped generic value) -- so chematic's current
//!   Dfsb-masked output is almost certainly ALSO wrong for these 172, just
//!   for a different, out-of-scope reason.
//! - **These same 172 (molecule, atom-triple) instances are, set-for-set
//!   IDENTICAL (172 = 172, zero symmetric difference, verified directly),
//!   to the (molecule, triple) instances flagged as routing candidates in
//!   BOTH `Angle` (277 total) and `StretchBend` (427 total)
//!   simultaneously** -- a real, confirmed shared root cause between the
//!   two populations (both consume the same `angle_type_for` output, which
//!   is itself downstream of the same misperceived bond order), but the
//!   shared mechanism is this aromaticity-perception gap, not a missing
//!   equivalence-class fallback ladder.
//!
//! A separate real bug, found while tracing this and independently
//! measured (see `angle_offset_bug_contributes` / `all_triples_*` in the
//! stderr summary): `angle_type_for`'s ring-offset arithmetic for the
//! `bt_sum=2` (3-ring) and `bt_sum in {1,2}` (4-ring) cases disagrees with
//! RDKit's real `getMMFFAngleType` formula (documented above,
//! self-consistency-proven from the shipped `MMFF94_STBN` data). This bug
//! is REAL but measured as LATENT on this 265-molecule corpus: of 10,107
//! total angle triples, 113 are ring-embedded 3-/4-membered angles (both
//! ring-detection mechanisms agree), and 0/113 hit the diverging
//! `bt_sum>=1` branches -- every ring-embedded angle triple in this corpus
//! happens to have `bt_sum=0` (the one case where chematic's table and
//! RDKit's formula still agree). It contributes NOTHING measured to either
//! the 427 or the 277 populations here; filed as a separate, independently
//! real but unproven-impact finding, not part of this diagnostic's
//! resolution count.
//!
//! A small, separate "false hit" population, invisible to
//! `mmff94_term_coverage_audit.rs` entirely (it only logs misses): of 8,000
//! triples where chematic's CURRENT `mmff94_stbn_type_only` already returns
//! a value (not "missing"), 4 use a classification code that disagrees with
//! this file's `rdkit_classification`, and all 4 of those also return a
//! numerically different `(kbaIJK, kbaKJI)`. Tiny, but a real silent-wrong-
//! parameter class distinct from the 427 "missing" population -- worth a
//! follow-up issue, not chased further here.
//!
//! Run:
//! ```text
//! cargo run --release -p chematic-3d --example mmff94_stbn_equivalence_diagnostic_227 \
//!   > validation/results/mmff94_stbn_equivalence_diagnostic_227.jsonl \
//!   2> validation/results/mmff94_stbn_equivalence_diagnostic_227_stderr.log
//! ```
//! Oracle validation:
//! ```text
//! .venv/bin/python scripts/mmff94_stbn_oracle_validate_227.py \
//!   validation/results/mmff94_stbn_equivalence_diagnostic_227.jsonl
//! ```

use chematic_core::{AtomIdx, Molecule};
use chematic_ff::mmff94_energy::MMFF94_STBN;
use chematic_ff::{
    angle_type_for, assign_mmff94_numeric_types, bond_type_for, mmff94_stbn_type_only,
};
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

// ── RDKit-faithful ports (independent of chematic-ff's production code) ────

/// Port of RDKit's `isAngleInRingOfSize3or4` (`AtomTyper.cpp:357-395`,
/// pinned commit). Purely local bond-adjacency, NOT SSSR-based: 3-ring iff
/// i-k are directly bonded; 4-ring iff i and k (excluding j) share a common
/// neighbour. Returns 0 if neither.
fn rdkit_ring_size_3_or_4(mol: &Molecule, i: AtomIdx, j: AtomIdx, k: AtomIdx) -> u8 {
    if mol.bond_between(i, k).is_some() {
        return 3;
    }
    let nbrs_i: std::collections::BTreeSet<AtomIdx> = mol
        .neighbors(i)
        .map(|(n, _)| n)
        .filter(|&n| n != j)
        .collect();
    let has_common = mol
        .neighbors(k)
        .map(|(n, _)| n)
        .filter(|&n| n != j)
        .any(|n| nbrs_i.contains(&n));
    if has_common { 4 } else { 0 }
}

/// Port of RDKit's `getMMFFAngleType` (`AtomTyper.cpp:2412-2447`).
/// `angleType = size; if (bondTypeSum) angleType += (bondTypeSum + size - 2)`.
fn rdkit_angle_type(bond_type_sum: u8, ring_size: u8) -> u8 {
    if ring_size == 0 {
        return bond_type_sum;
    }
    if bond_type_sum == 0 {
        ring_size
    } else {
        ring_size + (bond_type_sum + ring_size - 2)
    }
}

/// Port of RDKit's `getMMFFStretchBendType` (`AtomTyper.cpp:2480-2508`).
fn rdkit_stretch_bend_type(angle_type: u8, bond_type1: u8, bond_type2: u8) -> u8 {
    match angle_type {
        1 => {
            if bond_type1 != 0 || bond_type1 == bond_type2 {
                1
            } else {
                2
            }
        }
        2 => 3,
        4 => 4,
        3 => 5,
        5 => {
            if bond_type1 != 0 || bond_type1 == bond_type2 {
                6
            } else {
                7
            }
        }
        6 => 8,
        7 => {
            if bond_type1 != 0 || bond_type1 == bond_type2 {
                9
            } else {
                10
            }
        }
        8 => 11,
        _ => 0,
    }
}

/// `(selected_stbn_row_key, (kba_ijk, kba_kji))`.
type StbnLookupResult = ((u8, u8, u8, u8), (f64, f64));

/// Port of RDKit's `MMFFStbnCollection::getMMFFStbnParams` (`Params.h:601-`
/// `663`): single exact lookup into the type table, no fallback ladder, no
/// equivalence-class step. `bond_type_ij`/`bond_type_jk` here are the RAW
/// (unreordered) individual bond types, matching RDKit's own call
/// (`bondType[0]`, `bondType[1]`, not the reordered args used to compute
/// `stretch_bend_type`). Returns `(swap, (kba_ijk, kba_kji))` in i,j,k
/// order (swap already applied), or `None` on a genuine miss.
fn rdkit_stbn_exact_lookup(
    stretch_bend_type: u8,
    bond_type_ij: u8,
    bond_type_jk: u8,
    ti: u8,
    tj: u8,
    tk: u8,
) -> Option<StbnLookupResult> {
    let swap = if ti > tk {
        true
    } else if ti == tk {
        bond_type_ij < bond_type_jk
    } else {
        false
    };
    let (can_i, can_k) = if swap { (tk, ti) } else { (ti, tk) };
    let idx = MMFF94_STBN
        .binary_search_by_key(
            &(stretch_bend_type, can_i, tj, can_k),
            |&(a, i, j, k, _, _)| (a, i, j, k),
        )
        .ok()?;
    let (a, i, j, k, kba_ijk, kba_kji) = MMFF94_STBN[idx];
    let (out_ijk, out_kji) = if swap {
        (kba_kji, kba_ijk)
    } else {
        (kba_ijk, kba_kji)
    };
    Some(((a, i, j, k), (out_ijk, out_kji)))
}

/// Independent 29-row copy of RDKit's `defaultMMFFDfsb`
/// (`scripts/mmff94_provenance/rdkit_defaultMMFFDfsb.txt`), NOT a call into
/// chematic-ff's private `mmff94_dfsb_stbn` -- this diagnostic must not
/// share code with the production path it is cross-checking.
const RDKIT_DFSB: &[(u8, u8, u8, f64, f64)] = &[
    (0, 1, 0, 0.15, 0.15),
    (0, 1, 1, 0.10, 0.30),
    (0, 1, 2, 0.05, 0.35),
    (0, 1, 3, 0.05, 0.35),
    (0, 1, 4, 0.05, 0.35),
    (0, 2, 0, 0.00, 0.00),
    (0, 2, 1, 0.00, 0.15),
    (0, 2, 2, 0.00, 0.15),
    (0, 2, 3, 0.00, 0.15),
    (0, 2, 4, 0.00, 0.15),
    (1, 1, 1, 0.30, 0.30),
    (1, 1, 2, 0.30, 0.50),
    (1, 1, 3, 0.30, 0.50),
    (1, 1, 4, 0.30, 0.50),
    (2, 1, 2, 0.50, 0.50),
    (2, 1, 3, 0.50, 0.50),
    (2, 1, 4, 0.50, 0.50),
    (3, 1, 3, 0.50, 0.50),
    (3, 1, 4, 0.50, 0.50),
    (4, 1, 4, 0.50, 0.50),
    (1, 2, 1, 0.30, 0.30),
    (1, 2, 2, 0.25, 0.25),
    (1, 2, 3, 0.25, 0.25),
    (1, 2, 4, 0.25, 0.25),
    (2, 2, 2, 0.25, 0.25),
    (2, 2, 3, 0.25, 0.25),
    (2, 2, 4, 0.25, 0.25),
    (3, 2, 3, 0.25, 0.25),
    (3, 2, 4, 0.25, 0.25),
];

/// Port of RDKit's `getPeriodicTableRow` (`AtomTyper.cpp:251-264`).
fn rdkit_periodic_table_row(atomic_number: u8) -> u8 {
    match atomic_number {
        3..=10 => 1,
        11..=18 => 2,
        19..=36 => 3,
        37..=54 => 4,
        _ => 0,
    }
}

/// `(selected_dfsb_row_key, (f_ijk, f_kji))`.
type DfsbLookupResult = ((u8, u8, u8), (f64, f64));

/// Port of RDKit's `MMFFDfsbCollection::getMMFFDfsbParams` (`Params.h:690-`
/// `714`): swap iff `row_i > row_k` (simple, no tie-break -- unlike the
/// type-table lookup). `RDKit's own `isDoubleZero(kbaIJK) &&
/// isDoubleZero(kbaKJI)` exclusion applied by the caller
/// (`getMMFFStretchBendParams`) is replicated here too.
fn rdkit_dfsb_lookup(z_i: u8, z_j: u8, z_k: u8) -> Option<DfsbLookupResult> {
    let (row_i, row_j, row_k) = (
        rdkit_periodic_table_row(z_i),
        rdkit_periodic_table_row(z_j),
        rdkit_periodic_table_row(z_k),
    );
    let swap = row_i > row_k;
    let (can_i, can_k) = if swap { (row_k, row_i) } else { (row_i, row_k) };
    let &(ri, rj, rk, f_ijk, f_kji) = RDKIT_DFSB
        .iter()
        .find(|&&(r1, r2, r3, ..)| r1 == can_i && r2 == row_j && r3 == can_k)?;
    if f_ijk == 0.0 && f_kji == 0.0 {
        return None;
    }
    let (out_ijk, out_kji) = if swap { (f_kji, f_ijk) } else { (f_ijk, f_kji) };
    Some(((ri, rj, rk), (out_ijk, out_kji)))
}

// ── Full RDKit-faithful resolution for one StretchBend triple ───────────────

struct RdkitResolution {
    bond_type_ij: u8,
    bond_type_jk: u8,
    angle_type_rdkit: u8,
    stretch_bend_type_rdkit: u8,
    outcome: &'static str, // "exact" | "found_but_zero_dropped" | "dfsb_fallback" | "unresolved"
    selected_key: Option<Value>,
    selected_value: Option<(f64, f64)>,
}

#[allow(clippy::too_many_arguments)]
fn resolve_rdkit(
    mol: &Molecule,
    a: AtomIdx,
    b: AtomIdx,
    c: AtomIdx,
    ta: u8,
    tb: u8,
    tc: u8,
    za: u8,
    zb: u8,
    zc: u8,
) -> RdkitResolution {
    let order_ab = mol.bond_between(a, b).expect("angle bond a-b").1.order;
    let order_bc = mol.bond_between(b, c).expect("angle bond b-c").1.order;
    let bond_type_ij = bond_type_for(ta, tb, order_ab);
    let bond_type_jk = bond_type_for(tb, tc, order_bc);
    let ring_size_rdkit = rdkit_ring_size_3_or_4(mol, a, b, c);
    let angle_type_rdkit = rdkit_angle_type(bond_type_ij + bond_type_jk, ring_size_rdkit);

    // RDKit's exact arg-selection for getMMFFStretchBendType, canonicalized
    // on atomType[0] (=ta) vs atomType[2] (=tc) -- NOT the same swap rule
    // used inside getMMFFStbnParams itself (that one also tie-breaks on
    // bond type when ta==tc). See AtomTyper.cpp:3598-3600.
    let (sbt_arg1, sbt_arg2) = if ta <= tc {
        (
            bond_type_ij,
            if ta < tc { bond_type_jk } else { bond_type_ij },
        )
    } else {
        (bond_type_jk, bond_type_ij)
    };
    let stretch_bend_type_rdkit = rdkit_stretch_bend_type(angle_type_rdkit, sbt_arg1, sbt_arg2);

    if let Some((key, (v1, v2))) = rdkit_stbn_exact_lookup(
        stretch_bend_type_rdkit,
        bond_type_ij,
        bond_type_jk,
        ta,
        tb,
        tc,
    ) {
        if v1 == 0.0 && v2 == 0.0 {
            return RdkitResolution {
                bond_type_ij,
                bond_type_jk,
                angle_type_rdkit,
                stretch_bend_type_rdkit,
                outcome: "found_but_zero_dropped",
                selected_key: Some(json!(key)),
                selected_value: None,
            };
        }
        return RdkitResolution {
            bond_type_ij,
            bond_type_jk,
            angle_type_rdkit,
            stretch_bend_type_rdkit,
            outcome: "exact",
            selected_key: Some(json!(key)),
            selected_value: Some((v1, v2)),
        };
    }

    if let Some((key, val)) = rdkit_dfsb_lookup(za, zb, zc) {
        return RdkitResolution {
            bond_type_ij,
            bond_type_jk,
            angle_type_rdkit,
            stretch_bend_type_rdkit,
            outcome: "dfsb_fallback",
            selected_key: Some(json!(key)),
            selected_value: Some(val),
        };
    }

    RdkitResolution {
        bond_type_ij,
        bond_type_jk,
        angle_type_rdkit,
        stretch_bend_type_rdkit,
        outcome: "unresolved",
        selected_key: None,
        selected_value: None,
    }
}

// ── Main ─────────────────────────────────────────────────────────────────

fn main() {
    let corpus = load_corpus();
    eprintln!("corpus size: {}", corpus.len());

    let mut candidate_rows: Vec<Value> = Vec::new();
    let mut outcome_counts: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();
    let mut chematic_already_matches = 0usize; // dfsb_fallback AND chematic's current output is numerically identical
    let mut chematic_currently_wrong = 0usize; // exact OR found_but_zero_dropped (a real fixable discrepancy)
    let mut angle_offset_bug_contributes = 0usize; // computed angle_type != rdkit angle_type, same ring-membership verdict
    let mut ring_detection_disagrees = 0usize; // chematic ring membership (via SSSR) vs RDKit local-adjacency disagree

    // False-hit sweep: among ALL angle/stbn triples (not just the 427),
    // how many does chematic's type-only lookup currently return a value
    // for (a "hit") using the WRONG code relative to RDKit's real
    // stretch_bend_type? These never appear in the 427 (they're not
    // "missing"), so mmff94_term_coverage_audit.rs cannot see them.
    let mut total_stbn_hits = 0usize;
    let mut false_hits_wrong_code = 0usize;
    let mut false_hits_wrong_value = 0usize;

    // Unconditional (ALL angle triples in the corpus, regardless of hit or
    // miss status, StretchBend or not) reachability tally for the
    // `angle_type_for` ring-offset bug and for SSSR-vs-local-adjacency ring
    // detection disagreement -- settles whether either mechanism is latent
    // (never exercised by this corpus) or actually contributes, instead of
    // inferring from the StretchBend-only counters above.
    let mut all_triples = 0usize;
    let mut all_triples_ring_agrees_3or4 = 0usize;
    let mut all_triples_angle_type_diverges_ring_agrees = 0usize;

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

        for b_idx in 0..mol.atom_count() {
            let b = AtomIdx(b_idx as u32);
            let neighbors: Vec<AtomIdx> = mol.neighbors(b).map(|(nb, _)| nb).collect();
            if neighbors.len() < 2 {
                continue;
            }
            for i in 0..neighbors.len() {
                for j in (i + 1)..neighbors.len() {
                    let (a, c) = (neighbors[i], neighbors[j]);
                    let (ta, tc) = (types[a.0 as usize], types[c.0 as usize]);
                    let tb = types[b_idx];
                    let (za, zb, zc) = (
                        mol.atom(a).element.atomic_number(),
                        mol.atom(b).element.atomic_number(),
                        mol.atom(c).element.atomic_number(),
                    );

                    // Same classification chematic-ff's production code
                    // actually computes and feeds into the Stbn lookup.
                    let at_chematic =
                        angle_type_for(&mol, &rings, a.0 as usize, b_idx, c.0 as usize, &types);

                    let chematic_hit = mmff94_stbn_type_only(at_chematic, ta, tb, tc);

                    // Reachability measurement for the angle_type_for
                    // ring-offset bug, independent of whether this triple
                    // is a StretchBend routing candidate at all.
                    let ring_size_chematic_bit = if rings
                        .iter()
                        .any(|r| r.len() == 3 && [a, b, c].iter().all(|x| r.contains(x)))
                    {
                        3
                    } else if rings
                        .iter()
                        .any(|r| r.len() == 4 && [a, b, c].iter().all(|x| r.contains(x)))
                    {
                        4
                    } else {
                        0
                    };
                    let ring_size_rdkit_bit = rdkit_ring_size_3_or_4(&mol, a, b, c);
                    if ring_size_chematic_bit != ring_size_rdkit_bit {
                        ring_detection_disagrees += 1;
                    }

                    // Unconditional tally (every angle triple, hit or miss;
                    // `ring_detection_disagrees` above already covers the
                    // disagreement case for all triples, not just misses).
                    all_triples += 1;
                    if ring_size_chematic_bit == ring_size_rdkit_bit
                        && (ring_size_rdkit_bit == 3 || ring_size_rdkit_bit == 4)
                    {
                        all_triples_ring_agrees_3or4 += 1;
                        let order_ab = mol.bond_between(a, b).expect("a-b bond").1.order;
                        let order_bc = mol.bond_between(b, c).expect("b-c bond").1.order;
                        let bt_ij = bond_type_for(ta, tb, order_ab);
                        let bt_jk = bond_type_for(tb, tc, order_bc);
                        let angle_type_rdkit_all =
                            rdkit_angle_type(bt_ij + bt_jk, ring_size_rdkit_bit);
                        if at_chematic != angle_type_rdkit_all {
                            all_triples_angle_type_diverges_ring_agrees += 1;
                        }
                    }

                    if let Some((v1, v2)) = chematic_hit {
                        total_stbn_hits += 1;
                        let res = resolve_rdkit(&mol, a, b, c, ta, tb, tc, za, zb, zc);
                        if at_chematic != res.stretch_bend_type_rdkit {
                            false_hits_wrong_code += 1;
                            let matches_value = match (res.outcome, res.selected_value) {
                                ("exact", Some((r1, r2))) => (r1, r2) == (v1, v2),
                                ("found_but_zero_dropped", None) => false, // RDKit intends no term at all
                                ("dfsb_fallback", Some((r1, r2))) => (r1, r2) == (v1, v2),
                                _ => false,
                            };
                            if !matches_value {
                                false_hits_wrong_value += 1;
                            }
                        }
                        continue;
                    }

                    // Not a hit -- is it a routing_bug_candidate (row exists
                    // at SOME code 0..=8) or a genuine table_gap? Only
                    // routing candidates are this file's subject (matches
                    // mmff94_term_coverage_audit.rs's own definition).
                    let present_at =
                        (0..=8u8).find(|&at| mmff94_stbn_type_only(at, ta, tb, tc).is_some());
                    if present_at.is_none() {
                        continue; // genuine table_gap, not our population
                    }

                    let res = resolve_rdkit(&mol, a, b, c, ta, tb, tc, za, zb, zc);
                    *outcome_counts.entry(res.outcome).or_insert(0) += 1;

                    let ring_agrees = ring_size_chematic_bit == ring_size_rdkit_bit;
                    if ring_agrees && at_chematic != res.angle_type_rdkit {
                        angle_offset_bug_contributes += 1;
                    }

                    let (used_exact, used_dfsb) = match res.outcome {
                        "exact" | "found_but_zero_dropped" => (true, false),
                        "dfsb_fallback" => (false, true),
                        _ => (false, false),
                    };

                    // "Already numerically correct by accident": chematic's
                    // CURRENT production mmff94_stbn output for this exact
                    // (routing-candidate) triple already falls to Dfsb
                    // today (that's what makes it a routing_bug_candidate
                    // that Dfsb rescues in the first place -- see
                    // mmff94_term_coverage_audit.rs's `dfsb_resolved`); if
                    // the RDKit-correct answer ALSO lands on Dfsb with the
                    // same atomic numbers, chematic's current output is
                    // already numerically identical to the correct one --
                    // fixing the classification code wouldn't change the
                    // energy for this specific triple.
                    if res.outcome == "dfsb_fallback" {
                        chematic_already_matches += 1;
                    } else {
                        chematic_currently_wrong += 1;
                    }

                    candidate_rows.push(json!({
                        "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                        "atoms": [a.0, b.0, c.0],
                        "atom_types": [ta, tb, tc],
                        "atomic_numbers": [za, zb, zc],
                        "ring_size_chematic": if ring_size_chematic_bit == 0 { Value::Null } else { json!(ring_size_chematic_bit) },
                        "ring_size_rdkit": if ring_size_rdkit_bit == 0 { Value::Null } else { json!(ring_size_rdkit_bit) },
                        "bond_type_ij": res.bond_type_ij,
                        "bond_type_jk": res.bond_type_jk,
                        "computed_classification": at_chematic,
                        "angle_type_rdkit": res.angle_type_rdkit,
                        "rdkit_classification": res.stretch_bend_type_rdkit,
                        "present_at_different_classification_chematic_space": present_at,
                        "selected_parameter_kind": res.outcome,
                        "selected_parameter_key": res.selected_key,
                        "selected_parameter_value": res.selected_value.map(|(x, y)| json!([x, y])),
                        "used_exact": used_exact,
                        "used_equivalence": false,
                        "used_generic": false,
                        "used_dfsb": used_dfsb,
                        "note": match res.outcome {
                            "exact" => "RDKit's real algorithm hits a SPECIFIC MMFF94_STBN row at the correctly-derived stretch-bend type -- chematic's current production output for this triple uses the generic Dfsb periodic-row default instead. A real, fixable discrepancy.",
                            "found_but_zero_dropped" => "RDKit's real algorithm finds a row at the correct stretch-bend type whose kbaIJK/kbaKJI are BOTH zero -- RDKit treats this as 'no stretch-bend term' (isDoubleZero&&isDoubleZero short-circuits before Dfsb is ever tried). chematic's current production output supplies a nonzero Dfsb value instead -- RDKit intends zero contribution here, chematic currently contributes a nonzero one. A real, fixable discrepancy in the other direction.",
                            "dfsb_fallback" => "RDKit's real algorithm ALSO falls through to the periodic-row Dfsb default for this triple (no row exists at the correctly-derived stretch-bend type either) -- same fallback chematic's current (buggy-classification) code already reaches. Numerically already correct today by coincidence, not because the classification is right.",
                            _ => "Neither a specific row nor a Dfsb row exists under RDKit's real algorithm either -- genuinely unresolvable even with a fully correct classification (would need new parameter data, not a routing fix).",
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
    eprintln!("=== mmff94_stbn_equivalence_diagnostic_227 summary ===");
    eprintln!("routing_bug_candidate rows examined: {total} (audit's frozen count: 427)");
    for (k, v) in &outcome_counts {
        eprintln!("  outcome[{k}] = {v}");
    }
    eprintln!(
        "chematic_currently_wrong (exact OR found_but_zero_dropped -- a real, fixable discrepancy) = {chematic_currently_wrong}"
    );
    eprintln!(
        "chematic_already_matches (both land on Dfsb, numerically identical today) = {chematic_already_matches}"
    );
    eprintln!(
        "angle_offset_bug_contributes (ring-membership agrees, but angle_type_for's own bt_sum=2 ring-offset arithmetic disagrees with RDKit's formula) = {angle_offset_bug_contributes}"
    );
    eprintln!(
        "ring_detection_disagrees (chematic SSSR-based vs RDKit local-adjacency ring-of-3/4 verdict differs, ALL angle/stbn triples in corpus) = {ring_detection_disagrees}"
    );
    eprintln!(
        "all_triples = {all_triples}, all_triples_ring_agrees_3or4 (both mechanisms agree it's a 3- or 4-ring angle) = {all_triples_ring_agrees_3or4}, all_triples_angle_type_diverges_ring_agrees (of those, angle_type_for's own output != RDKit's correct formula's output) = {all_triples_angle_type_diverges_ring_agrees}"
    );
    eprintln!(
        "--- false-hit sweep (triples chematic's type-only lookup currently HITS, i.e. NOT in the 427 -- invisible to mmff94_term_coverage_audit.rs) ---"
    );
    eprintln!(
        "total_stbn_hits (all triples where mmff94_stbn_type_only(at_chematic,...) is Some) = {total_stbn_hits}"
    );
    eprintln!(
        "false_hits_wrong_code (at_chematic != stretch_bend_type_rdkit among those hits) = {false_hits_wrong_code}"
    );
    eprintln!(
        "false_hits_wrong_value (of those, chematic's returned (kba_ijk,kba_kji) also numerically differs from RDKit's real resolution) = {false_hits_wrong_value}"
    );
}
