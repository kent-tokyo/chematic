//! Phase 1A audit for issue #227 (MMFF94 strict coverage gap).
//!
//! Independently re-measures MMFF94 parameter coverage over the same 265-
//! molecule Wave 1 corpus (`validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json`)
//! used by PR #226, and for every missing bond/angle/torsion/out-of-plane/
//! stretch-bend term, dumps a rich per-term JSONL row: atom-level chemistry
//! context (element, aromaticity, formal charge, ring membership, smallest
//! ring size), the exact classification code chematic-ff computed
//! (bond_type/angle_type/torsion_type/stretch_bend_type), the lookup key
//! before/after chematic-ff's own internal normalization, whether the term
//! was found via a direct hit or a documented fallback tier, and -- the key
//! root-cause discriminator -- whether a row for the same atom-type tuple
//! exists *anywhere* in the underlying table at a *different* classification
//! code (a fixable classification/routing bug) or is absent at every code in
//! the term kind's own code space (a genuine parameter-table data gap, not a
//! lookup bug) -- 0..=8 for Bond/Angle/Torsion, 0..=11 for StretchBend
//! (issue #227 Priority 2C: StretchBend's table key is the *stretch-bend
//! type*, not the angle type — see `stretch_bend_type_for`'s doc).
//!
//! vdW and charge coverage are reported at molecule level (these are
//! per-atom-type / whole-molecule lookups, not per n-tuple like the other
//! four term kinds, and chematic-ff's coverage gate for
//! `ForceFieldPolicy::Mmff94BondAngleStrict` never checks them at all --
//! `mmff94_energy_breakdown`'s `stretch_bend_energy`/`vdw_energy` silently
//! contribute 0.0 for any missing term, never erroring. That silent-skip
//! behavior is itself an audit finding, not a bug being fixed here.
//!
//! Priority 2B (issue #227) update: RDKit's periodic-table-row stretch-bend
//! default (`MMFFDfsbCollection::getMMFFDfsbParams`) is now wired into
//! *production* `chematic_ff::mmff94_stbn` itself. This audit deliberately
//! keeps reporting the TYPE-ONLY diagnostic axis (`stbn_missing`,
//! `present_at_different_classification`, using `mmff94_stbn_type_only`)
//! *separately* from the final production resolution
//! (`stbn_final_unresolved`, `dfsb_resolved` per-row) -- coverage parity
//! (does *some* value get returned) and parameter-selection parity (is the
//! *correct* value being used) are different questions. A row whose
//! type-only lookup misses at its own classification code but hits at a
//! *different* one (`present_at_different_classification` is `Some`) is a
//! classification/routing-bug candidate regardless of whether Dfsb then
//! rescues it -- if Dfsb rescues it, that candidate's real, correctly-typed
//! parameter is now *masked* by RDKit's generic periodic-row default
//! instead, not fixed. Collapsing both axes into "0 missing" (an earlier,
//! incorrect version of this file did exactly that) would silently make
//! the 427-instance routing-candidate population undiscoverable from this
//! audit's own output.
//!
//! Run: `cargo run --release -p chematic-3d --example mmff94_term_coverage_audit \
//!   > validation/results/mmff94_coverage_227_term_audit.jsonl 2> validation/results/mmff94_coverage_227_stderr.log`

use std::collections::BTreeMap;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_ff::{
    OOP_SP2_TYPES, angle_type_for, assign_mmff94_numeric_types_with_view, bond_type_for,
    is_angle_in_ring_of_size_3_or_4, mmff94_angle_energy, mmff94_angle_energy_resolved,
    mmff94_bond_energy, mmff94_bond_energy_resolved, mmff94_charges_numeric, mmff94_oop,
    mmff94_stbn, mmff94_stbn_type_only, mmff94_torsion_energy, mmff94_vdw_combined,
    stretch_bend_type_for,
};
use chematic_perception::find_sssr;
use serde_json::{Value, json};

// ── Corpus loading (same manifests as PR #226 / issue #227) ─────────────────

fn load_manifest(path: &str) -> Value {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"))
}

struct CorpusMol {
    tier: String,
    name: String,
    smiles: String,
    primary_category: String,
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
        let manifest = load_manifest(path);
        for m in manifest["molecules"].as_array().expect("molecules array") {
            out.push(CorpusMol {
                tier: tier.to_string(),
                name: m["name"].as_str().unwrap().to_string(),
                smiles: m["smiles"].as_str().unwrap().to_string(),
                primary_category: m["primary_category"]
                    .as_str()
                    .unwrap_or("unknown")
                    .to_string(),
            });
        }
    }
    out
}

// ── Shared per-atom chemistry context ────────────────────────────────────────

fn ring_sizes_containing(rings: &[Vec<AtomIdx>], atom: AtomIdx) -> Vec<usize> {
    rings
        .iter()
        .filter(|r| r.contains(&atom))
        .map(|r| r.len())
        .collect()
}

fn smallest_ring_size(rings: &[Vec<AtomIdx>], atom: AtomIdx) -> Option<usize> {
    ring_sizes_containing(rings, atom).into_iter().min()
}

/// Lightweight derived heuristic (element + bond-order/aromaticity pattern),
/// NOT a stored chematic-core property -- chematic has no `Hybridization`
/// type. Diagnostic label only, not used by any lookup.
fn crude_hybridization(mol: &Molecule, atom: AtomIdx) -> &'static str {
    let a = mol.atom(atom);
    if a.aromatic {
        return "sp2(aromatic)";
    }
    let has_triple = mol
        .neighbors(atom)
        .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Triple);
    let has_double = mol
        .neighbors(atom)
        .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);
    if has_triple {
        "sp"
    } else if has_double {
        "sp2"
    } else {
        "sp3_or_nonC"
    }
}

fn atom_context(mol: &Molecule, rings: &[Vec<AtomIdx>], types: &[u8], a: AtomIdx) -> Value {
    let atom = mol.atom(a);
    json!({
        "index": a.0,
        "element": atom.element.symbol(),
        "aromatic": atom.aromatic,
        "formal_charge": atom.charge,
        "hybridization_heuristic": crude_hybridization(mol, a),
        "in_ring": !ring_sizes_containing(rings, a).is_empty(),
        "smallest_ring_size": smallest_ring_size(rings, a),
        "mmff94_numeric_type": types[a.0 as usize],
    })
}

// ── "present at any classification" scan -------------------------------------
// The key root-cause discriminator: does a row exist for this exact
// atom-type tuple (either order) at *some* classification code the actual
// classifier didn't produce? If yes -> classification/routing bug (fixable
// without new parameters). If absent at every code 0..=8 (Angle/Torsion/
// Bond) -> genuine parameter-table data gap (needs new, properly-sourced
// parameters). StretchBend's own classification code space is 0..=11 (the
// *stretch-bend type*, from `getMMFFStretchBendType` -- see
// `stretch_bend_type_for`'s doc, issue #227 Priority 2C), NOT 0..=8 like the
// other three term kinds -- Angle genuinely uses angle_type 0..=8 as its own
// table key and is unaffected by this distinction.

fn angle_present_at_any_type(ti: u8, tj: u8, tk: u8) -> Option<u8> {
    (0..=8u8).find(|&at| mmff94_angle_energy(at, ti, tj, tk).is_some())
}

fn torsion_present_at_any_type(ti: u8, tj: u8, tk: u8, tl: u8) -> Option<u8> {
    (0..=8u8).find(|&tt| mmff94_torsion_energy(tt, ti, tj, tk, tl).is_some())
}

fn bond_present_at_any_type(ti: u8, tj: u8) -> Option<u8> {
    (0..=8u8).find(|&bt| mmff94_bond_energy(bt, ti, tj).is_some())
}

fn stbn_present_at_any_type(ti: u8, tj: u8, tk: u8) -> Option<u8> {
    (0..=11u8).find(|&sbt| mmff94_stbn_type_only(sbt, ti, tj, tk).is_some())
}

// ── Main audit -----------------------------------------------------------------

struct MolAgg {
    tier: String,
    name: String,
    smiles: String,
    category: String,
    parse_ok: bool,
    typing_ok: bool,
    typing_error: Option<String>,
    bonds_total: usize,
    bonds_missing: usize,
    bonds_final_unresolved: usize,
    angles_total: usize,
    angles_missing: usize,
    angles_final_unresolved: usize,
    torsions_total: usize,
    torsions_missing: usize,
    /// Issue #227 Phase 1: torsions where RDKit itself generates no term at
    /// all (linear central atom, `chematic_ff::torsion_no_term_by_design`)
    /// -- not a coverage gap. Included in `torsions_total`, excluded from
    /// `torsions_missing`.
    torsions_no_term_by_design: usize,
    oop_total: usize,
    oop_missing: usize,
    stbn_total: usize,
    stbn_missing: usize,
    stbn_final_unresolved: usize,
    vdw_types_total: usize,
    vdw_types_missing: usize,
    charges_ok: bool,
    strict_gate_would_fail: bool,
}

fn main() {
    let corpus = load_corpus();
    eprintln!("corpus size: {}", corpus.len());

    let mut term_rows: Vec<Value> = Vec::new();
    let mut mol_aggs: Vec<MolAgg> = Vec::new();

    for cm in &corpus {
        let mol = match chematic_smiles::parse(&cm.smiles) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("PARSE ERROR {}: {e}", cm.name);
                mol_aggs.push(MolAgg {
                    tier: cm.tier.clone(),
                    name: cm.name.clone(),
                    smiles: cm.smiles.clone(),
                    category: cm.primary_category.clone(),
                    parse_ok: false,
                    typing_ok: false,
                    typing_error: Some(format!("parse_failure: {e}")),
                    bonds_total: 0,
                    bonds_missing: 0,
                    bonds_final_unresolved: 0,
                    angles_total: 0,
                    angles_missing: 0,
                    angles_final_unresolved: 0,
                    torsions_total: 0,
                    torsions_missing: 0,
                    torsions_no_term_by_design: 0,
                    oop_total: 0,
                    oop_missing: 0,
                    stbn_total: 0,
                    stbn_missing: 0,
                    stbn_final_unresolved: 0,
                    vdw_types_total: 0,
                    vdw_types_missing: 0,
                    charges_ok: false,
                    strict_gate_would_fail: true,
                });
                continue;
            }
        };

        let (types, mol) = match assign_mmff94_numeric_types_with_view(&mol) {
            // Shadow `mol` with the MMFF-specific re-perceived view (issue
            // #227 Phase 1, torsion parameter gap root cause): every
            // classification call below (`bond_type_for`/`angle_type_for`/
            // `torsion_type_for`/`stretch_bend_type_for`, all of which read
            // `BondOrder` directly) must see the SAME bond orders the
            // numeric types were derived from, not chematic's general/SMILES
            // aromaticity perception -- see
            // `assign_mmff94_numeric_types_with_view`'s doc for why. Same
            // atom count/topology as the original `mol`, so every other use
            // below (ring detection, atom iteration, `ctx`) is unaffected
            // except where `BondOrder`/`atom.aromatic` is read.
            Ok((t, view)) => (t, view),
            Err(e) => {
                eprintln!("TYPING ERROR {}: {e}", cm.name);
                mol_aggs.push(MolAgg {
                    tier: cm.tier.clone(),
                    name: cm.name.clone(),
                    smiles: cm.smiles.clone(),
                    category: cm.primary_category.clone(),
                    parse_ok: true,
                    typing_ok: false,
                    typing_error: Some(e.to_string()),
                    bonds_total: 0,
                    bonds_missing: 0,
                    bonds_final_unresolved: 0,
                    angles_total: 0,
                    angles_missing: 0,
                    angles_final_unresolved: 0,
                    torsions_total: 0,
                    torsions_missing: 0,
                    torsions_no_term_by_design: 0,
                    oop_total: 0,
                    oop_missing: 0,
                    stbn_total: 0,
                    stbn_missing: 0,
                    stbn_final_unresolved: 0,
                    vdw_types_total: 0,
                    vdw_types_missing: 0,
                    charges_ok: false,
                    strict_gate_would_fail: true,
                });
                continue;
            }
        };

        let ring_set = find_sssr(&mol);
        let rings: Vec<Vec<AtomIdx>> = ring_set.rings().to_vec();

        let ctx = |a: AtomIdx| atom_context(&mol, &rings, &types, a);

        let mut bonds_total = 0usize;
        let mut bonds_missing = 0usize;
        let mut bonds_final_unresolved = 0usize;
        let mut angles_total = 0usize;
        let mut angles_missing = 0usize;
        let mut angles_final_unresolved = 0usize;
        let mut torsions_total = 0usize;
        let mut torsions_missing = 0usize;
        let mut torsions_no_term_by_design = 0usize;
        let mut oop_total = 0usize;
        let mut oop_missing = 0usize;
        let mut stbn_total = 0usize;
        let mut stbn_missing = 0usize;
        let mut stbn_final_unresolved = 0usize;

        // -- Bond --
        for (_, bond) in mol.bonds() {
            bonds_total += 1;
            let (a1, a2) = (bond.atom1, bond.atom2);
            let (t1, t2) = (types[a1.0 as usize], types[a2.0 as usize]);
            let bt = bond_type_for(t1, t2, bond.order);
            let (key_lo, key_hi) = if t1 <= t2 { (t1, t2) } else { (t2, t1) };
            let hit = mmff94_bond_energy(bt, t1, t2);
            if hit.is_none() {
                bonds_missing += 1;
                let present_at = bond_present_at_any_type(t1, t2);
                // Issue #227 Stage C: does the empirical rule (eq. 18-19)
                // resolve what the type-only table lookup above did not?
                let empirical_resolved = mmff94_bond_energy_resolved(bt, t1, t2).is_some();
                if !empirical_resolved {
                    bonds_final_unresolved += 1;
                }
                term_rows.push(json!({
                    "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                    "term_kind": "Bond",
                    "atoms": [ctx(a1), ctx(a2)],
                    "bond_order": format!("{:?}", bond.order),
                    "classified_type": bt,
                    "lookup_key_before_normalization": [bt, t1, t2],
                    "lookup_key_after_normalization": [bt, key_lo, key_hi],
                    "direct_table_hit": false,
                    "fallback_available": false,
                    "fallback_hit": false,
                    "final_lookup_result": if empirical_resolved { "resolved_via_empirical_rule" } else { "missing" },
                    "present_at_different_classification": present_at,
                    "empirical_resolved": empirical_resolved,
                    "note": "mmff94_bond_energy has NO wildcard/type-0 fallback at all -- any exact (bond_type,ti,tj) miss is unconditionally missing at the TYPE-ONLY level; empirical_resolved reports the Stage C eq.18-19 fallback's own, separate outcome",
                }));
            }
        }

        // -- Angle (+ Stretch-Bend, same triples, never gated by the strict policy) --
        for b_idx in 0..mol.atom_count() {
            let b = AtomIdx(b_idx as u32);
            let neighbors: Vec<AtomIdx> = mol.neighbors(b).map(|(nb, _)| nb).collect();
            if neighbors.len() < 2 {
                continue;
            }
            for i in 0..neighbors.len() {
                for j in (i + 1)..neighbors.len() {
                    let (a, c) = (neighbors[i], neighbors[j]);
                    angles_total += 1;
                    stbn_total += 1;
                    let (ta, tc) = (types[a.0 as usize], types[c.0 as usize]);
                    let tb = types[b_idx];
                    let at =
                        angle_type_for(&mol, &rings, a.0 as usize, b_idx, c.0 as usize, &types);

                    // Shared with the Angle check below AND the StretchBend
                    // check further down -- both need the two flanking bonds'
                    // types (issue #227 Priority 2C: StretchBend's table key
                    // is the *stretch-bend type*, not the angle type `at` --
                    // see `stretch_bend_type_for`'s doc). Hoisted here so both
                    // checks share one computation, matching the now-fixed
                    // production call sites in chematic-ff's
                    // mmff94_minimizer.rs / chematic-3d's own
                    // `compute_mmff94_coverage` exactly.
                    let order_ab = mol.bond_between(a, b).expect("a-b angle bond").1.order;
                    let order_cb = mol.bond_between(c, b).expect("c-b angle bond").1.order;
                    let bt_ab = bond_type_for(ta, tb, order_ab);
                    let bt_cb = bond_type_for(tc, tb, order_cb);

                    let angle_hit = mmff94_angle_energy(at, ta, tb, tc);
                    if angle_hit.is_none() {
                        angles_missing += 1;
                        let present_at = angle_present_at_any_type(ta, tb, tc);
                        // direct vs fallback: type-0 fallback is only tried
                        // when at != 0 (mirrors mmff94_angle_energy's own logic).
                        let type0_hit = if at != 0 {
                            mmff94_angle_energy(0, ta, tb, tc).is_some()
                        } else {
                            false
                        };
                        // Issue #227 Stage C: does the empirical rule (eq.
                        // 20) resolve what the type-only table lookup above
                        // did not? Needs both flanking bonds' r0 -- if
                        // either is itself unresolvable (extremely rare,
                        // only for elements missing from
                        // MMFF94_COV_RAD_PAU_ELE/MMFF94_HERSCHBACH_LAURIE),
                        // the angle term is also left unresolved, matching
                        // RDKit's own real `getMMFFAngleBendParams` (which
                        // requires both flanking `getMMFFBondStretchParams`
                        // calls to succeed before even attempting empirical).
                        let bond_ab = mmff94_bond_energy_resolved(bt_ab, ta, tb);
                        let bond_cb = mmff94_bond_energy_resolved(bt_cb, tc, tb);
                        let empirical_resolved = match (bond_ab, bond_cb) {
                            (Some((bab, _)), Some((bcb, _))) => {
                                let ring_size = is_angle_in_ring_of_size_3_or_4(
                                    &mol,
                                    a.0 as usize,
                                    b_idx,
                                    c.0 as usize,
                                );
                                mmff94_angle_energy_resolved(
                                    at, ta, tb, tc, bab.r0, bcb.r0, ring_size,
                                )
                                .is_some()
                            }
                            _ => false,
                        };
                        if !empirical_resolved {
                            angles_final_unresolved += 1;
                        }
                        term_rows.push(json!({
                            "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                            "term_kind": "Angle",
                            "atoms": [ctx(a), ctx(b), ctx(c)],
                            "classified_type": at,
                            "lookup_key_before_normalization": [at, ta, tb, tc],
                            "lookup_key_after_normalization": [at, ta.min(tc), tb, ta.max(tc)],
                            "direct_table_hit": false,
                            "fallback_available": at != 0,
                            "fallback_hit": type0_hit,
                            "final_lookup_result": if empirical_resolved { "resolved_via_empirical_rule" } else { "missing" },
                            "present_at_different_classification": present_at,
                            "empirical_resolved": empirical_resolved,
                        }));
                    }

                    let sbt = stretch_bend_type_for(at, ta, tc, bt_ab, bt_cb);

                    // Review-driven fix (Priority 2B follow-up): a row must
                    // be emitted whenever the TYPE-ONLY lookup misses,
                    // regardless of whether the Dfsb fallback then rescues
                    // it -- rows are the only place `present_at_different_
                    // classification` (a classification/routing-bug
                    // candidate, independent of Dfsb) is visible at all. The
                    // earlier version of this file only emitted a row when
                    // the FINAL (Dfsb-inclusive) lookup missed, which
                    // silently dropped the 427/2,107 type-routing
                    // candidates that Dfsb happens to also rescue --
                    // masking, not fixing, that population. See
                    // `dfsb_resolved` below to distinguish "Dfsb rescued a
                    // routing-bug candidate" (parameter-selection parity
                    // still open, tracked separately) from "Dfsb rescued a
                    // genuine table gap" (the only case actually closed).
                    let type_only_hit = mmff94_stbn_type_only(sbt, ta, tb, tc);
                    if type_only_hit.is_none() {
                        stbn_missing += 1;
                        let present_at = stbn_present_at_any_type(ta, tb, tc);
                        let final_hit = mmff94_stbn(
                            sbt,
                            ta,
                            tb,
                            tc,
                            mol.atom(a).element.atomic_number(),
                            mol.atom(b).element.atomic_number(),
                            mol.atom(c).element.atomic_number(),
                        );
                        let dfsb_resolved = final_hit.is_some();
                        if !dfsb_resolved {
                            stbn_final_unresolved += 1;
                        }
                        term_rows.push(json!({
                            "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                            "term_kind": "StretchBend",
                            "atoms": [ctx(a), ctx(b), ctx(c)],
                            "classified_type": sbt,
                            "angle_type": at,
                            "lookup_key_before_normalization": [sbt, ta, tb, tc],
                            "final_lookup_result": if dfsb_resolved { "resolved_via_dfsb_fallback" } else { "missing" },
                            "present_at_different_classification": present_at,
                            "dfsb_resolved": dfsb_resolved,
                            "note": if dfsb_resolved {
                                if present_at.is_some() {
                                    "Type-only lookup missed at this triple's own classification code, but a row EXISTS at a different code (present_at_different_classification) -- a classification/routing-bug candidate, NOT a genuine table gap. chematic_ff::mmff94_stbn's RDKit-Dfsb fallback (Priority 2B) resolved this triple anyway (coverage parity achieved) -- but that means it is using RDKit's GENERIC periodic-row default, not the SPECIFIC parameter a correctly-routed classification would have used (parameter-selection parity still open). Never gated by ForceFieldPolicy::Mmff94BondAngleStrict's coverage check by default (gate_mmff94_stretch_bend=false)."
                                } else {
                                    "Absent at every classification code (a genuine type-table gap, not a routing-bug candidate) -- resolved by chematic_ff::mmff94_stbn's RDKit-Dfsb fallback (Priority 2B), matching RDKit's own real behavior exactly (this is the case Dfsb was designed to close). Never gated by ForceFieldPolicy::Mmff94BondAngleStrict's coverage check by default (gate_mmff94_stretch_bend=false)."
                                }
                            } else {
                                "NEVER gated by ForceFieldPolicy::Mmff94BondAngleStrict's coverage check by default (gate_mmff94_stretch_bend=false) -- silently contributes 0.0 energy in stretch_bend_energy, does not cause a typed failure. Also unresolved by chematic_ff::mmff94_stbn's RDKit-Dfsb fallback (Priority 2B) -- a genuine residual gap under RDKit's own complete stretch-bend algorithm, not just chematic's."
                            },
                        }));
                    }
                }
            }
        }

        // -- Torsion --
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
                    let (ti_, tj_, tk_, tl_) = (
                        types[i.0 as usize],
                        types[j.0 as usize],
                        types[k.0 as usize],
                        types[l.0 as usize],
                    );
                    let tt = chematic_ff::torsion_type_for(
                        &mol,
                        i.0 as usize,
                        j.0 as usize,
                        k.0 as usize,
                        l.0 as usize,
                        ti_,
                        tj_,
                        tk_,
                        tl_,
                    );
                    let hit = mmff94_torsion_energy(tt, ti_, tj_, tk_, tl_);
                    if hit.is_none() {
                        // Issue #227 Phase 1: a linear central atom means
                        // RDKit itself generates no torsion term here either
                        // (chematic_ff::torsion_no_term_by_design) -- correct
                        // behavior, not a coverage gap. Reported separately
                        // so it never inflates torsions_missing.
                        let no_term_by_design = chematic_ff::torsion_no_term_by_design(tj_, tk_);
                        if no_term_by_design {
                            torsions_no_term_by_design += 1;
                        } else {
                            torsions_missing += 1;
                        }
                        let present_at = torsion_present_at_any_type(ti_, tj_, tk_, tl_);
                        term_rows.push(json!({
                            "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                            "term_kind": "Torsion",
                            "atoms": [ctx(i), ctx(j), ctx(k), ctx(l)],
                            "classified_type": tt,
                            "lookup_key_before_normalization": [tt, ti_, tj_, tk_, tl_],
                            "direct_table_hit": false,
                            "fallback_available": true,
                            "fallback_hit": false,
                            "final_lookup_result": if no_term_by_design { "no_term_by_design" } else { "missing" },
                            "present_at_different_classification": present_at,
                            "note": "mmff94_torsion_energy already tries exact+reverse+2 single-wildcards+double-wildcard+type0-generic (7 tiers) before returning None",
                        }));
                    }
                }
            }
        }

        // -- Out-of-plane --
        for j_idx in 0..mol.atom_count() {
            let tj = types[j_idx];
            if OOP_SP2_TYPES.binary_search(&tj).is_err() {
                continue;
            }
            let j = AtomIdx(j_idx as u32);
            let neighbors: Vec<AtomIdx> = mol.neighbors(j).map(|(nb, _)| nb).collect();
            if neighbors.len() != 3 {
                continue;
            }
            oop_total += 1;
            let [i, k, l] = [neighbors[0], neighbors[1], neighbors[2]];
            let hit = mmff94_oop(
                tj,
                types[i.0 as usize],
                types[k.0 as usize],
                types[l.0 as usize],
            );
            if hit.is_none() {
                oop_missing += 1;
                term_rows.push(json!({
                    "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                    "term_kind": "Oop",
                    "atoms": [ctx(j), ctx(i), ctx(k), ctx(l)],
                    "final_lookup_result": "missing",
                    "note": "mmff94_oop already tries 6 substituent orderings + 4 wildcard tiers before returning None",
                }));
            }
        }

        // -- vdW (per distinct numeric type present in this molecule) --
        let mut distinct_types: Vec<u8> = types.clone();
        distinct_types.sort_unstable();
        distinct_types.dedup();
        let mut vdw_types_total = 0usize;
        let mut vdw_types_missing = 0usize;
        for &t in &distinct_types {
            vdw_types_total += 1;
            if mmff94_vdw_combined(t, t).is_none() {
                vdw_types_missing += 1;
                term_rows.push(json!({
                    "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                    "term_kind": "VanDerWaals",
                    "mmff94_numeric_type": t,
                    "final_lookup_result": "missing",
                    "note": "NEVER gated by ForceFieldPolicy::Mmff94BondAngleStrict's coverage check",
                }));
            }
        }

        // -- Charges (whole-molecule) --
        let charges_ok = mmff94_charges_numeric(&mol).is_ok();
        if !charges_ok {
            term_rows.push(json!({
                "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                "term_kind": "Charge",
                "final_lookup_result": "missing",
                "note": "mmff94_charges_numeric failed for the whole molecule",
            }));
        }

        // Issue #227 Stage C: the REAL strict gate (chematic-3d's
        // `compute_mmff94_coverage`, and the actual minimizer) now resolves
        // through the eq.18-20 empirical rule too -- so what actually
        // predicts gate failure is the FINAL (post-empirical) unresolved
        // count, not the type-only-table miss count. `bonds_missing`/
        // `angles_missing` stay as the type-only diagnostic axis (same
        // "coverage parity vs. parameter-selection parity" distinction this
        // file already draws for StretchBend/Dfsb -- see the module doc).
        let strict_gate_would_fail = bonds_final_unresolved > 0 || angles_final_unresolved > 0;

        mol_aggs.push(MolAgg {
            tier: cm.tier.clone(),
            name: cm.name.clone(),
            smiles: cm.smiles.clone(),
            category: cm.primary_category.clone(),
            parse_ok: true,
            typing_ok: true,
            typing_error: None,
            bonds_total,
            bonds_missing,
            bonds_final_unresolved,
            angles_total,
            angles_missing,
            angles_final_unresolved,
            torsions_total,
            torsions_missing,
            torsions_no_term_by_design,
            oop_total,
            oop_missing,
            stbn_total,
            stbn_missing,
            stbn_final_unresolved,
            vdw_types_total,
            vdw_types_missing,
            charges_ok,
            strict_gate_would_fail,
        });
    }

    // Emit per-term rows.
    for row in &term_rows {
        println!("{row}");
    }

    // Emit per-molecule aggregate rows (prefixed so they can be filtered
    // separately from term rows in the same JSONL stream).
    for m in &mol_aggs {
        let row = json!({
            "row_type": "molecule_summary",
            "tier": m.tier, "molecule_id": m.name, "smiles": m.smiles, "category": m.category,
            "parse_ok": m.parse_ok, "typing_ok": m.typing_ok, "typing_error": m.typing_error,
            "bonds_total": m.bonds_total, "bonds_missing": m.bonds_missing,
            "bonds_final_unresolved": m.bonds_final_unresolved,
            "angles_total": m.angles_total, "angles_missing": m.angles_missing,
            "angles_final_unresolved": m.angles_final_unresolved,
            "torsions_total": m.torsions_total, "torsions_missing": m.torsions_missing,
            "torsions_no_term_by_design": m.torsions_no_term_by_design,
            "oop_total": m.oop_total, "oop_missing": m.oop_missing,
            "stbn_total": m.stbn_total, "stbn_missing": m.stbn_missing,
            "stbn_final_unresolved": m.stbn_final_unresolved,
            "vdw_types_total": m.vdw_types_total, "vdw_types_missing": m.vdw_types_missing,
            "charges_ok": m.charges_ok,
            "strict_bond_angle_gate_would_fail": m.strict_gate_would_fail,
        });
        println!("{row}");
    }

    // Aggregate summary to stderr for a quick human sanity check.
    let n_total = mol_aggs.len();
    let n_fail = mol_aggs.iter().filter(|m| m.strict_gate_would_fail).count();
    let n_bonds_missing: usize = mol_aggs.iter().map(|m| m.bonds_missing).sum();
    let n_bonds_final_unresolved: usize = mol_aggs.iter().map(|m| m.bonds_final_unresolved).sum();
    let n_angles_missing: usize = mol_aggs.iter().map(|m| m.angles_missing).sum();
    let n_angles_final_unresolved: usize = mol_aggs.iter().map(|m| m.angles_final_unresolved).sum();
    let n_torsions_missing: usize = mol_aggs.iter().map(|m| m.torsions_missing).sum();
    let n_torsions_no_term_by_design: usize =
        mol_aggs.iter().map(|m| m.torsions_no_term_by_design).sum();
    let n_oop_missing: usize = mol_aggs.iter().map(|m| m.oop_missing).sum();
    let n_stbn_missing: usize = mol_aggs.iter().map(|m| m.stbn_missing).sum();
    let n_stbn_final_unresolved: usize = mol_aggs.iter().map(|m| m.stbn_final_unresolved).sum();
    eprintln!(
        "=== summary: total={n_total} bond+angle-gate-would-fail={n_fail} bonds_missing(type-only)={n_bonds_missing} bonds_final_unresolved(after Stage C empirical)={n_bonds_final_unresolved} angles_missing(type-only)={n_angles_missing} angles_final_unresolved(after Stage C empirical)={n_angles_final_unresolved} torsions_missing={n_torsions_missing} torsions_no_term_by_design(linear central atom, RDKit also has none)={n_torsions_no_term_by_design} oop_missing={n_oop_missing} stbn_type_only_missing(never gated, incl. Dfsb-masked routing candidates)={n_stbn_missing} stbn_final_unresolved(after Dfsb fallback)={n_stbn_final_unresolved} ==="
    );

    // Distinct missing-tuple pattern counts (angle), to see concentration.
    let mut angle_patterns: BTreeMap<(u8, u8, u8, u8), usize> = BTreeMap::new();
    for row in &term_rows {
        if row["term_kind"] == "Angle" {
            let key = row["lookup_key_before_normalization"].as_array().unwrap();
            let k = (
                key[0].as_u64().unwrap() as u8,
                key[1].as_u64().unwrap() as u8,
                key[2].as_u64().unwrap() as u8,
                key[3].as_u64().unwrap() as u8,
            );
            *angle_patterns.entry(k).or_insert(0) += 1;
        }
    }
    eprintln!(
        "distinct missing angle (type,i,j,k) patterns: {}",
        angle_patterns.len()
    );
    let mut sorted: Vec<_> = angle_patterns.into_iter().collect();
    sorted.sort_by_key(|b| std::cmp::Reverse(b.1));
    for (k, count) in sorted.iter().take(20) {
        eprintln!("  {k:?}: {count}");
    }
}
