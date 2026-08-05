//! Phase 1A audit for issue #227 (MMFF94 strict coverage gap).
//!
//! Independently re-measures MMFF94 parameter coverage over the same 265-
//! molecule Wave 1 corpus (`validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json`)
//! used by PR #226, and for every missing bond/angle/torsion/out-of-plane/
//! stretch-bend term, dumps a rich per-term JSONL row: atom-level chemistry
//! context (element, aromaticity, formal charge, ring membership, smallest
//! ring size), the exact classification code chematic-ff computed
//! (bond_type/angle_type/torsion_type), the lookup key before/after
//! chematic-ff's own internal normalization, whether the term was found via
//! a direct hit or a documented fallback tier, and -- the key root-cause
//! discriminator -- whether a row for the same atom-type tuple exists
//! *anywhere* in the underlying table at a *different* classification code
//! (a fixable classification/routing bug) or is absent at every code from
//! 0..=8 (a genuine parameter-table data gap, not a lookup bug).
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
//! *production* `chematic_ff::mmff94_stbn` itself (no longer a
//! diagnostic-only side-table here) -- so `stbn_missing` below already
//! reflects it, and a "would Dfsb resolve this" field would be meaningless
//! (there is no more hypothetical "would", it already does). The
//! `present_at_different_classification` discriminator still uses the
//! type-only lookup (`mmff94_stbn_type_only`), independent of Dfsb, since
//! it specifically asks about classification-code routing bugs.
//!
//! Run: `cargo run --release -p chematic-3d --example mmff94_term_coverage_audit \
//!   > validation/results/mmff94_coverage_227_term_audit.jsonl 2> validation/results/mmff94_coverage_227_stderr.log`

use std::collections::BTreeMap;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_ff::{
    OOP_SP2_TYPES, angle_type_for, assign_mmff94_numeric_types, bond_type_for, mmff94_angle_energy,
    mmff94_bond_energy, mmff94_charges_numeric, mmff94_oop, mmff94_stbn, mmff94_stbn_type_only,
    mmff94_torsion_energy, mmff94_vdw_combined,
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
// without new parameters). If absent at every code 0..=8 -> genuine
// parameter-table data gap (needs new, properly-sourced parameters).

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
    (0..=8u8).find(|&at| mmff94_stbn_type_only(at, ti, tj, tk).is_some())
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
    angles_total: usize,
    angles_missing: usize,
    torsions_total: usize,
    torsions_missing: usize,
    oop_total: usize,
    oop_missing: usize,
    stbn_total: usize,
    stbn_missing: usize,
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
                    angles_total: 0,
                    angles_missing: 0,
                    torsions_total: 0,
                    torsions_missing: 0,
                    oop_total: 0,
                    oop_missing: 0,
                    stbn_total: 0,
                    stbn_missing: 0,
                    vdw_types_total: 0,
                    vdw_types_missing: 0,
                    charges_ok: false,
                    strict_gate_would_fail: true,
                });
                continue;
            }
        };

        let types = match assign_mmff94_numeric_types(&mol) {
            Ok(t) => t,
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
                    angles_total: 0,
                    angles_missing: 0,
                    torsions_total: 0,
                    torsions_missing: 0,
                    oop_total: 0,
                    oop_missing: 0,
                    stbn_total: 0,
                    stbn_missing: 0,
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
        let mut angles_total = 0usize;
        let mut angles_missing = 0usize;
        let mut torsions_total = 0usize;
        let mut torsions_missing = 0usize;
        let mut oop_total = 0usize;
        let mut oop_missing = 0usize;
        let mut stbn_total = 0usize;
        let mut stbn_missing = 0usize;

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
                    "final_lookup_result": "missing",
                    "present_at_different_classification": present_at,
                    "note": "mmff94_bond_energy has NO wildcard/type-0 fallback at all -- any exact (bond_type,ti,tj) miss is unconditionally missing",
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
                            "final_lookup_result": "missing",
                            "present_at_different_classification": present_at,
                        }));
                    }

                    // Priority 2B: mmff94_stbn now includes RDKit's Dfsb
                    // periodic-row fallback in production, so a row only
                    // appears here for a triple genuinely unresolved even
                    // after that fallback (not merely "chematic's
                    // specific/generic table misses" -- that population is
                    // considerably smaller than pre-Priority-2B).
                    let stbn_hit = mmff94_stbn(
                        at,
                        ta,
                        tb,
                        tc,
                        mol.atom(a).element.atomic_number(),
                        mol.atom(b).element.atomic_number(),
                        mol.atom(c).element.atomic_number(),
                    );
                    if stbn_hit.is_none() {
                        stbn_missing += 1;
                        let present_at = stbn_present_at_any_type(ta, tb, tc);
                        term_rows.push(json!({
                            "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                            "term_kind": "StretchBend",
                            "atoms": [ctx(a), ctx(b), ctx(c)],
                            "classified_type": at,
                            "lookup_key_before_normalization": [at, ta, tb, tc],
                            "final_lookup_result": "missing",
                            "present_at_different_classification": present_at,
                            "note": "NEVER gated by ForceFieldPolicy::Mmff94BondAngleStrict's coverage check by default (gate_mmff94_stretch_bend=false) -- silently contributes 0.0 energy in stretch_bend_energy, does not cause a typed failure. Also unresolved by chematic_ff::mmff94_stbn's RDKit-Dfsb fallback (Priority 2B) -- a genuine residual gap under RDKit's own complete stretch-bend algorithm, not just chematic's.",
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
                        &rings,
                        i.0 as usize,
                        j.0 as usize,
                        k.0 as usize,
                        l.0 as usize,
                        tj_,
                        tk_,
                    );
                    let hit = mmff94_torsion_energy(tt, ti_, tj_, tk_, tl_);
                    if hit.is_none() {
                        torsions_missing += 1;
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
                            "final_lookup_result": "missing",
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

        let strict_gate_would_fail = bonds_missing > 0 || angles_missing > 0;

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
            angles_total,
            angles_missing,
            torsions_total,
            torsions_missing,
            oop_total,
            oop_missing,
            stbn_total,
            stbn_missing,
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
            "angles_total": m.angles_total, "angles_missing": m.angles_missing,
            "torsions_total": m.torsions_total, "torsions_missing": m.torsions_missing,
            "oop_total": m.oop_total, "oop_missing": m.oop_missing,
            "stbn_total": m.stbn_total, "stbn_missing": m.stbn_missing,
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
    let n_angles_missing: usize = mol_aggs.iter().map(|m| m.angles_missing).sum();
    let n_torsions_missing: usize = mol_aggs.iter().map(|m| m.torsions_missing).sum();
    let n_oop_missing: usize = mol_aggs.iter().map(|m| m.oop_missing).sum();
    let n_stbn_missing: usize = mol_aggs.iter().map(|m| m.stbn_missing).sum();
    eprintln!(
        "=== summary: total={n_total} bond+angle-gate-would-fail={n_fail} bonds_missing={n_bonds_missing} angles_missing={n_angles_missing} torsions_missing={n_torsions_missing} oop_missing={n_oop_missing} stbn_missing(never gated)={n_stbn_missing} ==="
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
