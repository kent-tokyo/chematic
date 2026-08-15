//! Phase 1B/1C diagnostic for issue #227's torsion empirical-rule work
//! (Phase 1 of the "raise chematic to 95/100 vs RDKit" roadmap).
//!
//! DIAGNOSTIC ONLY -- touches no production code. Re-derives, for every
//! Torsion instance `mmff94_term_coverage_audit.rs` reports as
//! `torsions_missing` on the 265-molecule Wave 1 corpus
//! (`validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json`),
//! the REAL MMFF bond-type/registry-flag inputs an empirical torsion rule
//! would need -- not a guess from atom-level `aromatic`/`in_ring` context
//! fields (those are chematic's separate ring/aromaticity-perception
//! output, not what `bond_type_for`/`torsion_type_for` actually consume).
//!
//! **Conclusion (issue #227 Phase 1, 2026-08-15, both findings
//! oracle-validated against all 254 real instances, not a sample): NO
//! Halgren empirical torsion rule was implemented.** Two hypotheses this
//! tool was built to test were both falsified:
//! 1. A from-scratch Halgren "aromatic b-c bond" empirical rule (derived
//!    from OpenBabel's `forcefieldmmff94.cpp` comments citing Halgren Part V
//!    pp. 631-632 Table X, since the primary paper is paywalled -- see
//!    `scripts/mmff94_provenance/PROVENANCE.md`) predicts a UNIFORM
//!    `V1=0,V2=6.0,V3=0` for every aromatic-central-bond case in this
//!    corpus. Live oracle: 0/254 match; real values cluster into 7 distinct
//!    tuples that vary with the TERMINAL atoms, which a central-bond-only
//!    empirical formula cannot produce by construction.
//! 2. `eq_level_probe` below (RDKit's real eqLevel canonical-type-
//!    substitution ladder, ported from Angle's already-production
//!    `MMFF94_EQ_LEVEL` table, applied to the terminal ti/tl positions with
//!    tj/tk/tors_type fixed): 0/254 hits.
//!
//! The REAL root cause (confirmed via `other_classification_codes` below,
//! which scans chematic's own unmodified 926-row `MMFF94_TORSION_ENERGY`
//! table at every OTHER classification code 0..=8): 254/254 already have a
//! real row in chematic's table whose value matches the oracle exactly, at
//! a torsion_type code the classification formula didn't reach --
//! because `torsion_type_for` was being fed the WRONG bond order for the
//! central j-k bond (chematic's general aromaticity perception says
//! `BondOrder::Aromatic`; RDKit's real sanitizer Kekulizes the same bond to
//! `Single`/`Double`, oracle-confirmed 254/254 via
//! `GetBondBetweenAtoms(j,k).GetIsAromatic() == False`). Fixed in
//! `crates/chematic-ff/src/mmff94_numeric.rs`
//! (`assign_mmff94_numeric_types_with_view`) -- a classification/bond-order-
//! source fix, not a new resolution tier. See `PROVENANCE.md`'s Torsion
//! entry for the full writeup.
//!
//! Run: `cargo run --release -p chematic-3d --example mmff94_torsion_empirical_diagnostic_227 \
//!   > validation/results/mmff94_torsion_empirical_diagnostic_227.jsonl`

use chematic_core::AtomIdx;
use chematic_ff::{
    assign_mmff94_numeric_types_with_view, bond_type_for, mmff94_numeric_type_info,
    mmff94_torsion_energy, torsion_type_for,
};
use serde_json::{Value, json};

/// TEMPORARY hypothesis probe (not the production table): a byte-for-byte
/// copy of `mmff94_energy::angle::MMFF94_EQ_LEVEL` (the SAME general-purpose
/// MMFF equivalence-class table Angle's Stage B already uses, not an
/// angle-specific one) used here only to test whether Torsion's real
/// resolution ladder needs the identical eqLevel substitution mechanism
/// applied to the TERMINAL atoms (ti/tl), central atoms/tors_type held
/// fixed -- exactly as PROVENANCE.md's Torsion entry already documented as
/// a real (if previously unexercised) RDKit mechanism. If this hypothesis
/// is confirmed against the live oracle, this gets promoted into a shared
/// `pub(crate)` table in `mmff94_energy/mod.rs` instead of being duplicated.
static EQ_LEVEL: &[(u8, [u8; 4])] = &[
    (1, [1, 1, 1, 0]),
    (2, [2, 2, 1, 0]),
    (3, [3, 3, 1, 0]),
    (4, [4, 4, 1, 0]),
    (5, [5, 5, 5, 0]),
    (6, [6, 6, 6, 0]),
    (7, [7, 7, 6, 0]),
    (8, [8, 8, 8, 0]),
    (9, [9, 9, 8, 0]),
    (10, [10, 10, 8, 0]),
    (11, [11, 11, 11, 0]),
    (12, [12, 12, 12, 0]),
    (13, [13, 13, 13, 0]),
    (14, [14, 14, 14, 0]),
    (15, [15, 15, 15, 0]),
    (16, [16, 16, 15, 0]),
    (17, [17, 17, 15, 0]),
    (18, [18, 18, 15, 0]),
    (19, [19, 19, 19, 0]),
    (20, [20, 1, 1, 0]),
    (21, [21, 21, 5, 0]),
    (22, [22, 22, 1, 0]),
    (23, [23, 23, 5, 0]),
    (24, [24, 24, 5, 0]),
    (25, [25, 25, 25, 0]),
    (26, [26, 26, 25, 0]),
    (27, [27, 28, 5, 0]),
    (28, [28, 28, 5, 0]),
    (29, [29, 29, 5, 0]),
    (30, [30, 2, 1, 0]),
    (31, [31, 31, 31, 0]),
    (32, [32, 7, 6, 0]),
    (33, [33, 21, 5, 0]),
    (34, [34, 8, 8, 0]),
    (35, [35, 6, 6, 0]),
    (36, [36, 36, 5, 0]),
    (37, [37, 2, 1, 0]),
    (38, [38, 9, 8, 0]),
    (39, [39, 10, 8, 0]),
    (40, [40, 10, 8, 0]),
    (41, [41, 3, 1, 0]),
    (42, [42, 42, 8, 0]),
    (43, [43, 10, 8, 0]),
    (44, [44, 16, 15, 0]),
    (45, [45, 10, 8, 0]),
    (46, [46, 9, 8, 0]),
    (47, [47, 42, 8, 0]),
    (48, [48, 9, 8, 0]),
    (49, [49, 6, 6, 0]),
    (50, [50, 21, 5, 0]),
    (51, [51, 7, 6, 0]),
    (52, [52, 21, 5, 0]),
    (53, [53, 42, 8, 0]),
    (54, [54, 9, 8, 0]),
    (55, [55, 10, 8, 0]),
];

fn eq_level(t: u8, stage: usize) -> u8 {
    EQ_LEVEL
        .binary_search_by_key(&t, |&(x, _)| x)
        .map(|idx| EQ_LEVEL[idx].1[stage])
        .unwrap_or(t)
}

/// `(stage, (ti, tj, tk, tl) resolved key, params)`.
type EqLevelHit = (usize, (u8, u8, u8, u8), chematic_ff::TorsionEnergyParams);

/// Probe: does substituting ti/tl through the eqLevel ladder (tj/tk/tors_type
/// fixed, matching Angle's real mechanism) find a real table row this
/// diagnostic's own `mmff94_torsion_energy` (exact+reverse+wildcard, no
/// eqLevel) does not? Tries stages 0..4 (Level 3/4/5, Level2==identity
/// already covered by the exact/reverse tiers) and both traversal
/// directions. Returns `(stage, resolved_type_tuple, params)` on first hit.
fn eq_level_probe(tt: u8, ti: u8, tj: u8, tk: u8, tl: u8) -> Option<EqLevelHit> {
    for stage in 0..4usize {
        let (si, sl) = (eq_level(ti, stage), eq_level(tl, stage));
        if let Some(p) = mmff94_torsion_energy(tt, si, tj, tk, sl) {
            return Some((stage, (si, tj, tk, sl), p));
        }
        if let Some(p) = mmff94_torsion_energy(tt, sl, tk, tj, si) {
            return Some((stage, (sl, tk, tj, si), p));
        }
    }
    None
}

fn load_manifest(path: &str) -> Value {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"))
}

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
        let manifest = load_manifest(path);
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

fn atom_flags(t: u8) -> Value {
    match mmff94_numeric_type_info(t) {
        Some(info) => json!({
            "type": t,
            "symbol": info.symbol,
            "element": info.element.symbol(),
            "atomic_number": info.atomic_number,
            "crd": info.coordination,
            "val": info.valence,
            "pilp": info.has_pi_lone_pair,
            "mltb": info.multiple_bond_count,
            "arom": info.aromatic,
            "lin": info.linear,
            "sbmb": info.single_bond_multiple_bond,
        }),
        None => json!({"type": t, "registry_entry": "MISSING"}),
    }
}

fn main() {
    let corpus = load_corpus();
    let mut n_rows = 0usize;
    for cm in &corpus {
        let mol = match chematic_smiles::parse(&cm.smiles) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("PARSE ERROR {}: {e}", cm.name);
                continue;
            }
        };
        let (types, mol) = match assign_mmff94_numeric_types_with_view(&mol) {
            Ok((t, view)) => (t, view),
            Err(e) => {
                eprintln!("TYPING ERROR {}: {e}", cm.name);
                continue;
            }
        };

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
                    let (ti, tj, tk, tl) = (
                        types[i.0 as usize],
                        types[j.0 as usize],
                        types[k.0 as usize],
                        types[l.0 as usize],
                    );
                    let tt = torsion_type_for(
                        &mol,
                        i.0 as usize,
                        j.0 as usize,
                        k.0 as usize,
                        l.0 as usize,
                        ti,
                        tj,
                        tk,
                        tl,
                    );
                    if mmff94_torsion_energy(tt, ti, tj, tk, tl).is_some() {
                        continue; // not missing -- out of scope for this diagnostic
                    }
                    n_rows += 1;

                    let order_ij = mol.bond_between(i, j).expect("i-j bond").1.order;
                    let order_jk = bond.order;
                    let order_kl = mol.bond_between(k, l).expect("k-l bond").1.order;
                    let bt_ij = bond_type_for(ti, tj, order_ij);
                    let bt_jk = bond_type_for(tj, tk, order_jk);
                    let bt_kl = bond_type_for(tk, tl, order_kl);

                    let eq_probe = eq_level_probe(tt, ti, tj, tk, tl).map(|(stage, key, p)| {
                        json!({"stage": stage, "resolved_key": [tt, key.0, key.1, key.2, key.3], "v": [p.v1, p.v2, p.v3]})
                    });

                    // Advisor-directed check: exactly which OTHER classification
                    // code(s) 0..=8 carry a real row for this same (ti,tj,tk,tl)
                    // (or reversed) tuple, and what value do they carry?
                    let mut other_codes: Vec<Value> = Vec::new();
                    for alt in 0..=8u8 {
                        if alt == tt {
                            continue;
                        }
                        if let Some(p) = mmff94_torsion_energy(alt, ti, tj, tk, tl) {
                            other_codes
                                .push(json!({"tt": alt, "dir": "fwd", "v": [p.v1, p.v2, p.v3]}));
                        }
                        if let Some(p) = mmff94_torsion_energy(alt, tl, tk, tj, ti) {
                            other_codes
                                .push(json!({"tt": alt, "dir": "rev", "v": [p.v1, p.v2, p.v3]}));
                        }
                    }

                    let row = json!({
                        "molecule_id": cm.name, "smiles": cm.smiles, "tier": cm.tier,
                        "atoms": [i.0, j.0, k.0, l.0],
                        "atomic_numbers": [
                            mol.atom(i).element.atomic_number(),
                            mol.atom(j).element.atomic_number(),
                            mol.atom(k).element.atomic_number(),
                            mol.atom(l).element.atomic_number(),
                        ],
                        "torsion_type": tt,
                        "flags": [atom_flags(ti), atom_flags(tj), atom_flags(tk), atom_flags(tl)],
                        "bond_order_ij": format!("{:?}", order_ij),
                        "bond_order_jk": format!("{:?}", order_jk),
                        "bond_order_kl": format!("{:?}", order_kl),
                        "bond_type_ij": bt_ij,
                        "bond_type_jk": bt_jk,
                        "bond_type_kl": bt_kl,
                        "eq_level_probe": eq_probe,
                        "other_classification_codes": other_codes,
                    });
                    println!("{row}");
                }
            }
        }
    }
    eprintln!("total missing-torsion rows dumped: {n_rows}");
}
