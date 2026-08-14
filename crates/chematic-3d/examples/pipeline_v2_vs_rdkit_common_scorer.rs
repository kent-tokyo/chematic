//! Wave 1 follow-up: a single, independent geometry + stereo scorer applied
//! identically to BOTH chematic's and RDKit's already-saved heavy-atom
//! coordinates (`validation/results/pipeline_v2_vs_rdkit_{chematic,rdkit}_rows.jsonl`,
//! produced by the Wave 1 dump/oracle -- neither is re-run here).
//!
//! Why a separate scorer, not `pipeline_v2`'s own `final_validation.sound`:
//! that field is the PRODUCTION pipeline's internal judgment, only ever
//! computed for chematic's own output. Comparing chematic's "sound" against
//! nothing on the RDKit side is not a comparison. This scorer:
//!   - re-derives ideal bond length from `Element::covalent_radius()` (public
//!     API), never touching pipeline_v2's own `pub(crate)` constants/thresholds
//!     or its private `compute_final_validation` -- fully independent, and
//!     the exact same formula the Wave 1 dump already used for the legacy arm.
//!   - is applied to EVERY successful row from EITHER engine, heavy atoms
//!     only (RDKit's saved coords are already heavy-atom-only; chematic's
//!     model has no explicit hydrogens to begin with).
//!   - also runs `chematic_3d::stereo_constraints::verify_stereo` -- the
//!     SAME judge chematic's own pipeline uses internally -- against RDKit's
//!     geometry, using the molecule reparsed from the same SMILES (atom
//!     mapping already verified 265/265 in the Wave 1 aggregate).
//!
//! Pipeline-internal judgment (`sound`, `final_stereo_*`) from the Wave 1
//! dump is NOT overwritten -- this script's output is a separate JSONL,
//! joined against it downstream by the report generator, so "independent
//! benchmark judgment" and "pipeline's own internal judgment" stay two
//! distinct, never-conflated fields.
//!
//! No pipeline_v2/production algorithm code is touched by this file.
//!
//! Run: `cargo run --release -p chematic-3d --example pipeline_v2_vs_rdkit_common_scorer
//!   > validation/results/pipeline_v2_vs_rdkit_common_scored_rows.jsonl`

use chematic_3d::coords::{Coords3D, Point3};
use chematic_3d::stereo_constraints::verify_stereo;
use chematic_core::{AtomIdx, Molecule};
use serde_json::{Value, json};
use std::collections::HashMap;

const CLASH_THRESHOLD_ANGSTROM: f64 = 1.2;
const SOUND_MAX_BOND_RATIO: f64 = 3.0;

struct GeometryCheck {
    atom_count_expected: usize,
    atom_count_actual: usize,
    atom_count_match: bool,
    all_finite: bool,
    coincident_atom_pairs: usize,
    worst_bond_length_ratio: f64,
    bond_violation_rate_15pct: f64,
    bond_violation_rate_50pct: f64,
    gross_clash_count: usize,
    independently_sound: bool,
}

fn coords_from_json(coords_json: &Value, n: usize) -> Option<Coords3D> {
    let arr = coords_json.as_array()?;
    if arr.len() != n {
        return None;
    }
    let mut coords = Coords3D::new_zeroed(n);
    for (i, entry) in arr.iter().enumerate() {
        let triple = entry.as_array()?;
        if triple.len() != 3 {
            return None;
        }
        let x = triple[0].as_f64()?;
        let y = triple[1].as_f64()?;
        let z = triple[2].as_f64()?;
        coords.set(AtomIdx(i as u32), Point3::new(x, y, z));
    }
    Some(coords)
}

fn independent_geometry_check(mol: &Molecule, coords: &Coords3D) -> GeometryCheck {
    let n = mol.atom_count();
    let atom_count_actual = coords.atom_count();
    let atom_count_match = atom_count_actual == n;
    let all_finite = coords.is_finite();

    let mut worst_ratio = 0.0f64;
    let mut violations_15 = 0usize;
    let mut violations_50 = 0usize;
    let bond_count = mol.bond_count().max(1);
    if atom_count_match {
        for (_bidx, bond) in mol.bonds() {
            let a = bond.atom1;
            let b = bond.atom2;
            let actual = coords.get(a).distance(&coords.get(b));
            let ideal = (mol.atom(a).element.covalent_radius() as f64
                + mol.atom(b).element.covalent_radius() as f64)
                .max(0.3);
            let rel_error = (actual / ideal - 1.0).abs();
            if rel_error > 0.15 {
                violations_15 += 1;
            }
            if rel_error > 0.50 {
                violations_50 += 1;
            }
            worst_ratio = worst_ratio.max(rel_error);
        }
    }

    let mut clashes = 0usize;
    let mut coincident = 0usize;
    if atom_count_match {
        for i in 0..n {
            for j in (i + 1)..n {
                let a = AtomIdx(i as u32);
                let b = AtomIdx(j as u32);
                let d = coords.get(a).distance(&coords.get(b));
                if d < 1e-3 {
                    coincident += 1;
                }
                // "Gross clash" means non-bonded overlap. Excluding directly-bonded
                // pairs matters: some legitimate bonds (C=O ~1.21A, C#N ~1.16A,
                // C#C ~1.20A) are shorter than CLASH_THRESHOLD_ANGSTROM.  Matches
                // pipeline_v2's own `gross_clash_count` (which does the same
                // exclusion) so both engines are judged by the same definition.
                if mol.bond_between(a, b).is_some() {
                    continue;
                }
                if d < CLASH_THRESHOLD_ANGSTROM {
                    clashes += 1;
                }
            }
        }
    }

    let independently_sound =
        atom_count_match && all_finite && coincident == 0 && worst_ratio <= SOUND_MAX_BOND_RATIO;

    GeometryCheck {
        atom_count_expected: n,
        atom_count_actual,
        atom_count_match,
        all_finite,
        coincident_atom_pairs: coincident,
        worst_bond_length_ratio: worst_ratio,
        bond_violation_rate_15pct: violations_15 as f64 / bond_count as f64,
        bond_violation_rate_50pct: violations_50 as f64 / bond_count as f64,
        gross_clash_count: clashes,
        independently_sound,
    }
}

fn stereo_to_json(v: &chematic_3d::stereo_constraints::StereoVerification) -> Value {
    let tet_satisfied = v.tetrahedral.iter().filter(|r| r.status.is_satisfied()).count();
    let tet_violated = v.tetrahedral.iter().filter(|r| r.status.is_violated()).count();
    let tet_declared = v.tetrahedral.len();
    let db_satisfied = v.double_bond.iter().filter(|r| r.status.is_satisfied()).count();
    let db_violated = v.double_bond.iter().filter(|r| r.status.is_violated()).count();
    let db_declared = v.double_bond.len();
    json!({
        "declared": v.n_declared(),
        "satisfied": v.n_satisfied(),
        "violated": v.n_violations(),
        "unevaluable": v.n_unevaluable(),
        "tetrahedral": {
            "declared": tet_declared,
            "satisfied": tet_satisfied,
            "violated": tet_violated,
            "unevaluable": tet_declared - tet_satisfied - tet_violated,
        },
        "double_bond": {
            "declared": db_declared,
            "satisfied": db_satisfied,
            "violated": db_violated,
            "unevaluable": db_declared - db_satisfied - db_violated,
        },
    })
}

fn load_jsonl(path: &str) -> Vec<Value> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("bad json line in {path}: {e}")))
        .collect()
}

fn load_manifest_smiles(path: &str) -> HashMap<String, String> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    let v: Value =
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("bad json in {path}: {e}"));
    v["molecules"]
        .as_array()
        .unwrap()
        .iter()
        .map(|m| {
            (
                m["name"].as_str().unwrap().to_string(),
                m["smiles"].as_str().unwrap().to_string(),
            )
        })
        .collect()
}

fn score_rows(
    rows: &[Value],
    engine: &str,
    smiles_by_name: &HashMap<String, HashMap<String, String>>,
) {
    for row in rows {
        if row["status"].as_str() != Some("success") {
            continue;
        }
        let tier = row["tier"].as_str().unwrap_or("?").to_string();
        let name = row["name"].as_str().unwrap_or("?").to_string();
        let arm = row["arm"].as_str().unwrap_or("?").to_string();

        let smiles = match smiles_by_name.get(&tier).and_then(|m| m.get(&name)) {
            Some(s) => s.clone(),
            None => {
                println!(
                    "{}",
                    json!({
                        "tier": tier, "name": name, "arm": arm, "engine": engine,
                        "status": "integrity_error", "reason": "smiles_not_found_in_manifest",
                    })
                );
                continue;
            }
        };

        let mol = match chematic_smiles::parse(&smiles) {
            Ok(m) => m,
            Err(e) => {
                println!(
                    "{}",
                    json!({
                        "tier": tier, "name": name, "arm": arm, "engine": engine,
                        "status": "integrity_error", "reason": format!("reparse_failed: {e}"),
                    })
                );
                continue;
            }
        };

        let coords_json = &row["coords"];
        if coords_json.is_null() {
            println!(
                "{}",
                json!({
                    "tier": tier, "name": name, "arm": arm, "engine": engine,
                    "status": "integrity_error", "reason": "success_row_missing_coords",
                })
            );
            continue;
        }

        let coords = match coords_from_json(coords_json, mol.atom_count()) {
            Some(c) => c,
            None => {
                println!(
                    "{}",
                    json!({
                        "tier": tier, "name": name, "arm": arm, "engine": engine,
                        "status": "integrity_error", "reason": "coords_count_mismatch_or_malformed",
                    })
                );
                continue;
            }
        };

        let geom = independent_geometry_check(&mol, &coords);
        let stereo = verify_stereo(&mol, &coords);

        println!(
            "{}",
            json!({
                "tier": tier,
                "name": name,
                "arm": arm,
                "engine": engine,
                "status": "scored",
                "atom_count_expected": geom.atom_count_expected,
                "atom_count_actual": geom.atom_count_actual,
                "atom_count_match": geom.atom_count_match,
                "all_finite": geom.all_finite,
                "coincident_atom_pairs": geom.coincident_atom_pairs,
                "worst_bond_length_ratio": geom.worst_bond_length_ratio,
                "bond_violation_rate_15pct": geom.bond_violation_rate_15pct,
                "bond_violation_rate_50pct": geom.bond_violation_rate_50pct,
                "gross_clash_count": geom.gross_clash_count,
                "independently_sound": geom.independently_sound,
                "stereo": stereo_to_json(&stereo),
            })
        );
    }
}

fn main() {
    let tier_a_smiles =
        load_manifest_smiles("validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json");
    let tier_b_smiles =
        load_manifest_smiles("validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json");
    let mut smiles_by_name: HashMap<String, HashMap<String, String>> = HashMap::new();
    smiles_by_name.insert("A".to_string(), tier_a_smiles);
    smiles_by_name.insert("B".to_string(), tier_b_smiles);

    let chematic_rows = load_jsonl("validation/results/pipeline_v2_vs_rdkit_chematic_rows.jsonl");
    let rdkit_rows = load_jsonl("validation/results/pipeline_v2_vs_rdkit_rdkit_rows.jsonl");

    score_rows(&chematic_rows, "chematic", &smiles_by_name);
    score_rows(&rdkit_rows, "rdkit", &smiles_by_name);
}
