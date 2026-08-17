//! Issue #227 Phase 2 Step 6 (BCI derived-formal-charge fix): surveys the
//! full 264-molecule Wave 1 corpus for every heavy atom whose MMFF94
//! numeric type is one RDKit's real `computeMMFFCharges` special-cases in
//! its "set formal charges upfront" switch (types 32/72 O2CM/SM,
//! 34/49/51/54/58/92/93/94/97 simple+1, 87/95/96/98/99 simple+2, 88
//! simple+3, 35/62/89/90/91 simple-1, 76 N5M, 55/56/81 conjugated-cation,
//! 61 diazonium) -- regardless of the atom's RAW formal charge, since
//! RDKit's derived charge for these types can be nonzero even when the raw
//! charge is zero (e.g. a carboxylate `=O` with raw charge 0 still gets a
//! shared -0.5 under O2CM redistribution). This is the accurate "what could
//! regress or needs corpus-exercised verification" survey a plain
//! `raw_charge != 0` filter would undercount (neighbor derived-charge also
//! feeds into an atom's own computation via `sum_fc`/the anionic leak).
//!
//! Committed (not deleted) so `mmff_derived_formal_charge`'s and
//! `o2cm_sm_formal_charge`'s "zero corpus exposure" claims in their doc
//! comments and in `scripts/mmff94_provenance/PROVENANCE.md` are
//! independently re-runnable, not just asserted -- same precedent as
//! `mmff94_bci_charges_dump_227.rs` and
//! `mmff94_bci_stereo_drift_diagnostic_227.rs` in this same directory.
//! Measurement-only: no `chematic-ff` production code is touched by this
//! file.
//!
//! Run: `cargo run --release -p chematic-3d --example mmff94_fchg_type_exposure_survey_227`

use chematic_core::AtomIdx;
use serde_json::Value;
use std::collections::BTreeSet;

fn load_manifest(path: &str) -> Vec<(String, String)> {
    let text = std::fs::read_to_string(path).unwrap();
    let v: Value = serde_json::from_str(&text).unwrap();
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

fn main() {
    let special: BTreeSet<u8> = [
        32, 72, 34, 49, 51, 54, 58, 92, 93, 94, 97, 87, 95, 96, 98, 99, 88, 35, 62, 89, 90, 91, 76,
        55, 56, 81, 61,
    ]
    .into_iter()
    .collect();

    let manifests: Vec<(&str, &str)> = vec![
        (
            "A",
            "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
        ),
        (
            "B",
            "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
        ),
    ];
    let mut hits: Vec<String> = Vec::new();
    let mut molecules_with_special: BTreeSet<String> = BTreeSet::new();
    let mut per_type_count: std::collections::BTreeMap<u8, usize> = Default::default();
    for (tier, path) in manifests {
        for (name, smiles) in load_manifest(path) {
            let mol = match chematic_smiles::parse(&smiles) {
                Ok(m) => m,
                Err(_) => continue,
            };
            let types = match chematic_ff::assign_mmff94_numeric_types(&mol) {
                Ok(t) => t,
                Err(_) => continue,
            };
            for i in 0..mol.atom_count() {
                let idx_i = AtomIdx(i as u32);
                let a = mol.atom(idx_i);
                if special.contains(&types[i]) {
                    molecules_with_special.insert(format!("{tier}:{name}"));
                    *per_type_count.entry(types[i]).or_insert(0) += 1;
                    let neighbor_desc: Vec<String> = mol
                        .neighbors(idx_i)
                        .map(|(nidx, _)| {
                            let n = mol.atom(nidx);
                            format!(
                                "#{}(type={},raw_q={},elem={:?})",
                                nidx.0, types[nidx.0 as usize], n.charge, n.element
                            )
                        })
                        .collect();
                    hits.push(format!(
                        "{tier}:{name}#{i} elem={:?} type={} raw_charge={} neighbors=[{}]",
                        a.element,
                        types[i],
                        a.charge,
                        neighbor_desc.join(", ")
                    ));
                }
            }
        }
    }
    println!(
        "=== {} atoms with a 'special' fChg-switch type, across {} molecules ===",
        hits.len(),
        molecules_with_special.len()
    );
    for h in &hits {
        println!("  {h}");
    }
    println!("=== per-type count ===");
    for (t, c) in &per_type_count {
        println!("  type {t}: {c} atoms");
    }
    eprintln!(
        "total_special_atoms={} molecules_with_special={} types_present={:?}",
        hits.len(),
        molecules_with_special.len(),
        per_type_count.keys().collect::<Vec<_>>()
    );
}
