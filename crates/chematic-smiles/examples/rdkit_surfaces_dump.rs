//! Dumps RDKit-port SMILES surfaces for a SMILES file, one JSON object per
//! line: `MolToCXSmiles`, `MolToRandomSmilesVect(5, 42)` and
//! `MolFragmentToSmiles` of the first half of the atoms.
//!
//!     cargo run --release -p chematic-smiles --example rdkit_surfaces_dump -- FILE.smi

use chematic_smiles::{
    RdkitSmilesParams, parse, rdkit_cx_smiles, rdkit_fragment_smiles, rdkit_random_smiles,
};
use serde_json::{Value, json};

fn wrap<T: Into<Value>>(r: Result<T, chematic_smiles::RdkitSmilesError>) -> Value {
    match r {
        Ok(v) => json!({"ok": v.into()}),
        Err(e) => json!({"err": e.to_string()}),
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("SMILES file");
    let text = std::fs::read_to_string(path).expect("read");
    let p = RdkitSmilesParams::default();
    let noncanon = RdkitSmilesParams {
        canonical: false,
        ..p
    };
    let kek = RdkitSmilesParams {
        kekule: true,
        all_hs_explicit: true,
        ..p
    };
    for line in text.lines() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() {
            continue;
        }
        let prob = match chematic_smiles::parse_template(smi) {
            Ok(m) => wrap(
                chematic_smiles::rdkit_detect_chemistry_problems(&m).map(|v| {
                    json!(
                        v.into_iter()
                            .map(|p| json!([p.kind, p.atoms]))
                            .collect::<Vec<_>>()
                    )
                }),
            ),
            Err(e) => json!({"err": e.to_string()}),
        };
        let Ok(mol) = parse(smi) else {
            println!("{}", json!({"smiles": smi, "parse": false, "prob": prob}));
            continue;
        };
        let n = chematic_smiles::rdkit_num_atoms(&mol).unwrap_or(0);
        let atoms: Vec<usize> = (0..(n / 2).max(1)).collect();
        println!(
            "{}",
            json!({
                "smiles": smi,
                "prob": prob,
                "probs": wrap(chematic_smiles::rdkit_detect_chemistry_problems_sanitized(&mol).map(|v| {
                    json!(v.into_iter().map(|p| json!([p.kind, p.atoms])).collect::<Vec<_>>())
                })),
                "cx": wrap(rdkit_cx_smiles(&mol, &p)),
                "rand": wrap(rdkit_random_smiles(&mol, 5, 42, &p)),
                "frag": wrap(rdkit_fragment_smiles(&mol, &atoms, None, &p)),
                "frag2": wrap(rdkit_fragment_smiles(&mol, &(0..n).step_by(2).collect::<Vec<_>>(), None, &p)),
                "frag3": wrap(rdkit_fragment_smiles(&mol, &(n / 3..n).collect::<Vec<_>>(), None, &noncanon)),
                "randk": wrap(rdkit_random_smiles(&mol, 3, 7, &kek)),
                "dm": wrap(chematic_smiles::rdkit_distance_matrix(&mol, false, false).map(|v| json!(v))),
                "dmbo": wrap(chematic_smiles::rdkit_distance_matrix(&mol, true, true).map(|v| json!(v))),
                "dm3d": wrap(chematic_smiles::rdkit_distance_matrix_3d(&mol, &fake_coords(n), true).map(|v| json!(v))),
                "msc": json!({"ok": pairs(chematic_fp::rdkit_morgan_fingerprint(&mol, &Default::default())
                    .map(|r| { let mut v: Vec<(u64, u32)> = r.sparse_counts.into_iter().map(|(k, c)| (u64::from(k), c)).collect(); v.sort_unstable(); v })
                    .unwrap_or_default())}),
                "apsc": json!({"ok": pairs(chematic_fp::rdkit_atom_pair_sparse_counts(&mol).into_iter().map(|(k, c)| (u64::from(k), c)).collect())}),
                "ttsc": json!({"ok": pairs(chematic_fp::rdkit_torsion_sparse_counts(&mol))}),
                "ttleg": json!({"ok": chematic_fp::rdkit_legacy_torsion_counts(&mol).into_iter().map(|(k, c)| json!([k, c])).collect::<Vec<_>>()}),
            })
        );
    }
}

fn pairs(v: Vec<(u64, u32)>) -> Value {
    v.into_iter().map(|(k, c)| json!([k, c])).collect()
}

/// Deterministic test coordinates (the reference dump uses the same).
fn fake_coords(n: usize) -> Vec<[f64; 3]> {
    (0..n)
        .map(|i| {
            [
                i as f64 * 0.5,
                ((i * i) % 7) as f64 * 0.3,
                (i % 3) as f64 * 1.1,
            ]
        })
        .collect()
}
