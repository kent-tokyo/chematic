//! Evaluate exposed development-set boron centers; run from the workspace root.
//! Output includes parse failures, per-center angles, bond lengths, and stereo.

use chematic_3d::{
    generate_coords, generate_coords_etkdg, get_bond_angle_deg, stereo_constraints::verify_stereo,
};
use serde_json::{Value, json};
fn main() {
    let input: Value = serde_json::from_slice(
        &std::fs::read(
            "validation/results/2026-11-27-selenium-crystal-stereo/source-development-input.json",
        )
        .unwrap(),
    )
    .unwrap();
    let mut out = vec![];
    let mut parse_failures = vec![];
    for row in input["results"].as_array().unwrap() {
        let smiles = row["smiles"].as_str().unwrap();
        let mol = match chematic_smiles::parse(smiles) {
            Ok(m) => m,
            Err(e) => {
                parse_failures.push(json!({"index":row["index"], "error":e.to_string()}));
                continue;
            }
        };
        let centers: Vec<_> = mol
            .atoms()
            .filter_map(|(idx, a)| {
                (a.element.atomic_number() == 5
                    && a.charge == 0
                    && mol.degree(idx) == 3
                    && mol.implicit_hydrogen_count(idx) == 0)
                    .then_some(idx)
            })
            .collect();
        if centers.is_empty() {
            continue;
        }
        let mut stages = vec![];
        for coords in [generate_coords(&mol), generate_coords_etkdg(&mol)] {
            let angles: Vec<_> = centers
                .iter()
                .map(|&c| {
                    let n: Vec<_> = mol.neighbors(c).map(|(x, _)| x).collect();
                    vec![
                        get_bond_angle_deg(&coords, n[0], c, n[1]),
                        get_bond_angle_deg(&coords, n[0], c, n[2]),
                        get_bond_angle_deg(&coords, n[1], c, n[2]),
                    ]
                })
                .collect();
            let stereo = verify_stereo(&mol, &coords);
            let tetra: Vec<_> = stereo
                .tetrahedral
                .iter()
                .map(|r| format!("{:?}", r.status))
                .collect();
            let double: Vec<_> = stereo
                .double_bond
                .iter()
                .map(|r| format!("{:?}", r.status))
                .collect();
            let mut clashes = 0;
            let mut bond_lengths = Vec::new();
            for i in 0..mol.atom_count() {
                for j in i + 1..mol.atom_count() {
                    let a = chematic_core::AtomIdx(i as u32);
                    let b = chematic_core::AtomIdx(j as u32);
                    let distance = coords.get(a).distance(&coords.get(b));
                    if mol.bond_between(a, b).is_some() {
                        bond_lengths.push(distance);
                    } else if distance < 0.8 {
                        clashes += 1;
                    }
                }
            }
            stages.push(json!({"angles":angles,"bond_lengths":bond_lengths,"tetra":tetra,"double":double,"clashes":clashes}));
        }
        out.push(json!({"index":row["index"],"centers":centers.len(),"stages":stages}));
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({"parse_failures":parse_failures,"results":out}))
            .unwrap()
    );
}
