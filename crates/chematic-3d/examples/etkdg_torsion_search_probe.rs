//! Compare the opt-in torsion search with its unchanged quality-profile base.

use std::io::{self, BufRead};

use chematic_3d::rdkit_embed::{
    RdkitEmbedOptions, RdkitMmffConformer, etkdg_embed_multiple_confs_mmff_best_gradient,
    etkdg_embed_multiple_confs_mmff_torsion_search,
};
use serde_json::{Value, json};

fn best(rows: &[RdkitMmffConformer]) -> Value {
    let best = rows
        .iter()
        .filter(|row| row.status == 0 && row.energy.is_finite())
        .min_by(|left, right| left.energy.total_cmp(&right.energy));
    json!({
        "returned": rows.len(),
        "converged": rows.iter().filter(|row| row.status == 0).count(),
        "best_converged_energy": best.map(|row| row.energy),
        "best_converged_coords": best.map(|row| &row.coords),
    })
}

fn main() {
    let options = RdkitEmbedOptions::default();
    for line in io::stdin().lock().lines() {
        let line = line.expect("read input row");
        if line.is_empty() {
            continue;
        }
        let row: Value = serde_json::from_str(&line).expect("JSON input row");
        let name = row["name"].as_str().expect("name");
        let smiles = row["smiles"].as_str().expect("smiles");
        let num_confs = row["num_confs"].as_u64().unwrap_or(20) as usize;
        let max_iters = row["max_iters"].as_u64().unwrap_or(200) as u32;
        let molecule =
            chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).expect("SMILES"));
        let atomic_numbers: Vec<_> = molecule
            .atoms()
            .map(|(_, atom)| atom.element.atomic_number())
            .collect();
        let bonds: Vec<_> = molecule
            .bonds()
            .map(|(_, bond)| {
                json!([
                    bond.atom1.0.min(bond.atom2.0),
                    bond.atom1.0.max(bond.atom2.0),
                    format!("{:?}", bond.order).to_uppercase()
                ])
            })
            .collect();
        let original = etkdg_embed_multiple_confs_mmff_best_gradient(
            &molecule, num_confs, &options, max_iters,
        )
        .expect("quality ensemble");
        let searched = etkdg_embed_multiple_confs_mmff_torsion_search(
            &molecule, num_confs, &options, max_iters,
        )
        .expect("torsion search ensemble");
        println!(
            "{}",
            json!({"name":name,"smiles":smiles,"atomic_numbers":atomic_numbers,
                "bonds":bonds,"num_confs":num_confs,
                "max_iters":max_iters,"quality":best(&original),"torsion":best(&searched)})
        );
    }
}
