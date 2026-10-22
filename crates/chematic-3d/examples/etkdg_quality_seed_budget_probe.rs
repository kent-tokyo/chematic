//! Read JSONL `{ "name": ..., "smiles": ... }` rows from stdin and report
//! the compatibility and quality-profile ETKDG/MMFF94 results. The paired
//! RDKit scorer is `scripts/benchmarks/etkdg_quality_seed_budget_probe.py`.

use std::io::{self, BufRead};

use chematic_3d::rdkit_embed::{
    RdkitEmbedOptions, RdkitMmffConformer, etkdg_embed_multiple_confs_mmff_best_gradient,
    rdkit_embed_multiple_confs_mmff,
};
use serde_json::{Value, json};

fn summary(conformers: &[RdkitMmffConformer]) -> Value {
    let best = conformers
        .iter()
        .filter(|conformer| conformer.energy.is_finite())
        .min_by(|a, b| a.energy.total_cmp(&b.energy));
    let best_converged = conformers
        .iter()
        .filter(|conformer| conformer.status == 0 && conformer.energy.is_finite())
        .min_by(|a, b| a.energy.total_cmp(&b.energy));
    json!({
        "generated": conformers.len(),
        "converged": conformers.iter().filter(|conformer| conformer.status == 0).count(),
        "best_energy": best.map(|conformer| conformer.energy),
        "best_coords": best.map(|conformer| &conformer.coords),
        "best_converged_energy": best_converged.map(|conformer| conformer.energy),
        "best_converged_coords": best_converged.map(|conformer| &conformer.coords),
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
        let view = chematic_smiles::rdkit_mol_view(&molecule).expect("RDKit atom view");
        let torsions = chematic_3d::rdkit_embed::etk::get_experimental_torsions(
            &view, true, false, true, true, 2,
        );
        let high_order_gradient_terms = torsions
            .exp_torsion_angles
            .iter()
            .any(|(_, constants)| constants[4] != 0.0 || constants[5] != 0.0);
        let atoms: Vec<_> = molecule
            .atoms()
            .map(|(_, atom)| atom.element.atomic_number())
            .collect();
        let bonds: Vec<_> = molecule
            .bonds()
            .map(|(_, bond)| {
                json!([
                    bond.atom1.0.min(bond.atom2.0),
                    bond.atom1.0.max(bond.atom2.0),
                    format!("{:?}", bond.order).to_uppercase(),
                ])
            })
            .collect();
        let compatibility =
            rdkit_embed_multiple_confs_mmff(&molecule, num_confs, &options, max_iters)
                .expect("compatibility ensemble");
        let quality = etkdg_embed_multiple_confs_mmff_best_gradient(
            &molecule, num_confs, &options, max_iters,
        )
        .expect("quality ensemble");
        println!(
            "{}",
            json!({
                "name": name,
                "smiles": smiles,
                "high_order_gradient_terms": high_order_gradient_terms,
                "atomic_numbers": atoms,
                "bonds": bonds,
                "compatibility": summary(&compatibility),
                "quality": summary(&quality),
                "quality_candidates": quality.iter().map(|conformer| json!({
                    "status": conformer.status,
                    "energy": conformer.energy,
                    "coords": &conformer.coords,
                })).collect::<Vec<_>>(),
            })
        );
    }
}
