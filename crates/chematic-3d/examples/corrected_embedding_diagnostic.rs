//! Emit corrected ETKDG starting coordinates for paired minimizer diagnostics.

use std::io::{self, BufRead};

use chematic_3d::rdkit_embed::{RdkitEmbedOptions, etkdg_embed_multiple_confs_corrected};
use serde_json::{Value, json};

fn main() {
    for line in io::stdin().lock().lines() {
        let line = line.expect("read input row");
        if line.is_empty() {
            continue;
        }
        let row: Value = serde_json::from_str(&line).expect("JSON input row");
        let smiles = row["smiles"].as_str().expect("smiles");
        let molecule =
            chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).expect("SMILES"));
        let mut options = RdkitEmbedOptions::default();
        if let Some(enabled) = row["macrocycle_torsions"].as_bool() {
            options.use_macrocycle_torsions = enabled;
        }
        if let Some(enabled) = row["macrocycle_14"].as_bool() {
            options.use_macrocycle14config = enabled;
        }
        let conformers = etkdg_embed_multiple_confs_corrected(&molecule, 20, &options)
            .expect("corrected embedding");
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
        println!(
            "{}",
            json!({"atomic_numbers": atoms, "bonds": bonds, "coords": conformers})
        );
    }
}
