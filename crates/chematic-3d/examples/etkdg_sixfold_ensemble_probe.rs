//! Exposed development probe for molecules with a nonzero sixth-order ETKDG
//! torsion coefficient. Compare with
//! `scripts/benchmarks/etkdg_sixfold_ensemble_probe.py`.

use chematic_3d::rdkit_embed::{
    RdkitEmbedOptions, etkdg_embed_multiple_confs_mmff_corrected, rdkit_embed_multiple_confs_mmff,
};
use serde_json::json;

const SMILES: &[&str] = &[
    "CS(=O)(=O)NCC",
    "CCS(=O)(=O)NCC",
    "CS(=O)(=O)NCCO",
    "CS(=O)(=O)NCCc1ccccc1",
    "CS(=O)(=O)NCCC",
    "CS(=O)(=O)NCCF",
];

fn main() {
    let options = RdkitEmbedOptions::default();
    for &smiles in SMILES {
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
        let legacy = rdkit_embed_multiple_confs_mmff(&mol, 20, &options, 200).unwrap();
        let corrected = etkdg_embed_multiple_confs_mmff_corrected(&mol, 20, &options, 200).unwrap();
        println!(
            "{}",
            json!({
                "smiles": smiles,
                "legacy": legacy.iter().map(|c| json!({"status": c.status, "energy": c.energy, "coords": c.coords})).collect::<Vec<_>>(),
                "corrected": corrected.iter().map(|c| json!({"status": c.status, "energy": c.energy, "coords": c.coords})).collect::<Vec<_>>(),
            })
        );
    }
}
