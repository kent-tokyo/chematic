//! Source-built MMFF94 ensemble probe for a fixed manifest.

use chematic_3d::rdkit_embed::{
    RdkitEmbedOptions, etkdg_embed_multiple_confs_mmff_best_gradient,
    rdkit_embed_multiple_confs_mmff,
};
use serde_json::{Value, json};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let manifest = std::env::args().nth(1).ok_or("expected manifest path")?;
    let limit: usize = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "24".into())
        .parse()?;
    let start: usize = std::env::args()
        .nth(3)
        .unwrap_or_else(|| "0".into())
        .parse()?;
    let data: Value = serde_json::from_str(&std::fs::read_to_string(manifest)?)?;
    let options = RdkitEmbedOptions::default();
    for row in data["molecules"]
        .as_array()
        .ok_or("missing molecules")?
        .iter()
        .skip(start)
        .take(limit)
    {
        let name = row["name"].as_str().ok_or("missing name")?;
        let smiles = row["smiles"].as_str().ok_or("missing smiles")?;
        let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles)?);
        for (profile, results) in [
            (
                "compat",
                rdkit_embed_multiple_confs_mmff(&molecule, 24, &options, 1000),
            ),
            (
                "quality",
                etkdg_embed_multiple_confs_mmff_best_gradient(&molecule, 24, &options, 1000),
            ),
        ] {
            let results = results?;
            let best = results
                .iter()
                .filter(|r| r.status == 0 && r.energy.is_finite())
                .map(|r| r.energy)
                .min_by(f64::total_cmp);
            let record = json!({"name":name,"profile":profile,"requested":24,
                "retained":results.len(),"converged":results.iter().filter(|r| r.status==0).count(),
                "best_converged_energy":best});
            println!("{record}");
        }
    }
    Ok(())
}
