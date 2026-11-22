//! Diagnostic only: separate corrected embedding from MMFF stereo screening.
use chematic_core::{Coords3D, Molecule, Point3};
use chematic_3d::rdkit_embed::{RdkitEmbedOptions, etkdg_embed_multiple_confs_corrected};
use chematic_3d::stereo_constraints::verify_stereo;
use chematic_ff::rdkit_mmff::RdkitMmffField;
use serde_json::{json, Value};

fn stereo(mol: &Molecule, coords: &[[f64; 3]]) -> Value {
    let coords = Coords3D { points: coords.iter().map(|&[x,y,z]| Point3::new(x,y,z)).collect() };
    let report = verify_stereo(mol, &coords);
    json!({"declared": report.n_declared(), "violations": report.n_violations(),
           "unevaluable": report.n_unevaluable(), "details": format!("{report:?}")})
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let smiles = std::fs::read_to_string(std::env::args().nth(1).ok_or("missing SMILES file")?)?;
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles.trim())?);
    let candidates = etkdg_embed_multiple_confs_corrected(&mol, 20, &RdkitEmbedOptions::default())?;
    let mut rows = Vec::new();
    if let Some(first) = candidates.first() {
        let field = RdkitMmffField::new(&mol, first, 100.0, true)?;
        for (index, original) in candidates.iter().enumerate() {
            let mut optimized = original.clone();
            let result = field.optimize_detailed(optimized.as_flattened_mut(), 200);
            rows.push(json!({"index": index, "raw": stereo(&mol, original),
                "optimized": stereo(&mol, &optimized), "status": result.status,
                "iterations": result.iterations, "energy": result.energy}));
        }
    }
    println!("{}", serde_json::to_string_pretty(&json!({"atom_count": mol.atom_count(),
        "raw_corrected_count": candidates.len(), "candidates": rows}))?);
    Ok(())
}
