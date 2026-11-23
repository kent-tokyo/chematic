//! Compare native declared-stereo reports against RDKit-generated geometries.
use chematic_3d::stereo_constraints::verify_stereo;
use chematic_core::{Coords3D, Point3};
use serde_json::{json, Value};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input: Value = serde_json::from_slice(&std::fs::read(std::env::args().nth(1).ok_or("missing input")?)?)?;
    let mut rows = Vec::new();
    for row in input.as_array().ok_or("expected cases")? {
        let smiles = row["smiles"].as_str().ok_or("missing smiles")?;
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles)?);
        let mut conformers = Vec::new();
        for conformer in row["conformers"].as_array().ok_or("missing conformers")? {
            let coords: Vec<[f64; 3]> = serde_json::from_value(conformer["coords"].clone())?;
            assert_eq!(coords.len(), mol.atom_count());
            let geometry = Coords3D { points: coords.iter().map(|&[x,y,z]| Point3::new(x,y,z)).collect() };
            let stereo = verify_stereo(&mol, &geometry);
            conformers.push(json!({"external": conformer["external"], "optimizer_status": conformer["status"],
                "declared": stereo.n_declared(), "violations": stereo.n_violations(),
                "unevaluable": stereo.n_unevaluable(), "details": format!("{stereo:?}")}));
        }
        rows.push(json!({"name": row["name"], "conformers": conformers}));
    }
    println!("{}", serde_json::to_string_pretty(&rows)?);
    Ok(())
}
