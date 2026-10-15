//! Emit one RDKit-compatible SMILES result per indexed V3000 SDF record.
//! Used by validation/cosmolkit_comparison/check_sdf_v3000_accuracy.py.

use std::error::Error;

use chematic_mol::SdfDataset;
use chematic_smiles::rdkit_canonical_smiles;
use serde_json::json;

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: sdf_v3000_accuracy_dump PATH")?;
    let dataset = SdfDataset::open(path)?;
    let mut rows = Vec::with_capacity(dataset.len());
    for index in 0..dataset.len() {
        let row = match dataset.get(index) {
            Ok(Some(record)) => match rdkit_canonical_smiles(&record.mol) {
                Ok(smiles) => json!({"smiles": smiles, "atoms": record.mol.atom_count()}),
                Err(error) => json!({"error": error.to_string()}),
            },
            Ok(None) => json!({"error": "record missing"}),
            Err(error) => json!({"error": error.to_string()}),
        };
        rows.push(row);
    }
    println!("{}", serde_json::to_string(&rows)?);
    Ok(())
}
