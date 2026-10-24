use crate::{WASM_MAX_ATOMS, WASM_MAX_INPUT_BYTES, WASM_MAX_JSON_STRING_BYTES};
use wasm_bindgen::prelude::*;

/// Write rigid PDBQT or throw for malformed/unrepresentable input. Coordinates
/// are JSON [x,y,z] arrays in Å; charges are JSON numbers in e. Both arrays must
/// contain exactly one value per graph atom, including explicit H. No defaults,
/// charge calculation, atom removal or ligand torsion preparation are applied.
/// AutoDock types retain the legacy heuristic, with unknown `X` types rejected.
/// Isotopes are refused; smiles_to_pdbqt_isotope_projection records their omission.
#[wasm_bindgen]
pub fn smiles_to_pdbqt_checked(
    smiles: &str,
    coords_json: &str,
    charges_json: &str,
    name: &str,
) -> Result<String, JsValue> {
    let (mol, coords, charges) = pdbqt_inputs(smiles, coords_json, charges_json, name)?;
    chematic_mol::write_pdbqt_checked(&mol, &coords, &charges, name)
        .map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Explicit isotope projection. Returns JSON {pdbqt, isotopes_omitted} with
/// zero-based source_atom_index, one-based pdbqt_serial, atomic_number and
/// mass_number for each label omitted from PDBQT. Retains all explicit atoms,
/// without H merging, charge recalculation or ligand torsion preparation.
#[wasm_bindgen]
pub fn smiles_to_pdbqt_isotope_projection(
    smiles: &str,
    coords_json: &str,
    charges_json: &str,
    name: &str,
) -> Result<String, JsValue> {
    let (mol, coords, charges) = pdbqt_inputs(smiles, coords_json, charges_json, name)?;
    let report = chematic_mol::project_pdbqt_isotopes(&mol, &coords, &charges, name)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    serde_json::to_string(&report).map_err(|e| JsValue::from_str(&e.to_string()))
}

type PdbqtInputs = (chematic_core::Molecule, Vec<(f64, f64, f64)>, Vec<f64>);

fn pdbqt_inputs(
    smiles: &str,
    coords_json: &str,
    charges_json: &str,
    name: &str,
) -> Result<PdbqtInputs, JsValue> {
    let error = |message: String| JsValue::from_str(&message);
    if smiles.len() > WASM_MAX_INPUT_BYTES
        || coords_json.len() > WASM_MAX_JSON_STRING_BYTES
        || charges_json.len() > WASM_MAX_JSON_STRING_BYTES
        || name.len() > 1024
    {
        return Err(error("PDBQT input exceeds WASM size limits".into()));
    }
    let mol = chematic_smiles::parse(smiles).map_err(|e| error(e.to_string()))?;
    if mol.atom_count() > WASM_MAX_ATOMS {
        return Err(error(format!(
            "PDBQT atom count exceeds WASM limit {WASM_MAX_ATOMS}"
        )));
    }
    let coords: Vec<[f64; 3]> =
        serde_json::from_str(coords_json).map_err(|e| error(format!("PDBQT coordinates: {e}")))?;
    let coords: Vec<_> = coords.into_iter().map(|p| (p[0], p[1], p[2])).collect();
    let charges: Vec<f64> =
        serde_json::from_str(charges_json).map_err(|e| error(format!("PDBQT charges: {e}")))?;
    Ok((mol, coords, charges))
}
