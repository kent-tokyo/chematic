use crate::{MolHandle, bounded_json_string, enforce_wasm_input_len};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
impl MolHandle {
    /// Checked 3D matrix in prepared RDKit atom order. JSON has one [x,y,z]
    /// row per prepared atom; malformed/nonfinite/unrepresentable input refuses.
    pub fn rdkit_distance_matrix_3d_checked_json(
        &self,
        coords_json: &str,
        use_atom_wts: bool,
    ) -> Result<String, JsValue> {
        enforce_wasm_input_len("coords_json", coords_json)?;
        let coords: Vec<[f64; 3]> = serde_json::from_str(coords_json)
            .map_err(|e| JsValue::from_str(&format!("invalid 3D coordinates JSON: {e}")))?;
        let matrix =
            chematic_smiles::rdkit_distance_matrix_3d_checked(&self.inner, &coords, use_atom_wts)
                .map_err(|e| JsValue::from_str(&e.to_string()))?;
        bounded_json_string(&matrix)
    }
}
