use crate::{MolHandle, WASM_MAX_OUTPUT_BYTES};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
impl MolHandle {
    /// Checked rooted SMILES as {smiles, atom_order}; order maps parsed to original indices.
    /// Limited to 256 atoms and verified ordinary tetrahedral frames; see the Rust API.
    pub fn smiles_rooted_checked_json(&self, rooted_at_atom: f64) -> Result<String, JsValue> {
        if !rooted_at_atom.is_finite()
            || rooted_at_atom < 0.0
            || rooted_at_atom.fract() != 0.0
            || rooted_at_atom > u32::MAX as f64
        {
            return Err(JsValue::from_str(
                "root atom must be a finite non-negative integer",
            ));
        }
        let (text, order) =
            chematic_smiles::write_rooted_checked(&self.inner, rooted_at_atom as usize)
                .map_err(|e| JsValue::from_str(&e.to_string()))?;
        let order: Vec<_> = order.into_iter().map(|i| i.0).collect();
        let output = serde_json::json!({"smiles":text, "atom_order":order}).to_string();
        if output.len() > WASM_MAX_OUTPUT_BYTES {
            return Err(JsValue::from_str("rooted SMILES exceeds output limit"));
        }
        Ok(output)
    }
}
