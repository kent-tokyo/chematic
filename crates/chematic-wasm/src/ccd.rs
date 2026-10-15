use crate::{MolHandle, bounded_json_string, enforce_wasm_input_len};
use wasm_bindgen::prelude::*;

fn parse(text: &str) -> Result<chematic_mol::CcdComponent, JsValue> {
    enforce_wasm_input_len("CCD text", text)?;
    chematic_mol::parse_ccd_component(text).map_err(|e| JsValue::from_str(&e.to_string()))
}
/// One CCD template including source IDs, charges, explicit bonds, optional
/// coordinates and unassigned stereo labels. No atom_site bond perception.
#[wasm_bindgen]
pub fn parse_ccd_json(text: &str) -> Result<String, JsValue> {
    bounded_json_string(&parse(text)?)
}
/// Project explicit CCD topology. Unassigned stereo annotations require
/// explicit loss consent; coordinates remain in parse_ccd_json's record.
#[wasm_bindgen]
pub fn mol_from_ccd(text: &str, allow_stereo_loss: bool) -> Result<MolHandle, JsValue> {
    let inner = parse(text)?
        .to_molecule(allow_stereo_loss)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(MolHandle {
        inner: inner.into(),
    })
}
