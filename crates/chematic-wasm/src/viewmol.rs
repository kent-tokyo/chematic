use super::*;
fn limits() -> chematic_mol::ViewmolParseLimits {
    chematic_mol::ViewmolParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms: WASM_MAX_ATOMS,
        ..Default::default()
    }
}
fn failure(e: impl core::fmt::Display) -> JsValue {
    JsValue::from_str(&e.to_string())
}
fn json<T: serde::Serialize>(value: &T) -> Result<String, JsValue> {
    let text = serde_json::to_string(value).map_err(failure)?;
    if text.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("ViewMol JSON output limit exceeded"));
    }
    Ok(text)
}
/// Read geometry as editable JSON, with raw factor/coordinate strings and fixed axes.
#[wasm_bindgen]
pub fn parse_viewmol_json(text: &str) -> Result<String, JsValue> {
    json(&chematic_mol::parse_viewmol_with_limits(text, &limits()).map_err(failure)?)
}
/// Write current fields without rescaling; unknown or duplicate keys refuse.
#[wasm_bindgen]
pub fn write_viewmol_from_json(text: &str) -> Result<String, JsValue> {
    let r = chematic_mol::parse_viewmol_json(text, &limits()).map_err(failure)?;
    chematic_mol::write_viewmol_with_limits(&r, &limits()).map_err(failure)
}
/// Explicit checked binary64 Angstrom positions; raw parse/write is independent.
#[wasm_bindgen]
pub fn viewmol_positions_angstrom_json(text: &str) -> Result<String, JsValue> {
    let r = chematic_mol::parse_viewmol_json(text, &limits()).map_err(failure)?;
    json(
        &r.positions_angstrom_with_limits(&limits())
            .map_err(failure)?,
    )
}
