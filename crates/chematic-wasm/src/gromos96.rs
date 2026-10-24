use super::*;
fn limits() -> chematic_mol::Gromos96ParseLimits {
    chematic_mol::Gromos96ParseLimits {
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
        return Err(JsValue::from_str("GROMOS96 JSON output limit exceeded"));
    }
    Ok(text)
}
/// Read geometry as editable JSON, with all blocks, fixed-width payload and caller-owned length-unit context.
#[wasm_bindgen]
pub fn parse_gromos96_json(text: &str, length_unit: &str) -> Result<String, JsValue> {
    json(
        &chematic_mol::parse_gromos96_with_limits(
            text,
            length_unit.parse().map_err(failure)?,
            &limits(),
        )
        .map_err(failure)?,
    )
}
/// Write current fields without rescaling; unknown or duplicate keys refuse.
#[wasm_bindgen]
pub fn write_gromos96_from_json(text: &str) -> Result<String, JsValue> {
    let r = chematic_mol::parse_gromos96_json(text, &limits()).map_err(failure)?;
    chematic_mol::write_gromos96_with_limits(&r, &limits()).map_err(failure)
}
/// Explicit checked binary64 Angstrom positions; raw parse/write is independent.
#[wasm_bindgen]
pub fn gromos96_geometry_angstrom_json(text: &str) -> Result<String, JsValue> {
    let r = chematic_mol::parse_gromos96_json(text, &limits()).map_err(failure)?;
    json(
        &r.geometry_angstrom_with_limits(&limits())
            .map_err(failure)?,
    )
}
