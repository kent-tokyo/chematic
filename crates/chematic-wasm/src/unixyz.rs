use super::*;
fn limits() -> chematic_mol::UnixyzParseLimits {
    chematic_mol::UnixyzParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms: WASM_MAX_ATOMS,
        ..Default::default()
    }
}
fn failure(e: impl core::fmt::Display) -> JsValue {
    JsValue::from_str(&e.to_string())
}

/// Read a bounded UniChem XYZ record as editable JSON. Coordinate fields are decimal strings
/// retaining precision and signed zero in Angstrom.
#[wasm_bindgen]
pub fn parse_unixyz_json(text: &str) -> Result<String, JsValue> {
    let record = chematic_mol::parse_unixyz_with_limits(text, &limits()).map_err(failure)?;
    let json = serde_json::to_string(&record).map_err(failure)?;
    if json.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("unixyz JSON output limit exceeded"));
    }
    Ok(json)
}
/// Serialize validated edited fields; duplicate/unknown JSON keys refuse.
#[wasm_bindgen]
pub fn write_unixyz_from_json(json: &str) -> Result<String, JsValue> {
    let record = chematic_mol::parse_unixyz_json(json, &limits()).map_err(failure)?;
    chematic_mol::write_unixyz_with_limits(&record, &limits()).map_err(failure)
}
