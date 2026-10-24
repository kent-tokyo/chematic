use super::*;
fn limits() -> chematic_mol::GprParseLimits {
    chematic_mol::GprParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms: WASM_MAX_ATOMS,
        max_bonds: 20_000,
        ..Default::default()
    }
}
fn failure(e: impl core::fmt::Display) -> JsValue {
    JsValue::from_str(&e.to_string())
}
/// Read gpr 100 / one coordinate set as editable JSON. Coordinates stay in nm;
/// partial charges and all decimal tokens retain source representation.
#[wasm_bindgen]
pub fn parse_gpr_json(text: &str) -> Result<String, JsValue> {
    let record = chematic_mol::parse_gpr_with_limits(text, &limits()).map_err(failure)?;
    let json = serde_json::to_string(&record).map_err(failure)?;
    if json.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("GPR JSON output limit exceeded"));
    }
    Ok(json)
}
/// Write current fields. Duplicate/unknown keys and unsupported semantics refuse.
#[wasm_bindgen]
pub fn write_gpr_from_json(json: &str) -> Result<String, JsValue> {
    let record = chematic_mol::parse_gpr_json(json, &limits()).map_err(failure)?;
    chematic_mol::write_gpr_with_limits(&record, &limits()).map_err(failure)
}
