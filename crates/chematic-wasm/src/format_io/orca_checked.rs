use super::*;

/// Write represented ORCA input fields with exact finite Cartesian values.
/// Requires explicit fields, rejects unknown/duplicate JSON fields, unit
/// conflicts, %coords and any metadata that would change on structured reparse.
/// Accepts either kind or type coordinate tags, never both. No job validation.
#[wasm_bindgen]
pub fn write_orca_input_json_checked(input_json: &str) -> Result<String, JsValue> {
    let text = chematic_mol::write_orca_input_json_checked_with_limits(
        input_json,
        &chematic_mol::OrcaWriteLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            max_output_bytes: WASM_MAX_INPUT_BYTES,
            max_atoms: WASM_MAX_ATOMS,
        },
    )
    .map_err(|error| JsValue::from_str(&error.to_string()))?;
    // Keep composition with the public structured WASM reader within its
    // existing byte/physical-line/keyword/block budgets.
    chematic_mol::parse_orca_input_with_limits(&text, &wasm_orca_input_limits())
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    Ok(text)
}
