//! DL-POLY records retain source units and opaque force-field names.
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn dlpoly_config_record_json(text: &str) -> Result<String, JsValue> {
    let record = chematic_mol::parse_dlpoly_config_with_limits(
        text,
        &chematic_mol::DlpolyConfigLimits {
            max_input_bytes: crate::WASM_MAX_INPUT_BYTES,
            max_atoms: crate::WASM_MAX_ATOMS,
        },
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let json = serde_json::to_string(&record).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if json.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("CONFIG JSON exceeds output limit"));
    }
    Ok(json)
}

#[wasm_bindgen]
pub fn dlpoly_config_from_record_json(json: &str) -> Result<String, JsValue> {
    if json.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("CONFIG JSON exceeds input limit"));
    }
    let record = chematic_mol::dlpoly_config_from_json(json)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if record.atoms.len() > crate::WASM_MAX_ATOMS {
        return Err(JsValue::from_str("CONFIG atom count exceeds limit"));
    }
    let out = chematic_mol::write_dlpoly_config(&record)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if out.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("CONFIG output exceeds limit"));
    }
    Ok(out)
}
