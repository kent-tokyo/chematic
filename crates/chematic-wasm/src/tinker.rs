//! Editable Tinker records. JSON is the same strict schema as Rust/Python.
use wasm_bindgen::prelude::*;

/// Read one Tinker XYZ frame. Names/types and declared adjacency remain opaque.
#[wasm_bindgen]
pub fn tinker_xyz_record_json(text: &str) -> Result<String, JsValue> {
    let limits = chematic_mol::tinker::TinkerParseLimits {
        max_input_bytes: crate::WASM_MAX_INPUT_BYTES,
        max_atoms: crate::WASM_MAX_ATOMS,
        ..Default::default()
    };
    let record = chematic_mol::tinker::parse_tinker_xyz_with_limits(text, &limits)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let out = serde_json::to_string(&record).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if out.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("Tinker JSON exceeds output limit"));
    }
    Ok(out)
}

/// Write an edited record, refusing unsupported keys rather than dropping them.
#[wasm_bindgen]
pub fn tinker_xyz_from_record_json(json: &str) -> Result<String, JsValue> {
    // A compact text frame can expand past the text-input limit as JSON.
    // Accept the same bound the reader uses for its returned record.
    if json.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("Tinker JSON exceeds input limit"));
    }
    let record = chematic_mol::tinker::tinker_xyz_from_json(json)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if record.atoms.len() > crate::WASM_MAX_ATOMS {
        return Err(JsValue::from_str("Tinker atom count exceeds limit"));
    }
    let out = chematic_mol::tinker::write_tinker_xyz(&record)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if out.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("Tinker output exceeds limit"));
    }
    Ok(out)
}
