//! All-frame Classic HISTORY records, bounded for browser use.
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn dlpoly_history_json(text: &str) -> Result<String, JsValue> {
    let g = chematic_mol::parse_dlpoly_history_with_limits(
        text,
        &chematic_mol::DlpolyHistoryLimits {
            max_input_bytes: crate::WASM_MAX_INPUT_BYTES,
            max_atoms: crate::WASM_MAX_ATOMS,
            max_frames: 1000,
            max_atom_records: 100_000,
        },
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let json = serde_json::to_string(&g).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if json.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("HISTORY JSON output limit"));
    }
    Ok(json)
}
#[wasm_bindgen]
pub fn dlpoly_history_from_json(json: &str) -> Result<String, JsValue> {
    if json.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("HISTORY JSON input limit"));
    }
    let g = chematic_mol::dlpoly_history_from_json(json)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if g.atom_count > crate::WASM_MAX_ATOMS
        || g.frames.len() > 1000
        || g.frames
            .len()
            .checked_mul(g.atom_count)
            .is_none_or(|n| n > 100_000)
    {
        return Err(JsValue::from_str("HISTORY frame/atom/total record limit"));
    }
    let out =
        chematic_mol::write_dlpoly_history(&g).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if out.len() > crate::WASM_MAX_INPUT_BYTES {
        return Err(JsValue::from_str(
            "HISTORY text output exceeds reader byte limit",
        ));
    }
    Ok(out)
}
