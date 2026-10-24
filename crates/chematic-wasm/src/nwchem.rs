//! Cartesian NWChem GEOMETRY blocks, bounded for browser use.
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn nwchem_geometry_json(text: &str) -> Result<String, JsValue> {
    let record = chematic_mol::parse_nwchem_geometry_with_limits(
        text,
        &chematic_mol::NwchemLimits {
            max_input_bytes: crate::WASM_MAX_INPUT_BYTES,
            max_centers: crate::WASM_MAX_ATOMS,
            ..Default::default()
        },
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let json = serde_json::to_string(&record).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if json.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(
            "NWChem geometry JSON exceeds output limit",
        ));
    }
    Ok(json)
}

#[wasm_bindgen]
pub fn nwchem_geometry_from_json(json: &str) -> Result<String, JsValue> {
    if json.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(
            "NWChem geometry JSON exceeds input limit",
        ));
    }
    let record = chematic_mol::nwchem_geometry_from_json(json)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if record.centers.len() > crate::WASM_MAX_ATOMS {
        return Err(JsValue::from_str(
            "NWChem geometry atom count exceeds limit",
        ));
    }
    let out = chematic_mol::write_nwchem_geometry(&record)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if out.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("NWChem geometry output exceeds limit"));
    }
    Ok(out)
}
