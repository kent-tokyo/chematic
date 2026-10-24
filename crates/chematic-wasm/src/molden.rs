//! Molden geometry/vibration JSON, preserving explicit units and mode identity.
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn molden_record_json(text: &str) -> Result<String, JsValue> {
    let record = chematic_mol::parse_molden_with_limits(
        text,
        &chematic_mol::MoldenLimits {
            max_input_bytes: crate::WASM_MAX_INPUT_BYTES,
            max_atoms: crate::WASM_MAX_ATOMS,
            max_modes: 10_000,
            max_displacements: 100_000,
            ..Default::default()
        },
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let json = serde_json::to_string(&record).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if json.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("Molden JSON exceeds output limit"));
    }
    Ok(json)
}

#[wasm_bindgen]
pub fn molden_from_record_json(json: &str) -> Result<String, JsValue> {
    if json.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("Molden JSON exceeds input limit"));
    }
    let record = chematic_mol::molden_record_from_json(json)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    if record
        .atoms
        .as_ref()
        .is_some_and(|v| v.atoms.len() > crate::WASM_MAX_ATOMS)
        || record
            .fr_coords
            .as_ref()
            .is_some_and(|v| v.len() > crate::WASM_MAX_ATOMS)
        || record
            .frequencies
            .as_ref()
            .is_some_and(|v| v.len() > 10_000)
        || record
            .normal_modes
            .as_ref()
            .is_some_and(|m| m.iter().map(|m| m.displacements_bohr.len()).sum::<usize>() > 100_000)
    {
        return Err(JsValue::from_str(
            "Molden atom/mode/displacement count exceeds limit",
        ));
    }
    let out = chematic_mol::write_molden(&record).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if out.len() > crate::WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str("Molden output exceeds limit"));
    }
    Ok(out)
}
