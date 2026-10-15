use super::*;

fn record(text: &str) -> Result<chematic_mol::HinRecord, JsValue> {
    chematic_mol::parse_hin_record_with_limits(
        text,
        chematic_mol::HinParseLimits {
            max_bytes: WASM_MAX_INPUT_BYTES,
            ..Default::default()
        },
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Inspect the bounded HIN record, including literal aromatic bonds and charges.
#[wasm_bindgen]
pub fn parse_hin_record_json(text: &str) -> Result<String, JsValue> {
    let json = record(text)?.to_json().to_string();
    if json.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(
            "HIN record JSON exceeds maximum output size",
        ));
    }
    Ok(json)
}

/// Validate a single HIN record and retain exact source text and line endings.
#[wasm_bindgen]
pub fn roundtrip_hin_record(text: &str) -> Result<String, JsValue> {
    record(text).map(|r| chematic_mol::write_hin_record(&r))
}
