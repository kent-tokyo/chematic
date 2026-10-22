use super::*;

fn document(text: &str) -> Result<chematic_mol::GjfSourceDocument, JsValue> {
    chematic_mol::parse_gjf_source_document_with_limits(
        text,
        &chematic_mol::GjfDocumentLimits {
            gaussian: chematic_mol::GaussianParseLimits {
                max_input_bytes: WASM_MAX_INPUT_BYTES,
                max_atoms: WASM_MAX_ATOMS,
                ..Default::default()
            },
            ..Default::default()
        },
    )
    .map_err(|error| JsValue::from_str(&error.to_string()))
}

/// Retain all explicit Cartesian Link1 jobs as exact source; no job selection.
/// Budgets are cumulative. Opaque sections and external resources are not interpreted.
#[wasm_bindgen]
pub fn roundtrip_gjf_document(text: &str) -> Result<String, JsValue> {
    document(text).map(|document| chematic_mol::write_gjf_source_document(&document))
}

/// Inspect every ordered Cartesian job with UTF-8 byte ranges in the original source.
/// Coordinate units are unchanged; no inferred bonds or external resource reads.
#[wasm_bindgen]
pub fn parse_gjf_document_json(text: &str) -> Result<String, JsValue> {
    let json = document(text)?.to_json().to_string();
    if json.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(
            "GJF document JSON exceeds maximum output size",
        ));
    }
    Ok(json)
}
