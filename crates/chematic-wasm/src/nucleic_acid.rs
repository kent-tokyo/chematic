use wasm_bindgen::prelude::*;

/// Validate and normalize a bounded `chematic.nucleic-acid.v1` document.
///
/// Returns a stable JSON envelope with either `ok: true` and the normalized
/// document or `ok: false` and a typed error category.
#[wasm_bindgen]
pub fn nucleic_acid_validate_json(document_json: &str) -> String {
    serde_json::to_string(&chematic_mol::validate_nucleic_acid_json(
        document_json,
        &chematic_mol::NucleicAcidLimits::default(),
    ))
    .expect("nucleic-acid validation envelope is serializable")
}

/// Apply one bounded metadata edit without changing atom ownership or linkage
/// topology. Returns the same tagged envelope as validation.
#[wasm_bindgen]
pub fn nucleic_acid_apply_json_command(document_json: &str, command_json: &str) -> String {
    serde_json::to_string(&chematic_mol::apply_nucleic_acid_json_command(
        document_json,
        command_json,
        &chematic_mol::NucleicAcidLimits::default(),
    ))
    .expect("nucleic-acid edit envelope is serializable")
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{Value, json};

    #[test]
    fn typed_error_shape_matches_binding_contract() {
        let value: Value = serde_json::from_str(&nucleic_acid_validate_json(
            &json!({"schema":"wrong"}).to_string(),
        ))
        .unwrap();
        assert_eq!(value["ok"], false);
        assert_eq!(value["error"]["code"], "invalid_json");
        assert_eq!(value["schema"], "chematic.nucleic-acid-validation.v1");
    }
}
