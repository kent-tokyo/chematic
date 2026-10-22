use super::*;

fn limits() -> chematic_mol::SdfSourceLimits {
    chematic_mol::SdfSourceLimits {
        sdf: chematic_mol::SdfParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            max_record_bytes: WASM_MAX_INPUT_BYTES,
            ..chematic_mol::SdfSourceLimits::default().sdf
        },
        max_atoms_per_record: WASM_MAX_ATOMS,
        max_total_atoms: WASM_MAX_ATOMS,
    }
}

/// Validate every local-profile SDF record and retain exact UTF-8 source.
/// No partial output on record/limit failure. Opaque metadata/stereo is not
/// interpreted or edited. Bounds: 1 MB input/record, 10000 records/total atoms.
#[wasm_bindgen]
pub fn roundtrip_sdf_source(text: &str) -> Result<String, JsValue> {
    chematic_mol::roundtrip_sdf_source_with_limits(text, &limits())
        .map_err(|error| JsValue::from_str(&error.to_string()))
}

pub(super) fn single_block(text: &str) -> Result<&str, JsValue> {
    chematic_mol::single_sdf_mol_block_with_limits(text, &limits())
        .map_err(|error| JsValue::from_str(&error.to_string()))
}
