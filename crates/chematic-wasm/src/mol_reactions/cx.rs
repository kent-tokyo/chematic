use super::*;

/// Enhanced-atom-stereo reaction CXSMILES with group indices remapped after
/// component/atom ordering (RDKit #9732). Other CX fields are rejected.
#[wasm_bindgen]
pub fn normalize_reaction_cxsmiles(text: &str, canonical: bool) -> Result<String, JsValue> {
    let reaction = chematic_rxn::parse_reaction_cxsmiles_with_limits(
        text,
        &chematic_rxn::ReactionParseLimits {
            max_input_bytes: WASM_MAX_INPUT_BYTES,
            max_components_per_side: WASM_MAX_BATCH_ITEMS,
            max_atoms_per_molecule: WASM_MAX_ATOMS,
            max_bonds_per_molecule: WASM_MAX_ATOMS.saturating_mul(2),
        },
    )
    .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let output = chematic_rxn::write_reaction_cxsmiles(&reaction, canonical)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    bounded_json_output(output)
}
