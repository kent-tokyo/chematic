use crate::{WASM_MAX_ATOMS, WASM_MAX_INPUT_BYTES, WASM_MAX_OUTPUT_BYTES};
use wasm_bindgen::prelude::*;

/// Read PDBQT atom data as JSON {atomic_numbers, coords, partial_charges}.
/// Arrays share source ATOM/HETATM order. Coordinates are in Å and partial
/// charges in e; parsed finite binary64 values are not rounded/recalculated.
/// No bonds or formal charges are inferred. Source atom identities, AutoDock
/// types, remarks and torsion trees are not retained. Throws on invalid atoms
/// or WASM size limits. Use mol_from_pdbqt separately for an element-only handle.
#[wasm_bindgen]
pub fn pdbqt_atom_data_json(pdbqt: &str) -> Result<String, JsValue> {
    let limits = chematic_mol::PdbqtParseLimits {
        max_input_bytes: WASM_MAX_INPUT_BYTES,
        max_atoms: WASM_MAX_ATOMS,
        ..Default::default()
    };
    let (mol, coords, charges) = chematic_mol::parse_pdbqt_with_limits(pdbqt, &limits)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let numbers: Vec<_> = mol
        .atoms()
        .map(|(_, atom)| atom.element.atomic_number())
        .collect();
    let text = serde_json::json!({
        "atomic_numbers": numbers,
        "coords": coords,
        "partial_charges": charges,
    })
    .to_string();
    if text.len() > WASM_MAX_OUTPUT_BYTES {
        return Err(JsValue::from_str(
            "PDBQT atom data exceeds WASM output limit",
        ));
    }
    Ok(text)
}
