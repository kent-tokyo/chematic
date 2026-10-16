use crate::{MolHandle, WASM_MAX_ATOMS, WASM_MAX_INPUT_BYTES};
use wasm_bindgen::prelude::*;

fn read(
    input: &str,
    allow_extended_targets: bool,
) -> Result<chematic_mol::PdbConectFrame, JsValue> {
    let mut limits = chematic_mol::PdbConectLimits::default();
    limits.atoms.max_input_bytes = WASM_MAX_INPUT_BYTES;
    limits.atoms.max_atoms = WASM_MAX_ATOMS;
    chematic_mol::parse_pdb_conect_with_options_and_limits(
        input,
        &chematic_mol::PdbConectOptions {
            allow_extended_targets,
        },
        &limits,
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))
}
/// Explicit CONECT graph, without distance-based bonding. Only unmodelled,
/// single-location atoms with explicit element fields are accepted.
#[wasm_bindgen]
pub fn mol_from_pdb_conect(
    input: &str,
    allow_extended_targets: bool,
) -> Result<MolHandle, JsValue> {
    Ok(MolHandle {
        inner: std::rc::Rc::new(read(input, allow_extended_targets)?.molecule),
    })
}
/// Indexed atoms/charges, bonds, coordinates and source serials for the same
/// explicit CONECT projection. Molecular identity and coordinate order agree.
#[wasm_bindgen]
pub fn pdb_conect_frame_json(input: &str, allow_extended_targets: bool) -> Result<String, JsValue> {
    let f = read(input, allow_extended_targets)?;
    let atoms: Vec<_> = f
        .molecule
        .atoms()
        .map(|(_, a)| {
            serde_json::json!({"atomic_number": a.element.atomic_number(), "charge": a.charge})
        })
        .collect();
    let bonds: Vec<_> = f
        .molecule
        .bonds()
        .map(|(_, b)| serde_json::json!([b.atom1.0, b.atom2.0, b.order.order_int()]))
        .collect();
    let coords: Vec<_> = f.coords.points.iter().map(|p| [p.x, p.y, p.z]).collect();
    serde_json::to_string(&serde_json::json!({
        "atoms": atoms, "bonds": bonds, "coords": coords, "serials": f.serials
    }))
    .map_err(|e| JsValue::from_str(&e.to_string()))
}
