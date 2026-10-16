use crate::MolHandle;
use wasm_bindgen::prelude::*;

/// One proof result. The witness maps original atoms to a copy with only the
/// selected tetrahedral declaration inverted; it fixes that selected atom.
#[wasm_bindgen]
pub struct TetrahedralCleanupResult {
    result: chematic_chem::TetrahedralStereoCleanup,
}
#[wasm_bindgen]
impl TetrahedralCleanupResult {
    #[wasm_bindgen(getter)]
    pub fn removed(&self) -> bool {
        self.result.removed
    }
    #[wasm_bindgen(getter)]
    pub fn atom_index(&self) -> u32 {
        self.result.atom.0
    }
    #[wasm_bindgen(getter)]
    pub fn search_visits(&self) -> u32 {
        self.result.search_visits as u32
    }
    /// Independent molecule; all source labels and other configurations remain.
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> MolHandle {
        MolHandle {
            inner: std::rc::Rc::new(self.result.molecule.clone()),
        }
    }
    pub fn report_json(&self) -> String {
        serde_json::json!({"removed":self.result.removed,"atom_index":self.result.atom.0,
            "search_visits":self.result.search_visits,
            "witness":self.result.witness.as_ref().map(|v|v.iter().map(|a|a.0).collect::<Vec<_>>())}).to_string()
    }
}
fn integer(value: &JsValue) -> Result<u32, JsValue> {
    value
        .as_f64()
        .filter(|v| v.is_finite() && v.fract() == 0.0 && (0.0..=u32::MAX as f64).contains(v))
        .map(|v| v as u32)
        .ok_or_else(|| {
            JsValue::from_str("cleanup index/budget must be a nonnegative 32-bit integer")
        })
}
/// Opt-in exact tetrahedral cleanup. Both index and budget must be JS numbers
/// containing integers (bool/string/null are rejected). The kernel caps inputs
/// at 64 atoms/128 bonds and max_visits at 10000000. Errors throw without a result.
#[wasm_bindgen]
pub fn remove_redundant_tetrahedral_stereo(
    mol: &MolHandle,
    atom_index: JsValue,
    max_visits: JsValue,
) -> Result<TetrahedralCleanupResult, JsValue> {
    let result = chematic_chem::remove_redundant_tetrahedral_stereo(
        &mol.inner,
        chematic_core::AtomIdx(integer(&atom_index)?),
        u64::from(integer(&max_visits)?),
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(TetrahedralCleanupResult { result })
}
