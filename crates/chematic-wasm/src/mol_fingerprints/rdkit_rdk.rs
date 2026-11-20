//! Settings and result conversion for the configured RDKit path fingerprint.

use crate::{MolHandle, bounded_json_string, enforce_wasm_input_len, enforce_wasm_molecule_size};
use wasm_bindgen::prelude::*;

#[derive(serde::Deserialize)]
#[serde(default, deny_unknown_fields, rename_all = "camelCase")]
struct RdkOptions {
    min_path: usize,
    max_path: usize,
    nbits: usize,
    bits_per_feature: usize,
    use_hs: bool,
    branched_paths: bool,
    use_bond_order: bool,
    from_atoms: Vec<u32>,
    ignore_atoms: Vec<u32>,
    atom_invariants: Option<Vec<u32>>,
    count_bounds: Option<Vec<u32>>,
    max_visits: usize,
    max_paths: usize,
}

impl Default for RdkOptions {
    fn default() -> Self {
        let default = chematic_fp::RdkitRdkConfig::default();
        Self {
            min_path: default.min_path,
            max_path: default.max_path,
            nbits: default.nbits,
            bits_per_feature: default.bits_per_feature,
            use_hs: default.use_hs,
            branched_paths: default.branched_paths,
            use_bond_order: default.use_bond_order,
            from_atoms: Vec::new(),
            ignore_atoms: Vec::new(),
            atom_invariants: default.atom_invariants,
            count_bounds: default.count_bounds,
            max_visits: default.max_visits,
            max_paths: default.max_paths,
        }
    }
}

/// Packed RDK bits and folded and sparse counts, using input atom indices.
/// A supplied `countBounds` array enables count simulation for the bit vector.
/// Packed bytes are LSB-first and counts are sorted [identifier, count] pairs.
#[wasm_bindgen]
pub fn rdkit_rdk_config_detail_json(
    mol: &MolHandle,
    options_json: &str,
) -> Result<String, JsValue> {
    enforce_wasm_input_len("RDK options", options_json)?;
    enforce_wasm_molecule_size(&mol.inner)?;
    let options: RdkOptions =
        serde_json::from_str(options_json).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let indices = |atoms: Vec<u32>| atoms.into_iter().map(chematic_core::AtomIdx).collect();
    let config = chematic_fp::RdkitRdkConfig {
        min_path: options.min_path,
        max_path: options.max_path,
        nbits: options.nbits,
        bits_per_feature: options.bits_per_feature,
        use_hs: options.use_hs,
        branched_paths: options.branched_paths,
        use_bond_order: options.use_bond_order,
        from_atoms: indices(options.from_atoms),
        ignore_atoms: indices(options.ignore_atoms),
        atom_invariants: options.atom_invariants,
        count_bounds: options.count_bounds,
        max_visits: options.max_visits,
        max_paths: options.max_paths,
    };
    let result = chematic_fp::rdkit_rdk_config(&mol.inner, &config)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    bounded_json_string(&serde_json::json!({
        "fingerprint": result.fingerprint.to_le_bytes(),
        "foldedCounts": result.folded_counts,
        "sparseCounts": result.sparse_counts,
    }))
}
