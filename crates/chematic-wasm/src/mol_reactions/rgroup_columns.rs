//! Label-stable, column-oriented R-group decomposition for WASM/Node.

use crate::{
    WASM_MAX_ATOMS, enforce_wasm_input_len, enforce_wasm_molecule_size, parse_smiles_json_array,
};
use std::collections::{BTreeMap, BTreeSet};
use wasm_bindgen::prelude::*;

pub(crate) fn rgroup_error_code(error: &chematic_chem::RGroupError) -> &'static str {
    match error {
        chematic_chem::RGroupError::InvalidSmarts(_) => "rgroup_invalid_smarts",
        chematic_chem::RGroupError::Canonicalization(_) => "rgroup_canonicalization",
        chematic_chem::RGroupError::InvalidLabel(_) => "rgroup_invalid_label",
        chematic_chem::RGroupError::DuplicateLabel(_) => "rgroup_duplicate_label",
        chematic_chem::RGroupError::NonTerminalAttachmentLabel(_) => {
            "rgroup_nonterminal_attachment"
        }
        chematic_chem::RGroupError::MultipleAttachments { .. } => "rgroup_multiple_attachments",
        chematic_chem::RGroupError::MatchLimitExceeded { .. } => "rgroup_match_limit",
    }
}

fn rgroup_error_js(error: chematic_chem::RGroupError) -> JsValue {
    let code = rgroup_error_code(&error);
    JsValue::from_str(&format!("{code}: {error}"))
}

/// Label-stable R-group decomposition returned as RDKit-style table columns.
///
/// `core_smarts` accepts the bounded mapped-core contract used by Rust and
/// Python: terminal mapped wildcards such as `[*:2]`, or mapped core atoms.
/// The returned JSON object contains `Core`, `R1`, `R2`, ... arrays parallel
/// to the input rows. Unmatched rows and absent substituents are `null`.
/// Output fragments retain their mapped attachment wildcard. The legacy
/// [`super::rgroup_decompose_json`] function is unchanged.
#[wasm_bindgen]
pub fn rgroup_decompose_columns_json(
    smiles_json: &str,
    core_smarts: &str,
) -> Result<String, JsValue> {
    enforce_wasm_input_len("core_smarts", core_smarts)?;
    let query = chematic_smarts::parse_smarts(core_smarts)
        .map_err(|error| JsValue::from_str(&format!("rgroup_invalid_smarts: {error}")))?;
    if query.atoms.len() > WASM_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "rgroup_query_too_large: core_smarts exceeds maximum query atom count ({WASM_MAX_ATOMS})"
        )));
    }

    let smiles_list = parse_smiles_json_array(smiles_json)?;
    let mols: Vec<chematic_core::Molecule> = smiles_list
        .iter()
        .map(|smiles| {
            let molecule = chematic_smiles::parse(smiles)
                .map_err(|error| JsValue::from_str(&format!("rgroup_invalid_smiles: {error}")))?;
            enforce_wasm_molecule_size(&molecule)?;
            Ok::<_, JsValue>(molecule)
        })
        .collect::<Result<_, _>>()?;
    let refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    let rows =
        chematic_chem::rgroup_decompose_labeled(core_smarts, &refs).map_err(rgroup_error_js)?;
    let labels: BTreeSet<u16> = rows
        .iter()
        .flatten()
        .flat_map(|row| row.r_groups.keys().copied())
        .collect();

    let mut columns: BTreeMap<String, Vec<Option<String>>> = BTreeMap::new();
    columns.insert(
        "Core".to_string(),
        rows.iter()
            .map(|row| row.as_ref().map(|row| row.core_smiles.clone()))
            .collect(),
    );
    for label in labels {
        columns.insert(
            format!("R{label}"),
            rows.iter()
                .map(|row| {
                    row.as_ref()
                        .and_then(|row| row.r_groups.get(&label).cloned())
                })
                .collect(),
        );
    }
    serde_json::to_string(&columns).map_err(|error| JsValue::from_str(&error.to_string()))
}
