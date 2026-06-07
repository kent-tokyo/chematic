//! High-level workflow APIs for WASM.
//!
//! Provides single-molecule reports, multi-molecule comparisons, and batch screening
//! as JSON strings suitable for JavaScript/TypeScript integration.

use wasm_bindgen::prelude::*;

/// Generate a complete molecular report (JSON string) from a SMILES.
/// Returns the JSON representation of a `MoleculeReport` struct.
///
/// # Example (JS)
/// ```javascript
/// const json = module.molecule_report_json("CC(=O)Oc1ccccc1C(=O)O");
/// const report = JSON.parse(json);
/// console.log(report.canonical_smiles, report.descriptors.tpsa);
/// ```
#[wasm_bindgen]
pub fn molecule_report_json(smiles: &str) -> Result<String, JsValue> {
    let report = chematic_chem::molecule_report(smiles)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    serde_json::to_string(&report)
        .map_err(|e| JsValue::from_str(&format!("JSON serialization failed: {}", e)))
}

/// Compare two or more SMILES strings (JSON string output).
/// Returns the JSON representation of a `MoleculeComparison` struct.
///
/// # Example (JS)
/// ```javascript
/// const json = module.compare_molecules_json("c1ccccc1", "Cc1ccccc1");
/// const comparison = JSON.parse(json);
/// console.log(comparison.pairwise[0].similarities.ecfp4_tanimoto);
/// ```
#[wasm_bindgen]
pub fn compare_molecules_json(smiles1: &str, smiles2: &str) -> Result<String, JsValue> {
    let smiles = [smiles1, smiles2];
    let comparison = chematic_chem::compare_molecules(&smiles)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    serde_json::to_string(&comparison)
        .map_err(|e| JsValue::from_str(&format!("JSON serialization failed: {}", e)))
}

/// Compare multiple SMILES strings (up to 256 by default).
/// Accepts a delimiter-separated list (e.g., newline or comma).
///
/// # Example (JS)
/// ```javascript
/// const smilesList = "c1ccccc1\nCc1ccccc1\nCCc1ccccc1";
/// const json = module.compare_molecules_batch_json(smilesList, "\n");
/// const comparison = JSON.parse(json);
/// ```
#[wasm_bindgen]
pub fn compare_molecules_batch_json(
    smiles_batch: &str,
    delimiter: &str,
) -> Result<String, JsValue> {
    let smiles_vec: Vec<&str> = smiles_batch.split(delimiter).collect();
    let comparison = chematic_chem::compare_molecules(&smiles_vec)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    serde_json::to_string(&comparison)
        .map_err(|e| JsValue::from_str(&format!("JSON serialization failed: {}", e)))
}

/// Screen a batch of SMILES strings (JSON string output).
/// Returns per-record results including pass/fail with error details.
/// Includes MaxMin diversity picking and Butina clustering by default.
///
/// # Example (JS)
/// ```javascript
/// const smilesList = "c1ccccc1\nCC\nCCC";
/// const json = module.screen_smiles_json(smilesList, "\n");
/// const report = JSON.parse(json);
/// console.log(report.records); // Array of ScreeningRecord
/// console.log(report.maxmin_picks); // Diversity-selected indices
/// console.log(report.butina_clusters); // Clustering result
/// ```
#[wasm_bindgen]
pub fn screen_smiles_json(smiles_batch: &str, delimiter: &str) -> String {
    let smiles_vec: Vec<&str> = smiles_batch.split(delimiter).collect();
    let report = chematic_chem::screen_smiles(&smiles_vec);
    serde_json::to_string(&report)
        .unwrap_or_else(|_| "{\"error\": \"JSON serialization failed\"}".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn workflow_json_serialization_aspirin() {
        // Test that the underlying workflow API produces serializable data
        let report = chematic_chem::molecule_report("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let json = serde_json::to_string(&report).unwrap();
        assert!(json.contains("molecular_weight"), "JSON should contain molecular_weight");
        assert!(json.contains("tpsa"), "JSON should contain tpsa");
    }

    #[test]
    fn workflow_json_serialization_compare() {
        let comparison = chematic_chem::compare_molecules(&["c1ccccc1", "Cc1ccccc1"]).unwrap();
        let json = serde_json::to_string(&comparison).unwrap();
        assert!(json.contains("pairwise"), "JSON should contain pairwise");
        assert!(json.contains("ecfp4_tanimoto"), "JSON should have similarity metrics");
    }

    #[test]
    fn workflow_json_serialization_screen() {
        let report = chematic_chem::screen_smiles(&["c1ccccc1", "CC", "CCC"]);
        let json = serde_json::to_string(&report).unwrap();
        assert!(json.contains("records"), "JSON should contain records");
        assert!(json.contains("maxmin_picks"), "JSON should contain diversity picks");
    }
}
