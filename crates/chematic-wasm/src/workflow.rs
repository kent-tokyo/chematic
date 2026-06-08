//! High-level workflow APIs for WASM.
//!
//! Provides single-molecule reports, multi-molecule comparisons, and batch screening
//! as JSON strings suitable for JavaScript/TypeScript integration.

use wasm_bindgen::prelude::*;

const WORKFLOW_MAX_INPUT_BYTES: usize = 1_000_000;
const WORKFLOW_MAX_BATCH_ITEMS: usize = 1_024;

fn enforce_input_len(label: &str, input: &str) -> Result<(), JsValue> {
    if input.len() > WORKFLOW_MAX_INPUT_BYTES {
        return Err(JsValue::from_str(&format!(
            "{label} exceeds maximum input size ({} > {} bytes)",
            input.len(),
            WORKFLOW_MAX_INPUT_BYTES
        )));
    }
    Ok(())
}

fn split_bounded_batch<'a>(
    smiles_batch: &'a str,
    delimiter: &str,
    label: &str,
) -> Result<Vec<&'a str>, JsValue> {
    enforce_input_len(label, smiles_batch)?;
    if delimiter.is_empty() {
        return Err(JsValue::from_str("delimiter must not be empty"));
    }
    let smiles_vec: Vec<&str> = smiles_batch.split(delimiter).collect();
    if smiles_vec.len() > WORKFLOW_MAX_BATCH_ITEMS {
        return Err(JsValue::from_str(&format!(
            "{label} exceeds maximum item count ({} > {})",
            smiles_vec.len(),
            WORKFLOW_MAX_BATCH_ITEMS
        )));
    }
    Ok(smiles_vec)
}

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
    enforce_input_len("smiles", smiles)?;
    let report =
        chematic_chem::molecule_report(smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
    enforce_input_len("smiles1", smiles1)?;
    enforce_input_len("smiles2", smiles2)?;
    let smiles = [smiles1, smiles2];
    let comparison =
        chematic_chem::compare_molecules(&smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
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
    let smiles_vec = split_bounded_batch(smiles_batch, delimiter, "smiles_batch")?;
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
    let smiles_vec = match split_bounded_batch(smiles_batch, delimiter, "smiles_batch") {
        Ok(values) => values,
        Err(e) => {
            let msg = e
                .as_string()
                .unwrap_or_else(|| "invalid batch input".to_string());
            return format!("{{\"error\":\"{}\"}}", msg.replace('"', "\\\""));
        }
    };
    let report = chematic_chem::screen_smiles(&smiles_vec);
    serde_json::to_string(&report)
        .unwrap_or_else(|_| "{\"error\": \"JSON serialization failed\"}".to_string())
}

/// Generate 3D coordinates from SMILES (raw distance geometry, no minimization).
/// Returns PDB format string with atoms positioned in 3D space.
///
/// # Example (JS)
/// ```javascript
/// const pdbStr = module.generate_3d_from_smiles("c1ccccc1");
/// console.log(pdbStr);  // PDB file content
/// ```
#[wasm_bindgen]
pub fn generate_3d_from_smiles(smiles: &str) -> Result<String, JsValue> {
    enforce_input_len("smiles", smiles)?;
    let mol = chematic_smiles::parse(smiles)
        .map_err(|e| JsValue::from_str(&format!("SMILES parse error: {}", e)))?;
    let coords = chematic_3d::generate_coords(&mol);
    let pdb_str = chematic_3d::write_pdb(&mol, &coords);
    Ok(pdb_str)
}

/// Generate 3D coordinates and minimize from SMILES string.
/// Pipeline: distance geometry → DREIDING minimization.
/// Better geometry quality than raw DG; suitable for graphics.
///
/// # Example (JS)
/// ```javascript
/// const pdbStr = module.generate_3d_optimized_pdb("c1ccccc1");
/// console.log(pdbStr);  // PDB file with optimized geometry
/// ```
#[wasm_bindgen]
pub fn generate_3d_optimized_pdb(smiles: &str) -> Result<String, JsValue> {
    enforce_input_len("smiles", smiles)?;
    let mol = chematic_smiles::parse(smiles)
        .map_err(|e| JsValue::from_str(&format!("SMILES parse error: {}", e)))?;
    let coords = chematic_3d::generate_and_minimize_dreiding(&mol);
    let pdb_str = chematic_3d::write_pdb(&mol, &coords);
    Ok(pdb_str)
}

#[cfg(test)]
mod tests {
    #[test]
    fn workflow_json_serialization_aspirin() {
        let report = chematic_chem::molecule_report("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let json = serde_json::to_string(&report).unwrap();
        assert!(
            json.contains("molecular_weight"),
            "JSON should contain molecular_weight"
        );
        assert!(json.contains("tpsa"), "JSON should contain tpsa");
    }

    #[test]
    fn workflow_json_serialization_compare() {
        let comparison = chematic_chem::compare_molecules(&["c1ccccc1", "Cc1ccccc1"]).unwrap();
        let json = serde_json::to_string(&comparison).unwrap();
        assert!(json.contains("pairwise"), "JSON should contain pairwise");
        assert!(
            json.contains("ecfp4_tanimoto"),
            "JSON should have similarity metrics"
        );
    }

    #[test]
    fn workflow_json_serialization_screen() {
        let report = chematic_chem::screen_smiles(&["c1ccccc1", "CC", "CCC"]);
        let json = serde_json::to_string(&report).unwrap();
        assert!(json.contains("records"), "JSON should contain records");
        assert!(
            json.contains("maxmin_picks"),
            "JSON should contain diversity picks"
        );
    }

    #[test]
    fn distance_geometry_benzene_has_coords() {
        let mol = chematic_smiles::parse("c1ccccc1").unwrap();
        let coords = chematic_3d::generate_coords(&mol);
        assert_eq!(coords.atom_count(), 6, "benzene has 6 carbons");
    }

    #[test]
    fn distance_geometry_ethane_reasonable_distance() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let coords = chematic_3d::generate_coords(&mol);
        let p0 = coords.get(chematic_core::AtomIdx(0));
        let p1 = coords.get(chematic_core::AtomIdx(1));
        let dist = p0.distance(&p1);
        assert!((dist - 1.54).abs() < 0.15, "C-C distance should be ~1.54 Å");
    }
}
