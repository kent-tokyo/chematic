//! High-level workflow APIs for WASM.
//!
//! Provides single-molecule reports, multi-molecule comparisons, and batch screening
//! as JSON strings suitable for JavaScript/TypeScript integration.

use wasm_bindgen::prelude::*;

const WORKFLOW_MAX_INPUT_BYTES: usize = 1_000_000;
const WORKFLOW_MAX_BATCH_ITEMS: usize = 1_024;
const WORKFLOW_MAX_ATOMS: usize = 10_000;

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

fn enforce_batch_molecule_sizes(smiles: &[&str]) -> Result<(), String> {
    for (idx, value) in smiles.iter().enumerate() {
        if let Ok(mol) = chematic_smiles::parse(value.trim())
            && mol.atom_count() > WORKFLOW_MAX_ATOMS
        {
            return Err(format!(
                "smiles_batch[{idx}] exceeds maximum atom count ({} > {})",
                mol.atom_count(),
                WORKFLOW_MAX_ATOMS
            ));
        }
    }
    Ok(())
}

fn enforce_molecule_size(smiles: &str, label: &str) -> Result<(), JsValue> {
    let mol = chematic_smiles::parse(smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
    if mol.atom_count() > WORKFLOW_MAX_ATOMS {
        return Err(JsValue::from_str(&format!(
            "{label} exceeds maximum atom count ({} > {})",
            mol.atom_count(),
            WORKFLOW_MAX_ATOMS
        )));
    }
    Ok(())
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
    enforce_molecule_size(smiles, "smiles")?;
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
    enforce_molecule_size(smiles1, "smiles1")?;
    enforce_molecule_size(smiles2, "smiles2")?;
    let smiles = [smiles1, smiles2];
    let comparison =
        chematic_chem::compare_molecules(&smiles).map_err(|e| JsValue::from_str(&e.to_string()))?;
    serde_json::to_string(&comparison)
        .map_err(|e| JsValue::from_str(&format!("JSON serialization failed: {}", e)))
}

/// Compare multiple SMILES strings (up to 256 by default).
/// Accepts a delimiter-separated list (e.g., newline or comma).
/// The input is limited to 1 MiB, 1,024 records, and 10,000 atoms per parsed
/// molecule.
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
    enforce_batch_molecule_sizes(&smiles_vec).map_err(|error| JsValue::from_str(&error))?;
    let comparison = chematic_chem::compare_molecules(&smiles_vec)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    serde_json::to_string(&comparison)
        .map_err(|e| JsValue::from_str(&format!("JSON serialization failed: {}", e)))
}

/// Screen a batch of SMILES strings (JSON string output).
/// Returns per-record results including pass/fail with error details.
/// Includes MaxMin diversity picking and Butina clustering by default.
/// The input is limited to 1 MiB, 1,024 records, and 10,000 atoms per parsed
/// molecule.
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
    if let Err(error) = enforce_batch_molecule_sizes(&smiles_vec) {
        return format!("{{\"error\":\"{}\"}}", error.replace('"', "\\\""));
    }
    let report = chematic_chem::screen_smiles(&smiles_vec);
    serde_json::to_string(&report)
        .unwrap_or_else(|_| "{\"error\": \"JSON serialization failed\"}".to_string())
}

/// Canonicalize a bounded delimiter-separated SMILES batch.
///
/// Each result retains its input index and original text. Invalid records are
/// returned inline with `status: "rejected"`; later records are still
/// processed in deterministic input order. The manifest includes explicit
/// outcome accounting so callers do not infer full success from a zero error
/// list while inputs were skipped elsewhere.
#[wasm_bindgen]
pub fn canonicalize_smiles_batch_json(
    smiles_batch: &str,
    delimiter: &str,
) -> Result<String, JsValue> {
    let smiles_vec = split_bounded_batch(smiles_batch, delimiter, "smiles_batch")?;
    let canonicalizer = chematic_smiles::SmilesBatchCanonicalizer::default();
    let result = canonicalizer.canonicalize_with_result(smiles_vec);
    let all_succeeded = result.all_succeeded();
    let output: Vec<serde_json::Value> = result
        .records
        .into_iter()
        .map(|record| match record.result {
            chematic_smiles::BatchCanonicalization::Accepted { canonical_smiles } => {
                serde_json::json!({
                    "input_index": record.input_index,
                    "input": record.input,
                    "status": "accepted",
                    "canonical_smiles": canonical_smiles,
                    "error": null,
                    "error_stage": null,
                })
            }
            chematic_smiles::BatchCanonicalization::Rejected { error } => serde_json::json!({
                "input_index": record.input_index,
                "input": record.input,
                "status": "rejected",
                "error": error,
                "error_stage": "parse",
            }),
        })
        .collect();
    let manifest = serde_json::json!({
        "schema_version": 1,
        "operation": "canonicalize_smiles",
        "status": "complete",
        "record_count": output.len(),
        "input_count": result.input_count,
        "accepted_count": result.accepted_count,
        "rejected_count": result.rejected_count,
        "refused_count": result.refused_count,
        "skipped_count": result.skipped_count,
        "all_succeeded": all_succeeded,
        "records": output,
    });
    let json = serde_json::to_string(&manifest)
        .map_err(|error| JsValue::from_str(&format!("JSON serialization failed: {error}")))?;
    if json.len() > WORKFLOW_MAX_INPUT_BYTES {
        return Err(JsValue::from_str(
            "batch result exceeds maximum output size",
        ));
    }
    Ok(json)
}

/// Stateful, bounded adapter for an unknown-length SMILES source.
///
/// Callers distinguish observation from processing: this lets a Worker report
/// rows already read but not yet processed when it is cancelled.  `finish_json`
/// is the only normal EOF path; `stop_json` reports an unknown unread suffix
/// and never fabricates skipped or successful rows.
#[wasm_bindgen]
pub struct SmilesBatchStreamHandle {
    stream: Option<chematic_smiles::SmilesBatchStream>,
}

impl Default for SmilesBatchStreamHandle {
    fn default() -> Self {
        Self::new()
    }
}

fn stream_record_json(record: &chematic_smiles::BatchCanonicalRecord) -> serde_json::Value {
    match &record.result {
        chematic_smiles::BatchCanonicalization::Accepted { canonical_smiles } => {
            serde_json::json!({
                "input_index": record.input_index,
                "input": record.input,
                "status": "accepted",
                "canonical_smiles": canonical_smiles,
                "error": null,
                "error_stage": null,
            })
        }
        chematic_smiles::BatchCanonicalization::Rejected { error } => serde_json::json!({
            "input_index": record.input_index,
            "input": record.input,
            "status": "rejected",
            "error": error,
            "error_stage": "parse",
        }),
    }
}

fn stream_manifest_json(
    result: chematic_smiles::BatchCanonicalStreamResult,
) -> Result<String, JsValue> {
    let complete = result.stream_complete;
    let manifest = serde_json::json!({
        "schema_version": 1,
        "operation": "canonicalize_smiles_stream",
        "input_kind": "unknown_length_stream",
        "status": if complete { "complete" } else { "incomplete" },
        "observed_input_count": result.observed_input_count,
        "completed_count": result.records.len(),
        "accepted_count": result.accepted_count,
        "rejected_count": result.rejected_count,
        "refused_count": 0,
        "skipped_count": 0,
        "unprocessed_observed_count": result.unprocessed_observed_count,
        "unread_input": if complete { serde_json::json!(0) } else { serde_json::json!("unknown") },
        "all_succeeded": result.all_succeeded(),
        "terminal_reason": result.terminal_reason.map(|reason| reason.as_str()),
        "records": result.records.iter().map(stream_record_json).collect::<Vec<_>>(),
    });
    let json = serde_json::to_string(&manifest)
        .map_err(|error| JsValue::from_str(&format!("JSON serialization failed: {error}")))?;
    if json.len() > WORKFLOW_MAX_INPUT_BYTES {
        return Err(JsValue::from_str(
            "stream batch result exceeds maximum output size",
        ));
    }
    Ok(json)
}

#[wasm_bindgen]
impl SmilesBatchStreamHandle {
    /// Start an unknown-length batch.  A handle is terminal after `finish_json`
    /// or `stop_json` and cannot be reused for a second source.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            stream: Some(chematic_smiles::SmilesBatchCanonicalizer::default().stream()),
        }
    }

    /// Register one row observed from the source without processing it yet.
    pub fn observe(&mut self, smiles: &str) -> Result<usize, JsValue> {
        enforce_input_len("stream smiles", smiles)?;
        let stream = self
            .stream
            .as_mut()
            .ok_or_else(|| JsValue::from_str("stream is already terminal"))?;
        Ok(stream.observe(smiles))
    }

    /// Process one previously observed row, returning `null` when none remain.
    pub fn process_next_json(&mut self) -> Result<String, JsValue> {
        let stream = self
            .stream
            .as_mut()
            .ok_or_else(|| JsValue::from_str("stream is already terminal"))?;
        serde_json::to_string(&stream.process_next().map(stream_record_json))
            .map_err(|error| JsValue::from_str(&format!("JSON serialization failed: {error}")))
    }

    /// Mark EOF and process every observed pending row into a complete manifest.
    pub fn finish_json(&mut self) -> Result<String, JsValue> {
        let stream = self
            .stream
            .take()
            .ok_or_else(|| JsValue::from_str("stream is already terminal"))?;
        stream_manifest_json(stream.finish())
    }

    /// Stop before EOF with a typed reason and an explicitly unknown unread suffix.
    pub fn stop_json(&mut self, terminal_reason: &str) -> Result<String, JsValue> {
        let terminal_reason = terminal_reason
            .parse::<chematic_smiles::StreamTerminalReason>()
            .map_err(JsValue::from_str)?;
        let stream = self
            .stream
            .take()
            .ok_or_else(|| JsValue::from_str("stream is already terminal"))?;
        stream_manifest_json(stream.stop(terminal_reason))
    }
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
    if mol.atom_count() > WORKFLOW_MAX_ATOMS {
        return Err(JsValue::from_str(
            "smiles exceeds maximum atom count (10000)",
        ));
    }
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
    if mol.atom_count() > WORKFLOW_MAX_ATOMS {
        return Err(JsValue::from_str(
            "smiles exceeds maximum atom count (10000)",
        ));
    }
    let coords = chematic_3d::generate_and_minimize_dreiding(&mol);
    let pdb_str = chematic_3d::write_pdb(&mol, &coords);
    Ok(pdb_str)
}

#[cfg(test)]
mod tests {
    use super::{SmilesBatchStreamHandle, WORKFLOW_MAX_ATOMS, enforce_batch_molecule_sizes};

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
    fn workflow_batch_atom_limit_is_explicit() {
        let huge = "C".repeat(WORKFLOW_MAX_ATOMS + 1);
        let err = enforce_batch_molecule_sizes(&[&huge]).unwrap_err();
        assert!(err.contains("maximum atom count"));
    }

    #[test]
    fn stream_handle_preserves_observed_prefix_at_cancellation() {
        let mut stream = SmilesBatchStreamHandle::new();
        assert_eq!(stream.observe("CCO").unwrap(), 0);
        assert_eq!(stream.observe("C1CC").unwrap(), 1);
        assert_eq!(stream.observe("CCN").unwrap(), 2);
        let processed: serde_json::Value =
            serde_json::from_str(&stream.process_next_json().unwrap()).unwrap();
        assert_eq!(processed["input_index"], 0);
        assert_eq!(processed["status"], "accepted");

        let manifest: serde_json::Value =
            serde_json::from_str(&stream.stop_json("cancelled").unwrap()).unwrap();
        assert_eq!(manifest["input_kind"], "unknown_length_stream");
        assert_eq!(manifest["status"], "incomplete");
        assert_eq!(manifest["observed_input_count"], 3);
        assert_eq!(manifest["completed_count"], 1);
        assert_eq!(manifest["unprocessed_observed_count"], 2);
        assert_eq!(manifest["unread_input"], "unknown");
        assert_eq!(manifest["terminal_reason"], "cancelled");
        assert_eq!(manifest["all_succeeded"], false);
    }

    #[test]
    fn stream_handle_finish_processes_pending_rows() {
        let mut stream = SmilesBatchStreamHandle::new();
        stream.observe("CCO").unwrap();
        stream.observe("C1CC").unwrap();
        let manifest: serde_json::Value =
            serde_json::from_str(&stream.finish_json().unwrap()).unwrap();
        assert_eq!(manifest["status"], "complete");
        assert_eq!(manifest["completed_count"], 2);
        assert_eq!(manifest["accepted_count"], 1);
        assert_eq!(manifest["rejected_count"], 1);
        assert_eq!(manifest["unprocessed_observed_count"], 0);
        assert_eq!(manifest["unread_input"], 0);
        assert_eq!(manifest["terminal_reason"], serde_json::Value::Null);
        assert_eq!(manifest["all_succeeded"], false);
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
