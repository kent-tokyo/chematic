use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Validate and normalize a bounded ``chematic.nucleic-acid.v1`` document.
///
/// The returned JSON is a stable tagged envelope. ``ok`` is false for typed
/// validation failures; invalid documents are never flattened into molecules.
#[pyfunction]
fn nucleic_acid_validate_json(document_json: &str) -> PyResult<String> {
    serde_json::to_string(&chematic_mol::validate_nucleic_acid_json(
        document_json,
        &chematic_mol::NucleicAcidLimits::default(),
    ))
    .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Apply one bounded metadata edit to a nucleic-acid document.
///
/// This returns the same tagged validation envelope as
/// :func:`nucleic_acid_validate_json` and does not permit topology edits.
#[pyfunction]
fn nucleic_acid_apply_json_command(document_json: &str, command_json: &str) -> PyResult<String> {
    serde_json::to_string(&chematic_mol::apply_nucleic_acid_json_command(
        document_json,
        command_json,
        &chematic_mol::NucleicAcidLimits::default(),
    ))
    .map_err(|error| PyValueError::new_err(error.to_string()))
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(nucleic_acid_validate_json, module)?)?;
    module.add_function(wrap_pyfunction!(nucleic_acid_apply_json_command, module)?)?;
    Ok(())
}
