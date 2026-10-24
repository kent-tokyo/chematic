use super::*;

/// Read the formatted Classic DL-POLY HISTORY trajectory, preserving all frames and raw units.
#[pyfunction]
pub(super) fn parse_dlpoly_history<'py>(
    py: Python<'py>,
    #[pyo3(from_py_with=serialized_text)] text: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let record = chematic_mol::parse_dlpoly_history(text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let json = serde_json::to_string(&record).map_err(|e| PyValueError::new_err(e.to_string()))?;
    json_loads(py, &json)
}

/// Write an edited DL-POLY HISTORY trajectory from its current structured fields.
#[pyfunction]
pub(super) fn write_dlpoly_history(record: &Bound<'_, PyDict>) -> PyResult<String> {
    let kwargs = PyDict::new(record.py());
    kwargs.set_item("allow_nan", false)?;
    let json: String = record
        .py()
        .import("json")?
        .call_method("dumps", (record,), Some(&kwargs))?
        .extract()?;
    let record = chematic_mol::dlpoly_history_from_json(&json)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    chematic_mol::write_dlpoly_history(&record).map_err(|e| PyValueError::new_err(e.to_string()))
}
