use super::*;

/// Read one DL-POLY CONFIG/REVCON with optional velocities and forces.
/// Numbers retain source units; FIELD/CONTROL context is not inferred.
#[pyfunction]
pub(super) fn parse_dlpoly_config<'py>(
    py: Python<'py>,
    #[pyo3(from_py_with = serialized_text)] text: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let record = chematic_mol::parse_dlpoly_config(text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let json = serde_json::to_string(&record).map_err(|e| PyValueError::new_err(e.to_string()))?;
    json_loads(py, &json)
}

/// Serialize edited CONFIG data; refuse unknown fields or inconsistent levels.
#[pyfunction]
pub(super) fn write_dlpoly_config(record: &Bound<'_, PyDict>) -> PyResult<String> {
    let kwargs = PyDict::new(record.py());
    kwargs.set_item("allow_nan", false)?;
    let json: String = record
        .py()
        .import("json")?
        .call_method("dumps", (record,), Some(&kwargs))?
        .extract()?;
    let record = chematic_mol::dlpoly_config_from_json(&json)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    chematic_mol::write_dlpoly_config(&record).map_err(|e| PyValueError::new_err(e.to_string()))
}
