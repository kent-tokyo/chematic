use super::*;

/// Parse a single Tinker XYZ frame into an editable record, retaining opaque
/// names/types, explicit neighbors and optional cell parameters. No perception.
#[pyfunction]
pub(super) fn parse_tinker_xyz<'py>(
    py: Python<'py>,
    #[pyo3(from_py_with = serialized_text)] text: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let record = chematic_mol::tinker::parse_tinker_xyz(text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let json = serde_json::to_string(&record).map_err(|e| PyValueError::new_err(e.to_string()))?;
    json_loads(py, &json)
}

/// Serialize current record fields, refusing unknown keys or invalid topology.
#[pyfunction]
pub(super) fn write_tinker_xyz(record: &Bound<'_, PyDict>) -> PyResult<String> {
    let kwargs = PyDict::new(record.py());
    kwargs.set_item("allow_nan", false)?;
    let json: String = record
        .py()
        .import("json")?
        .call_method("dumps", (record,), Some(&kwargs))?
        .extract()?;
    let record = chematic_mol::tinker::tinker_xyz_from_json(&json)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    chematic_mol::tinker::write_tinker_xyz(&record)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
