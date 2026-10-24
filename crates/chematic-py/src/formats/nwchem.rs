use super::*;

/// Read the single Cartesian NWChem GEOMETRY block, preserving source units and center attributes.
#[pyfunction]
pub(super) fn parse_nwchem_geometry<'py>(
    py: Python<'py>,
    #[pyo3(from_py_with=serialized_text)] text: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let record = chematic_mol::parse_nwchem_geometry(text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let json = serde_json::to_string(&record).map_err(|e| PyValueError::new_err(e.to_string()))?;
    json_loads(py, &json)
}

/// Write an edited NWChem geometry from its current structured fields.
#[pyfunction]
pub(super) fn write_nwchem_geometry(record: &Bound<'_, PyDict>) -> PyResult<String> {
    let kwargs = PyDict::new(record.py());
    kwargs.set_item("allow_nan", false)?;
    let json: String = record
        .py()
        .import("json")?
        .call_method("dumps", (record,), Some(&kwargs))?
        .extract()?;
    let record = chematic_mol::nwchem_geometry_from_json(&json)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    chematic_mol::write_nwchem_geometry(&record).map_err(|e| PyValueError::new_err(e.to_string()))
}
