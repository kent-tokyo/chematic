use super::*;

/// Convert PQR into rigid PDBQT with explicit source-order AutoDock types.
/// Radii remain only in the returned sidecar; consent is required. This does
/// not perceive bonds, recalculate charges, remove H or build ligand torsions.
#[pyfunction]
#[pyo3(signature = (text, atom_types, *, allow_radius_loss=false))]
pub(super) fn pqr_to_pdbqt<'py>(
    py: Python<'py>,
    text: &str,
    atom_types: Vec<String>,
    allow_radius_loss: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let source = chematic_mol::parse_pqr(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let types: Vec<&str> = atom_types.iter().map(String::as_str).collect();
    let report = chematic_mol::pqr_to_pdbqt(&source, &types, allow_radius_loss)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let result = PyDict::new(py);
    result.set_item("pdbqt", report.pdbqt)?;
    result.set_item("source_charges", report.source_charges)?;
    result.set_item("source_radii", report.source_radii)?;
    result.set_item("source_coords", report.source_coords)?;
    result.set_item("radii_omitted", report.radii_omitted)?;
    result.set_item("charge_provenance", "pqr_atom_records")?;
    result.set_item("atom_type_provenance", "caller_supplied")?;
    result.set_item(
        "max_charge_rounding_error",
        report.max_charge_rounding_error,
    )?;
    result.set_item(
        "max_coordinate_rounding_error",
        report.max_coordinate_rounding_error,
    )?;
    Ok(result)
}
