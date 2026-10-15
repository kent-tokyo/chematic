use super::*;

/// Checked single-model PDB coordinates. Long names never truncate. The source
/// sidecar keeps original atom/residue aliases and precision; no bonds or TER.
#[pyfunction]
#[pyo3(signature=(text, *, namespace="auth", allow_metadata_loss=false))]
pub(super) fn mmcif_to_pdb_checked<'py>(
    py: Python<'py>,
    text: &str,
    namespace: &str,
    allow_metadata_loss: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let naming = match namespace {
        "auth" => chematic_mol::PdbNamespace::Author,
        "label" => chematic_mol::PdbNamespace::Label,
        _ => return Err(PyValueError::new_err("namespace must be 'auth' or 'label'")),
    };
    let source = chematic_mol::parse_mmcif_residues(text)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let report = chematic_mol::write_mmcif_pdb_checked(&source, naming, allow_metadata_loss)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let result = PyDict::new(py);
    result.set_item("pdb", report.pdb)?;
    result.set_item("namespace", namespace)?;
    result.set_item("source_text", text)?;
    result.set_item("losses", report.losses)?;
    result.set_item(
        "source",
        mmcif_residues::mmcif_residues_to_pydict(py, &source)?,
    )?;
    result.set_item(
        "max_coordinate_rounding_error",
        report.max_coordinate_rounding_error,
    )?;
    result.set_item(
        "max_occupancy_rounding_error",
        report.max_occupancy_rounding_error,
    )?;
    result.set_item("max_b_iso_rounding_error", report.max_b_iso_rounding_error)?;
    Ok(result)
}
