use super::*;
use pyo3::types::PyAnyMethods;
fn limits(
    bytes: Option<usize>,
    atoms: Option<usize>,
    line: Option<usize>,
) -> chematic_mol::ViewmolParseLimits {
    let d = chematic_mol::ViewmolParseLimits::default();
    chematic_mol::ViewmolParseLimits {
        max_input_bytes: bytes.unwrap_or(d.max_input_bytes),
        max_atoms: atoms.unwrap_or(d.max_atoms),
        max_line_len: line.unwrap_or(d.max_line_len),
    }
}
fn failure(e: chematic_mol::ViewmolError) -> PyErr {
    PyValueError::new_err(e.to_string())
}
fn record_from_py(
    py: Python<'_>,
    record: &Bound<'_, PyDict>,
    limits: &chematic_mol::ViewmolParseLimits,
) -> PyResult<chematic_mol::ViewmolRecord> {
    let kwargs = PyDict::new(py);
    kwargs.set_item("allow_nan", false)?;
    let json: String = py
        .import("json")?
        .call_method("dumps", (record,), Some(&kwargs))?
        .extract()?;
    chematic_mol::parse_viewmol_json(&json, limits).map_err(failure)
}
/// Read bounded ViewMol geometry, retaining raw scale, coordinates and fixed axes.
#[pyfunction]
#[pyo3(signature=(text, *, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_viewmol<'py>(
    py: Python<'py>,
    #[pyo3(from_py_with = serialized_text)] text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let value = chematic_mol::parse_viewmol_with_limits(
        text,
        &limits(max_input_bytes, max_atoms, max_line_len),
    )
    .map_err(failure)?;
    let json = serde_json::to_string(&value).map_err(|e| PyValueError::new_err(e.to_string()))?;
    materials::json_loads(py, &json)
}
/// Write current validated fields, without coordinate conversion or field loss.
#[pyfunction]
#[pyo3(signature=(record, *, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn write_viewmol(
    py: Python<'_>,
    record: &Bound<'_, PyDict>,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<String> {
    let l = limits(max_input_bytes, max_atoms, max_line_len);
    let r = record_from_py(py, record, &l)?;
    chematic_mol::write_viewmol_with_limits(&r, &l).map_err(failure)
}
/// Explicit binary64 coordinate projection; reject overflow/underflow to zero.
#[pyfunction]
#[pyo3(signature=(record, *, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn viewmol_positions_angstrom(
    py: Python<'_>,
    record: &Bound<'_, PyDict>,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Vec<[f64; 3]>> {
    let l = limits(max_input_bytes, max_atoms, max_line_len);
    record_from_py(py, record, &l)?
        .positions_angstrom_with_limits(&l)
        .map_err(failure)
}
