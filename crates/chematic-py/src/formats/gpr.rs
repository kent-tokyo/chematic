use super::*;
use pyo3::types::PyAnyMethods;
fn limits(
    bytes: Option<usize>,
    atoms: Option<usize>,
    bonds: Option<usize>,
    line: Option<usize>,
) -> chematic_mol::GprParseLimits {
    let defaults = chematic_mol::GprParseLimits::default();
    chematic_mol::GprParseLimits {
        max_input_bytes: bytes.unwrap_or(defaults.max_input_bytes),
        max_atoms: atoms.unwrap_or(defaults.max_atoms),
        max_bonds: bonds.unwrap_or(defaults.max_bonds),
        max_line_len: line.unwrap_or(defaults.max_line_len),
    }
}
fn failure(e: chematic_mol::GprError) -> PyErr {
    PyValueError::new_err(e.to_string())
}
/// Read gpr 100 / one coordinate set. Coordinates and partial charges retain decimal strings.
#[pyfunction]
#[pyo3(signature=(text, *, max_input_bytes=None, max_atoms=None, max_bonds=None, max_line_len=None))]
pub(super) fn parse_gpr<'py>(
    py: Python<'py>,
    #[pyo3(from_py_with = serialized_text)] text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_bonds: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let record = chematic_mol::parse_gpr_with_limits(
        text,
        &limits(max_input_bytes, max_atoms, max_bonds, max_line_len),
    )
    .map_err(failure)?;
    let json = serde_json::to_string(&record).map_err(|e| PyValueError::new_err(e.to_string()))?;
    materials::json_loads(py, &json)
}
/// Write current checked fields; absent partial charges remain absent.
#[pyfunction]
#[pyo3(signature=(record, *, max_input_bytes=None, max_atoms=None, max_bonds=None, max_line_len=None))]
pub(super) fn write_gpr(
    py: Python<'_>,
    record: &Bound<'_, PyDict>,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_bonds: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<String> {
    let kwargs = PyDict::new(py);
    kwargs.set_item("allow_nan", false)?;
    let json: String = py
        .import("json")?
        .call_method("dumps", (record,), Some(&kwargs))?
        .extract()?;
    let limits = limits(max_input_bytes, max_atoms, max_bonds, max_line_len);
    let record = chematic_mol::parse_gpr_json(&json, &limits).map_err(failure)?;
    chematic_mol::write_gpr_with_limits(&record, &limits).map_err(failure)
}
