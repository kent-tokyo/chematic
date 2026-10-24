use super::*;
use pyo3::types::PyAnyMethods;
fn limits(
    bytes: Option<usize>,
    atoms: Option<usize>,
    line: Option<usize>,
    frames: Option<usize>,
) -> chematic_mol::Gromos96ParseLimits {
    let d = chematic_mol::Gromos96ParseLimits::default();
    chematic_mol::Gromos96ParseLimits {
        max_input_bytes: bytes.unwrap_or(d.max_input_bytes),
        max_atoms: atoms.unwrap_or(d.max_atoms),
        max_line_len: line.unwrap_or(d.max_line_len),
        max_frames: frames.unwrap_or(d.max_frames),
    }
}
fn failure(e: chematic_mol::Gromos96Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}
fn record_from_py(
    py: Python<'_>,
    record: &Bound<'_, PyDict>,
    limits: &chematic_mol::Gromos96ParseLimits,
) -> PyResult<chematic_mol::Gromos96Document> {
    let kwargs = PyDict::new(py);
    kwargs.set_item("allow_nan", false)?;
    let json: String = py
        .import("json")?
        .call_method("dumps", (record,), Some(&kwargs))?
        .extract()?;
    chematic_mol::parse_gromos96_json(&json, limits).map_err(failure)
}
/// Read bounded GROMOS96 geometry, retaining validated blocks, atom labels, raw fixed-width fields and comments.
#[pyfunction]
#[pyo3(signature=(text, *, length_unit="unspecified", max_input_bytes=None, max_atoms=None, max_line_len=None, max_frames=None))]
pub(super) fn parse_gromos96<'py>(
    py: Python<'py>,
    #[pyo3(from_py_with = serialized_text)] text: &str,
    length_unit: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
    max_frames: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let value = chematic_mol::parse_gromos96_with_limits(
        text,
        length_unit.parse().map_err(failure)?,
        &limits(max_input_bytes, max_atoms, max_line_len, max_frames),
    )
    .map_err(failure)?;
    let json = serde_json::to_string(&value).map_err(|e| PyValueError::new_err(e.to_string()))?;
    materials::json_loads(py, &json)
}
/// Write current validated fields, without coordinate conversion or field loss.
#[pyfunction]
#[pyo3(signature=(record, *, max_input_bytes=None, max_atoms=None, max_line_len=None, max_frames=None))]
pub(super) fn write_gromos96(
    py: Python<'_>,
    record: &Bound<'_, PyDict>,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
    max_frames: Option<usize>,
) -> PyResult<String> {
    let l = limits(max_input_bytes, max_atoms, max_line_len, max_frames);
    let r = record_from_py(py, record, &l)?;
    chematic_mol::write_gromos96_with_limits(&r, &l).map_err(failure)
}
/// Explicit binary64 coordinate projection; reject overflow/underflow to zero.
#[pyfunction]
#[pyo3(signature=(record, *, max_input_bytes=None, max_atoms=None, max_line_len=None, max_frames=None))]
pub(super) fn gromos96_geometry_angstrom<'py>(
    py: Python<'py>,
    record: &Bound<'_, PyDict>,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
    max_frames: Option<usize>,
) -> PyResult<Bound<'py, PyAny>> {
    let l = limits(max_input_bytes, max_atoms, max_line_len, max_frames);
    let geometry = record_from_py(py, record, &l)?
        .geometry_angstrom_with_limits(&l)
        .map_err(failure)?;
    let json =
        serde_json::to_string(&geometry).map_err(|e| PyValueError::new_err(e.to_string()))?;
    py.import("json")?.call_method1("loads", (json,))
}
