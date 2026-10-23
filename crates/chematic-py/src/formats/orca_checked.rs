use super::*;
use pyo3::types::{PyBool, PyInt};

fn output_limit(value: &Bound<'_, PyAny>) -> PyResult<usize> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyValueError::new_err(
            "max_output_bytes must be an integer excluding bool",
        ));
    }
    let size: usize = value.extract()?;
    if size == 0 || size > 64 * 1024 * 1024 {
        return Err(PyValueError::new_err("max_output_bytes must be 1..=64 MiB"));
    }
    Ok(size)
}

/// Write a structured ORCA dictionary without Cartesian numeric rounding.
/// Requires comments/keywords/blocks/coords and explicit xyz atom/state fields;
/// refuses unknown fields, bool scalars, conflicting units and lossy metadata.
/// Opaque job settings are retained as represented, not chemically validated.
#[pyfunction]
#[pyo3(signature = (input, *, max_output_bytes=67_108_864))]
pub(super) fn write_orca_input_checked(
    py: Python<'_>,
    input: &Bound<'_, PyDict>,
    #[pyo3(from_py_with=output_limit)] max_output_bytes: usize,
) -> PyResult<String> {
    let options = PyDict::new(py);
    options.set_item("allow_nan", false)?;
    let json: String = py
        .import("json")?
        .call_method("dumps", (input,), Some(&options))?
        .extract()?;
    chematic_mol::write_orca_input_json_checked_with_limits(
        &json,
        &chematic_mol::OrcaWriteLimits {
            max_input_bytes: 16 * 1024 * 1024,
            max_output_bytes,
            max_atoms: 100_000,
        },
    )
    .map_err(|error| PyValueError::new_err(error.to_string()))
}
