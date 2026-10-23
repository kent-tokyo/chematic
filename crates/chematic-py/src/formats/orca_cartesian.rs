use super::*;

/// Read one Cartesian ORCA geometry with explicit coordinate units, without rescaling.
/// Angs/Bohrs are supported; missing units default to Angstrom. Conflicting or
/// unsupported units, %coords overrides and non-Cartesian/external geometries
/// are refused. Atom extras are opaque; no bonds or atom charges are inferred.
#[pyfunction]
#[pyo3(signature=(text,*,max_input_bytes=16_777_216,max_atoms=100_000))]
pub(super) fn parse_orca_cartesian<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: usize,
    max_atoms: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let geometry = chematic_mol::parse_orca_cartesian_with_limits(
        text,
        &chematic_mol::OrcaInputParseLimits {
            max_input_bytes,
            max_atoms,
            ..Default::default()
        },
    )
    .map_err(|error| PyValueError::new_err(error.to_string()))?;
    let json = geometry.to_json().to_string();
    if json.len() > 16_777_216 {
        return Err(PyValueError::new_err("ORCA Cartesian JSON exceeds 16 MiB"));
    }
    py.import("json")?
        .call_method1("loads", (json,))?
        .cast_into::<PyDict>()
        .map_err(Into::into)
}
