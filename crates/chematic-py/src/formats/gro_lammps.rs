use super::gro::parse_gro;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_gro, m)?)?;
    m.add_function(wrap_pyfunction!(gro_to_lammps_atomic, m)?)?;
    Ok(())
}

/// GRO to LAMMPS atomic data with explicit real/metal units, types and masses.
/// Cell and coordinates convert nm -> Angstroms together. Velocity conversion
/// uses the selected time unit. Original source/identity are retained separately.
#[pyfunction]
#[pyo3(signature = (text, atom_types, masses, *, units))]
pub(super) fn gro_to_lammps_atomic<'py>(
    py: Python<'py>,
    text: &str,
    atom_types: Vec<i64>,
    masses: Vec<f64>,
    units: &str,
) -> PyResult<Bound<'py, PyDict>> {
    use chematic_mol::gro_lammps::GroLammpsUnits;
    let units = match units {
        "real" => GroLammpsUnits::Real,
        "metal" => GroLammpsUnits::Metal,
        _ => return Err(PyValueError::new_err("units must be real or metal")),
    };
    let report = chematic_mol::gro_lammps::gro_to_lammps_atomic(text, &atom_types, &masses, units)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let identity = report
        .source_frame
        .atoms
        .iter()
        .enumerate()
        .map(|(i, atom)| {
            let d = PyDict::new(py);
            d.set_item("lammps_id", i + 1)?;
            d.set_item("source_row", i)?;
            d.set_item("atom_number", atom.atom_number)?;
            d.set_item("atom_name", &atom.atom_name)?;
            d.set_item("residue_number", atom.residue_number)?;
            d.set_item("residue_name", &atom.residue_name)?;
            Ok(d)
        })
        .collect::<PyResult<Vec<_>>>()?;
    let result = PyDict::new(py);
    result.set_item("data_text", chematic_mol::write_lammps_data(&report.data))?;
    result.set_item("units", units.as_str())?;
    result.set_item("atom_style", "atomic")?;
    result.set_item("length_scale", 10.0)?;
    result.set_item("velocity_scale", units.velocity_scale())?;
    result.set_item("atom_identity", identity)?;
    result.set_item("source_text", &report.source_frame.source_text)?;
    Ok(result)
}
