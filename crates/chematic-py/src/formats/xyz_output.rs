use super::*;
use pyo3::types::{PyBool, PyByteArray, PyBytes, PyInt, PyString};

fn integer(value: &Bound<'_, PyAny>) -> PyResult<usize> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyValueError::new_err(
            "expected a nonnegative integer, excluding bool",
        ));
    }
    value
        .extract()
        .map_err(|_| PyValueError::new_err("integer must fit usize"))
}
fn not_text(value: &Bound<'_, PyAny>) -> PyResult<()> {
    if value.is_instance_of::<PyString>()
        || value.is_instance_of::<PyBytes>()
        || value.is_instance_of::<PyByteArray>()
    {
        return Err(PyValueError::new_err(
            "geometry rows must be non-text iterables",
        ));
    }
    Ok(())
}
fn position(row: &Bound<'_, PyAny>) -> PyResult<[f64; 3]> {
    not_text(row)?;
    let mut values = row.try_iter()?;
    let mut out = [0.0; 3];
    for value in &mut out {
        let v = values
            .next()
            .transpose()?
            .ok_or_else(|| PyValueError::new_err("position requires exactly 3 coordinates"))?;
        if v.is_instance_of::<PyBool>() {
            return Err(PyValueError::new_err("coordinates exclude bool"));
        }
        *value = v.extract()?;
    }
    if values.next().transpose()?.is_some() {
        return Err(PyValueError::new_err(
            "position requires exactly 3 coordinates",
        ));
    }
    Ok(out)
}
fn atoms(
    elements: &Bound<'_, PyAny>,
    positions: &Bound<'_, PyAny>,
) -> PyResult<Vec<chematic_mol::XyzAtom>> {
    not_text(elements)?;
    not_text(positions)?;
    let mut positions = positions.try_iter()?;
    let mut atoms = Vec::new();
    for symbol in elements.try_iter()? {
        if atoms.len() >= chematic_mol::xyz_output::MAX_XYZ_OUTPUT_ATOMS {
            return Err(PyValueError::new_err("atom count exceeds 100000"));
        }
        let symbol = symbol?;
        let text: &str = symbol.extract()?;
        let element = chematic_core::Element::from_symbol(text)
            .filter(|e| e.symbol() == text)
            .ok_or_else(|| PyValueError::new_err("expected a canonical element symbol"))?;
        let row = positions
            .next()
            .transpose()?
            .ok_or_else(|| PyValueError::new_err("element and position counts differ"))?;
        let [x, y, z] = position(&row)?;
        atoms
            .try_reserve(1)
            .map_err(|e| PyValueError::new_err(format!("atom allocation failed: {e}")))?;
        atoms.push(chematic_mol::XyzAtom { element, x, y, z });
    }
    if positions.next().transpose()?.is_some() {
        return Err(PyValueError::new_err("element and position counts differ"));
    }
    Ok(atoms)
}

/// Write plain XYZ coordinates with explicit decimal precision and minimum width.
/// Finite input, canonical elements, bounded geometry/output; no bond inference.
#[pyfunction]
#[pyo3(signature=(elements,positions_angstrom,*,comment="",notation="fixed",decimal_places=12,min_width=20,max_output_bytes=67108864))]
#[allow(clippy::too_many_arguments)]
pub(super) fn write_xyz_checked(
    elements: &Bound<'_, PyAny>,
    positions_angstrom: &Bound<'_, PyAny>,
    comment: &str,
    notation: &str,
    #[pyo3(from_py_with=integer)] decimal_places: usize,
    #[pyo3(from_py_with=integer)] min_width: usize,
    #[pyo3(from_py_with=integer)] max_output_bytes: usize,
) -> PyResult<String> {
    let options = chematic_mol::XyzWriteOptions {
        decimal_places,
        min_width,
        max_output_bytes,
        notation: match notation {
            "fixed" => chematic_mol::XyzNotation::Fixed,
            "scientific" => chematic_mol::XyzNotation::Scientific,
            _ => {
                return Err(PyValueError::new_err(
                    "notation must be 'fixed' or 'scientific'",
                ));
            }
        },
    };
    options
        .validate(comment)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let atoms = atoms(elements, positions_angstrom)?;
    chematic_mol::write_xyz_checked(&atoms, comment, &options)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_pyfunctions!(
        m,
        from_xyz,
        from_extxyz,
        from_extxyz_all,
        to_extxyz,
        write_xyz_checked
    )
}
