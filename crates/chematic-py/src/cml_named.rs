use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

/// Write explicit unique CML atom IDs and optional Cartesian Angstrom triples.
/// IDs follow atom insertion order and are never renamed. Raises ValueError
/// for unsupported graph semantics, invalid/duplicate IDs, counts or numerics.
/// This exports a graph; source-only MOL2 metadata is outside this contract.
#[pyfunction]
#[pyo3(signature = (mol, atom_ids, coordinates=None))]
fn write_cml_named_checked(
    mol: &crate::Mol,
    atom_ids: &Bound<'_, PyAny>,
    coordinates: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let n = mol.inner.atom_count();
    if n == 0 || n > 4096 || mol.inner.bond_count() > 8192 {
        return Err(PyValueError::new_err(
            "named CML requires 1..=4096 atoms and at most 8192 bonds",
        ));
    }
    if atom_ids.is_instance_of::<PyString>() || atom_ids.is_instance_of::<PyBytes>() {
        return Err(PyTypeError::new_err(
            "atom_ids must be an iterable of str IDs",
        ));
    }
    let mut ids = Vec::new();
    for item in atom_ids.try_iter()? {
        let item = item?;
        if ids.len() == n {
            return Err(PyValueError::new_err("atom_ids must have one ID per atom"));
        }
        ids.push(item.extract::<String>()?);
    }
    let coords = coordinates
        .map(|rows| {
            let mut output = Vec::new();
            for row in rows.try_iter()? {
                let row = row?;
                if output.len() == n {
                    return Err(PyValueError::new_err(
                        "coordinates must have one 3D row per atom",
                    ));
                }
                output.push(row.extract::<[f64; 3]>()?);
            }
            Ok::<_, PyErr>(output)
        })
        .transpose()?;
    chematic_mol::write_cml_named_checked(&mol.inner, &ids, coords.as_deref())
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(write_cml_named_checked, m)?)
}
