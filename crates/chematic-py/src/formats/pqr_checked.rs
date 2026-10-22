use super::*;
use pyo3::types::{PyBool, PyByteArray, PyBytes, PyInt, PyString};

fn record(row: &Bound<'_, PyAny>) -> PyResult<chematic_mol::PqrAtomRecord> {
    let dict = row
        .cast::<PyDict>()
        .map_err(|_| PyValueError::new_err("PQR records must be dictionaries"))?;
    for (key, _) in dict.iter() {
        let key: &str = key.extract()?;
        if ![
            "group_pdb",
            "serial",
            "atom_name",
            "res_name",
            "chain_id",
            "res_seq",
            "icode",
            "x",
            "y",
            "z",
            "charge",
            "radius",
            "element",
        ]
        .contains(&key)
        {
            return Err(PyValueError::new_err("unknown PQR atom field"));
        }
    }
    for key in [
        "group_pdb",
        "serial",
        "atom_name",
        "res_name",
        "res_seq",
        "x",
        "y",
        "z",
        "charge",
        "radius",
    ] {
        let value = dict
            .get_item(key)?
            .filter(|v| !v.is_none())
            .ok_or_else(|| {
                PyValueError::new_err(format!("missing required PQR atom field '{key}'"))
            })?;
        if value.is_instance_of::<PyBool>() {
            return Err(PyValueError::new_err("PQR scalar fields exclude bool"));
        }
    }
    // Bound text before the legacy extraction helper copies Python strings.
    for (key, limit) in [
        ("group_pdb", 6),
        ("atom_name", 64),
        ("res_name", 64),
        ("chain_id", 64),
        ("icode", 1),
        ("element", 3),
    ] {
        if let Some(value) = dict.get_item(key)?.filter(|v| !v.is_none()) {
            let text = value.cast::<PyString>()?;
            if text.len()? > limit {
                return Err(PyValueError::new_err(format!(
                    "PQR text field '{key}' exceeds limit"
                )));
            }
        }
    }
    super::materials::pydict_to_pqr_atom(dict)
}

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

/// Write explicit PQR atom fields without coordinate/charge/radius rounding.
/// Requires identity and numeric fields; optional chain/icode/element are retained
/// or inferred. Refuses unknown fields, bool scalars and unrepresentable metadata.
/// Bounded iterables, 100000 atoms and at most 64 MiB output; no bond inference.
#[pyfunction]
#[pyo3(signature = (atoms, *, max_output_bytes=67_108_864))]
pub(super) fn write_pqr_checked(
    atoms: &Bound<'_, PyAny>,
    #[pyo3(from_py_with=output_limit)] max_output_bytes: usize,
) -> PyResult<String> {
    if atoms.is_instance_of::<PyString>()
        || atoms.is_instance_of::<PyBytes>()
        || atoms.is_instance_of::<PyByteArray>()
    {
        return Err(PyValueError::new_err(
            "PQR atoms must be a non-text iterable",
        ));
    }
    let mut records = Vec::new();
    for row in atoms.try_iter()? {
        if records.len() >= 100_000 {
            return Err(PyValueError::new_err("PQR output exceeds 100000 atoms"));
        }
        records.push(record(&row?)?);
    }
    chematic_mol::write_pqr_checked_with_limits(
        &records,
        &chematic_mol::PqrWriteLimits {
            max_atoms: 100_000,
            max_output_bytes,
        },
    )
    .map_err(|error| PyValueError::new_err(error.to_string()))
}
