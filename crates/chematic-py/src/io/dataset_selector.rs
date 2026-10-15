use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyInt, PyList, PySlice, PyTuple};

pub(super) fn positive_integer(value: &Bound<'_, PyAny>) -> PyResult<usize> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyTypeError::new_err("expected an integer, excluding bool"));
    }
    value
        .extract::<usize>()
        .ok()
        .filter(|v| *v > 0)
        .ok_or_else(|| PyValueError::new_err("expected a positive integer fitting usize"))
}
pub(super) fn scalar_index(value: &Bound<'_, PyAny>, len: usize) -> PyResult<usize> {
    if value.is_instance_of::<PyBool>() || !value.is_instance_of::<PyInt>() {
        return Err(PyTypeError::new_err(
            "record index must be an integer, excluding bool",
        ));
    }
    let index = value
        .extract::<isize>()
        .map_err(|_| PyIndexError::new_err("SDF record index out of range"))?;
    super::dataset_index(index, len)
}
pub(super) fn error_policy(errors: &str) -> PyResult<bool> {
    match errors {
        "keep" => Ok(false),
        "raise" => Ok(true),
        _ => Err(PyValueError::new_err("errors must be keep or raise")),
    }
}
fn limit_error(actual: usize, limit: usize) -> PyErr {
    super::sdf_dataset_error(chematic_mol::MolParseError::ResourceLimit {
        resource: "selected records",
        actual,
        limit,
    })
}
pub(super) fn indices(key: &Bound<'_, PyAny>, len: usize, limit: usize) -> PyResult<Vec<usize>> {
    if let Ok(slice) = key.cast::<PySlice>() {
        let len = isize::try_from(len)
            .map_err(|_| PyIndexError::new_err("dataset too large for Python indexing"))?;
        let s = slice.indices(len)?;
        let count = s.slicelength;
        if count > limit {
            return Err(limit_error(count, limit));
        }
        let mut out = Vec::with_capacity(count);
        let mut current = s.start;
        for k in 0..count {
            out.push(current as usize);
            if k + 1 < count {
                current = current
                    .checked_add(s.step)
                    .ok_or_else(|| PyIndexError::new_err("slice overflow"))?;
            }
        }
        return Ok(out);
    }
    if key.is_instance_of::<PyInt>() && !key.is_instance_of::<PyBool>() {
        return Ok(vec![scalar_index(key, len)?]);
    }
    if !key.is_instance_of::<PyList>() && !key.is_instance_of::<PyTuple>() {
        return Err(PyTypeError::new_err(
            "selector must be an int, slice, list or tuple",
        ));
    }
    let count = key.len()?;
    if count == 0 {
        return Ok(Vec::new());
    }
    if key.get_item(0)?.is_instance_of::<PyBool>() {
        return mask_indices(key, len, limit);
    }
    if count > limit {
        return Err(limit_error(count, limit));
    }
    key.try_iter()?
        .map(|item| scalar_index(&item?, len))
        .collect()
}
fn mask_indices(key: &Bound<'_, PyAny>, len: usize, limit: usize) -> PyResult<Vec<usize>> {
    if key.len()? != len {
        return Err(PyValueError::new_err(
            "boolean mask length must match dataset",
        ));
    }
    let mut out = Vec::new();
    for (i, item) in key.try_iter()?.enumerate() {
        let item = item?;
        if !item.is_instance_of::<PyBool>() {
            return Err(PyTypeError::new_err(
                "boolean masks cannot contain integer indices",
            ));
        }
        if item.extract::<bool>()? {
            if out.len() >= limit {
                return Err(limit_error(out.len() + 1, limit));
            }
            out.push(i);
        }
    }
    Ok(out)
}
