use chematic_chem::descriptor_distances::{
    DescriptorDistanceError as Error, MAX_DESCRIPTOR_FEATURES, MAX_DESCRIPTOR_ROWS,
};
use pyo3::exceptions::{PyMemoryError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyBytes, PyString};

fn error(value: Error) -> PyErr {
    if value == Error::Allocation {
        PyMemoryError::new_err(value.to_string())
    } else {
        PyValueError::new_err(value.to_string())
    }
}

fn numeric_iterable(value: &Bound<'_, PyAny>) -> PyResult<()> {
    if value.is_instance_of::<PyString>()
        || value.is_instance_of::<PyBytes>()
        || value.is_instance_of::<PyByteArray>()
    {
        return Err(PyTypeError::new_err(
            "descriptor matrix and rows must be numeric iterables, not text or bytes",
        ));
    }
    Ok(())
}

/// Checked condensed Euclidean distances from rectangular finite vectors.
/// Lower triangle: (1,0), (2,0), (2,1), ...; no normalization or unit changes.
/// Accepts nested numeric iterables/NumPy rows. Maximum 1024 rows, 256 features,
/// 50 million pair-feature operations. Iterators consume at most one excess
/// entry per dimension; source exceptions propagate. Invalid/range/limit
/// inputs raise ValueError; wrong types TypeError; allocation MemoryError.
/// Empty/singleton input returns []; zero features yield zero distances.
#[pyfunction]
fn euclidean_distance_matrix_checked(descriptor_matrix: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    numeric_iterable(descriptor_matrix)?;
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for row in descriptor_matrix.try_iter()? {
        let row = row?;
        if rows.len() == MAX_DESCRIPTOR_ROWS {
            return Err(error(Error::RowLimit));
        }
        numeric_iterable(&row)?;
        let mut values = Vec::new();
        for value in row.try_iter()? {
            let value = value?;
            if values.len() == MAX_DESCRIPTOR_FEATURES {
                return Err(error(Error::FeatureLimit));
            }
            values
                .try_reserve(1)
                .map_err(|_| error(Error::Allocation))?;
            values.push(value.extract::<f64>()?);
        }
        if let Some(first) = rows.first()
            && values.len() != first.len()
        {
            return Err(error(Error::Ragged {
                row: rows.len(),
                expected: first.len(),
                actual: values.len(),
            }));
        }
        rows.try_reserve(1).map_err(|_| error(Error::Allocation))?;
        rows.push(values);
    }
    chematic_chem::descriptor_distances::euclidean_distance_matrix_checked(&rows).map_err(error)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(euclidean_distance_matrix_checked, module)?)
}
