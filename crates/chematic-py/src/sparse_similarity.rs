use chematic_fp::sparse_matrix::{
    MAX_SPARSE_BITS_PER_FINGERPRINT, MAX_SPARSE_FINGERPRINTS, MAX_SPARSE_INPUT_IDS,
    SparseMatrixError as Error, SparseTanimotoMetric as Metric,
};
use pyo3::exceptions::{PyAttributeError, PyMemoryError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyByteArray, PyBytes, PyString};

fn error(value: Error) -> PyErr {
    if value == Error::Allocation {
        PyMemoryError::new_err(value.to_string())
    } else {
        PyValueError::new_err(value.to_string())
    }
}
fn positions(value: &Bound<'_, PyAny>) -> PyResult<()> {
    if value.is_instance_of::<PyString>()
        || value.is_instance_of::<PyBytes>()
        || value.is_instance_of::<PyByteArray>()
    {
        return Err(PyTypeError::new_err(
            "sparse fingerprints must be iterables of integer on-bit positions, not text or bytes",
        ));
    }
    Ok(())
}

/// Checked binary sparse-set Tanimoto matrix, with no dense bit expansion.
/// Unsigned 32-bit on-bit positions or GetOnBits() sources; duplicates/order are ignored. Returns
/// strict lower triangle (1,0), (2,0), (2,1), ... as a list. Empty/empty
/// similarity is 1; distance is 1-similarity. Not count-weighted similarity.
/// Limits: 1024 fingerprints, 65536 raw IDs per row, one million raw IDs
/// total, and 50 million pair-bit work units. Source errors propagate;
/// at most one excess item is consumed to detect limits. Invalid limits or
/// metric raise ValueError; types TypeError; out-of-u32 indices OverflowError;
/// allocation MemoryError. Python bool indices are refused.
#[pyfunction]
#[pyo3(signature = (fingerprints, *, metric="distance"))]
fn sparse_tanimoto_matrix_checked(
    fingerprints: &Bound<'_, PyAny>,
    metric: &str,
) -> PyResult<Vec<f64>> {
    let metric = match metric {
        "distance" => Metric::Distance,
        "similarity" => Metric::Similarity,
        _ => {
            return Err(PyValueError::new_err(
                "metric must be 'distance' or 'similarity'",
            ));
        }
    };
    positions(fingerprints)?;
    let mut rows = Vec::new();
    let mut total = 0;
    for row in fingerprints.try_iter()? {
        let row = row?;
        if rows.len() == MAX_SPARSE_FINGERPRINTS {
            return Err(error(Error::RowLimit));
        }
        // Keep the source row alive while consuming any method-produced view.
        let bits = match row.getattr("GetOnBits") {
            Ok(method) => method.call0()?,
            Err(e) if e.is_instance_of::<PyAttributeError>(row.py()) => row.clone(),
            Err(e) => return Err(e),
        };
        positions(&bits)?;
        let mut ids = Vec::new();
        for value in bits.try_iter()? {
            let value = value?;
            if ids.len() == MAX_SPARSE_BITS_PER_FINGERPRINT {
                return Err(error(Error::BitsPerFingerprint { row: rows.len() }));
            }
            if total == MAX_SPARSE_INPUT_IDS {
                return Err(error(Error::TotalIds));
            }
            if value.is_instance_of::<PyBool>() {
                return Err(PyTypeError::new_err(
                    "bit positions must be integers, not booleans",
                ));
            }
            ids.try_reserve(1).map_err(|_| error(Error::Allocation))?;
            ids.push(value.extract::<u32>()?);
            total += 1;
        }
        rows.try_reserve(1).map_err(|_| error(Error::Allocation))?;
        rows.push(ids);
    }
    chematic_fp::sparse_tanimoto_matrix_checked(&rows, metric).map_err(error)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(sparse_tanimoto_matrix_checked, module)?)
}
