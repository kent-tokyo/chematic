//! Packed fingerprint folding with explicit width and termination information.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyBytesMethods, PyDict};

/// OR-fold a packed fingerprint until the target density or minimum length.
/// Returns bytes, width, density, fold count and whether the target was reached.
#[pyfunction]
#[pyo3(signature = (fp, density=0.3, min_length=64, *, bit_width=None))]
fn fold_to_target_density<'py>(
    py: Python<'py>,
    fp: &Bound<'py, PyBytes>,
    density: f64,
    min_length: usize,
    bit_width: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let packed = fp.as_bytes();
    let width = match bit_width {
        Some(width) => width,
        None => packed
            .len()
            .checked_mul(8)
            .ok_or_else(|| PyValueError::new_err("fold density bit width overflow"))?,
    };
    let result = chematic_fp::fold_packed_to_target_density(packed, width, density, min_length)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    let out = PyDict::new(py);
    out.set_item(
        "fingerprint",
        PyBytes::new(py, &result.fingerprint.to_le_bytes()),
    )?;
    out.set_item("bit_width", result.fingerprint.bit_width())?;
    out.set_item("original_bits", result.original_bits)?;
    out.set_item("on_bit_count", result.fingerprint.popcount())?;
    out.set_item("density", result.density())?;
    out.set_item("folds", result.folds)?;
    out.set_item("target_reached", result.target_reached)?;
    out.set_item(
        "stop_reason",
        if result.target_reached {
            "target_reached"
        } else {
            "minimum_length"
        },
    )?;
    Ok(out)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_pyfunctions!(m, fold_to_target_density)
}
