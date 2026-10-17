use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

use crate::Mol;

/// Select a copy of a supplied Mol using the existing RDKit-port default score.
/// Equal scores use the smallest canonical SMILES; identical ties keep the
/// first candidate. Lists, tuples and generators are supported. No enumeration
/// or stereo cleanup is performed; the selected source's metadata is retained.
///
/// Empty input, invalid limits and unsupported candidates raise ValueError.
/// Non-Mol items raise TypeError; iterator exceptions propagate. At most
/// max_candidates + 1 items are consumed, with no partial winner on overflow.
#[pyfunction]
#[pyo3(signature = (candidates, *, max_candidates=1000))]
fn rdkit_pick_tautomer(candidates: &Bound<'_, PyAny>, max_candidates: i64) -> PyResult<Mol> {
    if !(1..=chematic_smiles::rdkit::RDKIT_TAUTOMER_PICK_MAX_CANDIDATES as i64)
        .contains(&max_candidates)
    {
        return Err(PyValueError::new_err("max_candidates must be in 1..=4096"));
    }
    let max_candidates = max_candidates as usize;
    let mut owned = Vec::new();
    for item in candidates.try_iter()? {
        let item = item?;
        if owned.len() == max_candidates {
            return Err(PyValueError::new_err("candidates exceed max_candidates"));
        }
        let mol = item.extract::<PyRef<'_, Mol>>().map_err(|_| {
            PyTypeError::new_err(format!("candidate {} must be a Mol", owned.len()))
        })?;
        // Own both the Arc and the property store before advancing the iterator.
        // A generator may discard its last molecule on the next call.
        owned.push(mol.clone());
    }
    let pick = chematic_smiles::rdkit::rdkit_pick_tautomer(
        owned.iter().map(|m| m.inner.as_ref()),
        max_candidates,
    )
    .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(owned.swap_remove(pick.index))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rdkit_pick_tautomer, m)?)
}
