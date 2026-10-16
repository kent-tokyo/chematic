use super::*;

fn checked_integer(value: &Bound<'_, PyAny>) -> PyResult<u32> {
    if value.is_instance_of::<pyo3::types::PyBool>() {
        return Err(PyValueError::new_err(
            "cleanup indices/budgets must be integers, not bool",
        ));
    }
    value.extract::<u32>().map_err(|_| {
        PyValueError::new_err("cleanup indices/budgets must be nonnegative 32-bit integers")
    })
}

#[pymethods]
impl Mol {
    /// Remove one tetrahedral declaration only with an exact, target-fixed
    /// graph witness preserving all other stereo and source labels. Returns
    /// molecule/removed/atom_index/search_visits/witness. No witness means
    /// conservative retention. Unsupported profiles or exhausted budgets raise
    /// ValueError; the source is unchanged and no partial cleanup is returned.
    #[pyo3(signature=(atom_idx, *, max_visits=1000000))]
    fn remove_redundant_tetrahedral_stereo<'py>(
        &self,
        py: Python<'py>,
        #[pyo3(from_py_with=checked_integer)] atom_idx: u32,
        #[pyo3(from_py_with=checked_integer)] max_visits: u32,
    ) -> PyResult<Bound<'py, PyDict>> {
        let r = chematic_chem::remove_redundant_tetrahedral_stereo(
            &self.inner,
            chematic_core::AtomIdx(atom_idx),
            u64::from(max_visits),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let out = PyDict::new(py);
        out.set_item(
            "molecule",
            Mol {
                inner: Arc::new(r.molecule),
                props: self.props.clone(),
            },
        )?;
        out.set_item("removed", r.removed)?;
        out.set_item("atom_index", r.atom.0)?;
        out.set_item("search_visits", r.search_visits)?;
        out.set_item(
            "witness",
            r.witness
                .map(|v| v.into_iter().map(|a| a.0).collect::<Vec<_>>()),
        )?;
        Ok(out)
    }
}
