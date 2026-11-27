use super::*;

fn status_parts(
    status: chematic_3d::stereo_constraints::StereoStatus,
) -> (&'static str, Option<String>) {
    use chematic_3d::stereo_constraints::StereoStatus;
    match status {
        StereoStatus::Satisfied => ("satisfied", None),
        StereoStatus::Violated => ("violated", None),
        StereoStatus::Unevaluable(reason) => ("unevaluable", Some(format!("{reason:?}"))),
    }
}

#[pymethods]
impl Mol {
    /// Report declared atom and double-bond stereo checks at supplied coordinates.
    /// Entries the verifier reports as unsupported are unevaluable. The
    /// verifier does not cover every possible stereochemical specification.
    fn stereo_3d_diagnostics<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
    ) -> PyResult<Bound<'py, PyDict>> {
        use chematic_3d::stereo_constraints::verify_stereo;
        if coords.len() != self.inner.atom_count()
            || coords.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(PyValueError::new_err(
                "coordinates must be finite and match the atom count",
            ));
        }
        let geometry = flat_to_coords3d(&coords);
        let mol = self.inner.clone();
        let report = py.detach(move || verify_stereo(&mol, &geometry));
        let result = PyDict::new(py);
        result.set_item("reported", report.n_declared())?;
        result.set_item("satisfied", report.n_satisfied())?;
        result.set_item("violated", report.n_violations())?;
        result.set_item("unevaluable", report.n_unevaluable())?;
        let atoms = pyo3::types::PyList::empty(py);
        for item in report.tetrahedral {
            let row = PyDict::new(py);
            let (status, reason) = status_parts(item.status);
            row.set_item("atom", item.atom.0)?;
            row.set_item("status", status)?;
            row.set_item("reason", reason)?;
            atoms.append(row)?;
        }
        result.set_item("atom_stereo", atoms)?;
        let bonds = pyo3::types::PyList::empty(py);
        for item in report.double_bond {
            let row = PyDict::new(py);
            let (status, reason) = status_parts(item.status);
            row.set_item("bond", item.bond.0)?;
            row.set_item("status", status)?;
            row.set_item("reason", reason)?;
            bonds.append(row)?;
        }
        result.set_item("double_bond_stereo", bonds)?;
        Ok(result)
    }
}
