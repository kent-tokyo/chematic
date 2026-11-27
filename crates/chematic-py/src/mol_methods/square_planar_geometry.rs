use super::*;

#[pymethods]
impl Mol {
    /// Measure declared square-planar ligand angles without judging validity.
    fn square_planar_geometry_diagnostics<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        use chematic_3d::square_planar_geometry::square_planar_angles;
        use chematic_core::SquarePlanarPermutation;
        if coords.len() != self.inner.atom_count()
            || coords.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(PyValueError::new_err(
                "coordinates must be finite and match the atom count",
            ));
        }
        let geometry = flat_to_coords3d(&coords);
        let reports = square_planar_angles(&self.inner, &geometry)
            .ok_or_else(|| PyValueError::new_err("invalid coordinates"))?;
        reports
            .into_iter()
            .map(|report| {
                let item = PyDict::new(py);
                item.set_item("center", report.center.0)?;
                let tag = report.permutation.map(|p| match p {
                    SquarePlanarPermutation::SP1 => "SP1",
                    SquarePlanarPermutation::SP2 => "SP2",
                    SquarePlanarPermutation::SP3 => "SP3",
                });
                item.set_item("permutation", tag)?;
                item.set_item("neighbors", report.neighbors)?;
                item.set_item(
                    "trans",
                    report
                        .trans
                        .into_iter()
                        .map(|a| (a.atoms[0], a.atoms[1], a.degrees))
                        .collect::<Vec<_>>(),
                )?;
                item.set_item(
                    "cis",
                    report
                        .cis
                        .into_iter()
                        .map(|a| (a.atoms[0], a.atoms[1], a.degrees))
                        .collect::<Vec<_>>(),
                )?;
                item.set_item("unavailable", report.unavailable)?;
                Ok(item)
            })
            .collect()
    }
}
