use super::*;

#[pymethods]
impl Mol {
    /// Checked SASA with explicit Angstrom coordinates and one positive radius
    /// per atom, including explicit H. Returns total/per_atom areas in Å².
    /// No coordinate generation or default radii. Empty inputs return zero;
    /// invalid input, ambiguous coincident spheres and work/numeric limits
    /// raise ValueError. This is native Shrake-Rupley, not FreeSASA parity.
    #[pyo3(signature = (coordinates, radii, *, probe_radius=1.4, sphere_points=960))]
    fn sasa_with_radii_checked<'py>(
        &self,
        py: Python<'py>,
        coordinates: &Bound<'py, PyAny>,
        radii: &Bound<'py, PyAny>,
        probe_radius: f64,
        sphere_points: i64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let n = self.inner.atom_count();
        if n > 2048 {
            return Err(PyValueError::new_err(
                "checked SASA supports at most 2048 atoms",
            ));
        }
        if !(32..=65536).contains(&sphere_points) {
            return Err(PyValueError::new_err("sphere_points must be in 32..=65536"));
        }
        let mut coords = Vec::new();
        for row in coordinates.try_iter()? {
            let row = row?;
            if coords.len() == n {
                return Err(PyValueError::new_err(
                    "coordinates must have one entry per atom",
                ));
            }
            coords.push(row.extract::<[f64; 3]>()?);
        }
        let mut values = Vec::new();
        for radius in radii.try_iter()? {
            let radius = radius?;
            if values.len() == n {
                return Err(PyValueError::new_err("radii must have one entry per atom"));
            }
            values.push(radius.extract::<f64>()?);
        }
        let result = chematic_3d::sasa_with_radii_checked(
            &self.inner,
            &coords,
            &values,
            probe_radius,
            sphere_points as usize,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let out = PyDict::new(py);
        out.set_item("total", result.total)?;
        out.set_item("per_atom", result.per_atom)?;
        Ok(out)
    }
}
