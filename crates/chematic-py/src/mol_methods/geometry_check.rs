use super::*;

#[pymethods]
impl Mol {
    /// Read-only geometric screen on caller coordinates in atom order (Angstrom).
    /// Returns bond-length outliers, nonbonded overlaps, skipped bond indices and
    /// declared stereo outcomes. No universal validity boolean or minimization.
    #[pyo3(signature = (coordinates, *, bond_relative_tolerance=0.5, nonbonded_min_distance=1.2))]
    fn check_geometry<'py>(
        &self,
        py: Python<'py>,
        coordinates: &Bound<'py, PyAny>,
        bond_relative_tolerance: f64,
        nonbonded_min_distance: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let n = self.inner.atom_count();
        if n == 0 || n > chematic_3d::MAX_GEOMETRY_CHECK_ATOMS {
            return Err(PyValueError::new_err(
                "geometry check requires 1..=2048 atoms",
            ));
        }
        let mut coords = Vec::with_capacity(n);
        for row in coordinates.try_iter()? {
            let row = row?;
            if coords.len() == n {
                return Err(PyValueError::new_err(
                    "coordinates must have one entry per atom",
                ));
            }
            coords.push(row.extract::<[f64; 3]>()?);
        }
        let config = chematic_3d::GeometryCheckConfig {
            bond_relative_tolerance,
            nonbonded_min_distance,
        };
        let report = chematic_3d::check_geometry(&self.inner, &flat_to_coords3d(&coords), &config)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        geometry_report_dict(py, &report)
    }
}

fn geometry_report_dict<'py>(
    py: Python<'py>,
    report: &chematic_3d::GeometryCheckReport,
) -> PyResult<Bound<'py, PyDict>> {
    let out = PyDict::new(py);
    out.set_item("atom_count", report.atom_count)?;
    out.set_item("bonds_checked", report.bonds_checked)?;
    out.set_item(
        "stereo_groups_not_evaluated",
        report.stereo_groups_not_evaluated,
    )?;
    out.set_item(
        "bond_relative_tolerance",
        report.config.bond_relative_tolerance,
    )?;
    out.set_item(
        "nonbonded_min_distance",
        report.config.nonbonded_min_distance,
    )?;
    out.set_item(
        "skipped_bonds",
        report
            .skipped_bonds
            .iter()
            .map(|idx| idx.0)
            .collect::<Vec<_>>(),
    )?;
    out.set_item(
        "max_relative_bond_deviation",
        report.max_relative_bond_deviation,
    )?;
    let bonds = PyList::empty(py);
    for v in &report.bond_violations {
        let item = PyDict::new(py);
        item.set_item("bond", v.bond.0)?;
        item.set_item("atoms", v.atoms.map(|idx| idx.0))?;
        item.set_item("distance", v.distance)?;
        item.set_item("reference_length", v.reference_length)?;
        item.set_item("relative_deviation", v.relative_deviation)?;
        bonds.append(item)?;
    }
    out.set_item("bond_violations", bonds)?;
    let clashes = PyList::empty(py);
    for v in &report.nonbonded_clashes {
        let item = PyDict::new(py);
        item.set_item("atoms", v.atoms.map(|idx| idx.0))?;
        item.set_item("distance", v.distance)?;
        clashes.append(item)?;
    }
    out.set_item("nonbonded_clashes", clashes)?;
    let stereo = crate::pipeline_v2::stereo_verification_dict(py, &report.stereo)?;
    // The pipeline's legacy boolean means only "zero violations" and can be
    // true for unevaluable entries. This API exposes explicit counts instead.
    stereo.del_item("is_fully_satisfied")?;
    out.set_item("stereo", stereo)?;
    Ok(out)
}
