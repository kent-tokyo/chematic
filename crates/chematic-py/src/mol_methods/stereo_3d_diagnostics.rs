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

fn stereo_element_dict<'py>(
    py: Python<'py>,
    element: chematic_3d::stereo_constraints::StereoElement,
) -> PyResult<Bound<'py, PyDict>> {
    use chematic_3d::stereo_constraints::StereoElement;
    let row = PyDict::new(py);
    match element {
        StereoElement::Tetrahedral(atom) => {
            row.set_item("kind", "atom")?;
            row.set_item("index", atom.0)?;
        }
        StereoElement::DoubleBond(bond) => {
            row.set_item("kind", "bond")?;
            row.set_item("index", bond.0)?;
        }
    }
    Ok(row)
}

#[pymethods]
impl Mol {
    /// Attempt to repair declared 3D stereo at supplied coordinates.
    ///
    /// Returns ``ok``, ``fully_verified``, ``coordinates``, ``repaired``,
    /// ``unevaluable``, and ``failures``. Failed repairs return no coordinates.
    /// ``fully_verified`` refers only to stereo elements reported by the verifier;
    /// it is false when any reported element was unevaluable. A molecule with no
    /// reported stereo is vacuously fully verified. This operation does not
    /// minimize the structure or guarantee the absence of steric clashes.
    /// Set ``clash_cutoff`` to a positive distance in Å to reject newly close
    /// nonbonded atom pairs and receive their identities in ``clash_check``.
    #[pyo3(signature = (coords, *, clash_cutoff=None))]
    fn repair_stereo_3d<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
        clash_cutoff: Option<f64>,
    ) -> PyResult<Bound<'py, PyDict>> {
        use chematic_3d::geometry_check::{
            GeometryCheckConfig, check_geometry, compare_nonbonded_clashes,
        };
        use chematic_3d::stereo_constraints::{repair_stereo, verify_stereo};
        if coords.len() != self.inner.atom_count()
            || coords.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(PyValueError::new_err(
                "coordinates must be finite and match the atom count",
            ));
        }
        let geometry = flat_to_coords3d(&coords);
        let mol = self.inner.clone();
        let (outcome, initial, comparison) = py
            .detach(move || {
                if let Some(cutoff) = clash_cutoff {
                    let config = GeometryCheckConfig {
                        nonbonded_min_distance: cutoff,
                        ..GeometryCheckConfig::default()
                    };
                    check_geometry(&mol, &geometry, &config)?;
                }
                let initial = verify_stereo(&mol, &geometry);
                let outcome = repair_stereo(&mol, &geometry);
                let comparison = match (clash_cutoff, &outcome) {
                    (Some(cutoff), Ok(repaired)) => Some(compare_nonbonded_clashes(
                        &mol,
                        &geometry,
                        &repaired.coords,
                        cutoff,
                    )?),
                    (Some(_), Err(_)) => None,
                    (None, _) => None,
                };
                Ok::<_, chematic_3d::geometry_check::GeometryCheckError>((
                    outcome, initial, comparison,
                ))
            })
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let result = PyDict::new(py);
        let repaired_rows = pyo3::types::PyList::empty(py);
        let failures = pyo3::types::PyList::empty(py);
        let unevaluable = pyo3::types::PyList::empty(py);
        let rejected_clashes = comparison
            .as_ref()
            .is_some_and(|c| !c.new_clashes.is_empty());
        if let Some(cutoff) = clash_cutoff {
            let report = PyDict::new(py);
            report.set_item("cutoff", cutoff)?;
            if let Some(comparison) = &comparison {
                report.set_item(
                    "status",
                    if rejected_clashes {
                        "rejected"
                    } else {
                        "passed"
                    },
                )?;
                report.set_item("before_count", comparison.before_count)?;
                report.set_item("after_count", comparison.after_count)?;
                for (key, clashes) in [
                    ("new", &comparison.new_clashes),
                    ("resolved", &comparison.resolved_clashes),
                ] {
                    let rows = pyo3::types::PyList::empty(py);
                    for clash in clashes {
                        let row = PyDict::new(py);
                        row.set_item("atoms", [clash.atoms[0].0, clash.atoms[1].0])?;
                        row.set_item("distance", clash.distance)?;
                        rows.append(row)?;
                    }
                    report.set_item(key, rows)?;
                }
            } else {
                report.set_item("status", "not_evaluated_repair_failed")?;
            }
            result.set_item("clash_check", report)?;
        } else {
            result.set_item("clash_check", py.None())?;
        }
        match outcome {
            Ok(outcome) if rejected_clashes => {
                result.set_item("ok", false)?;
                result.set_item("fully_verified", false)?;
                result.set_item("coordinates", py.None())?;
                let failure = PyDict::new(py);
                failure.set_item("kind", "geometry")?;
                failure.set_item("reason", "NewNonbondedClashes")?;
                failures.append(failure)?;
                for element in outcome.unevaluable {
                    unevaluable.append(stereo_element_dict(py, element)?)?;
                }
            }
            Ok(outcome) => {
                result.set_item("ok", true)?;
                result.set_item("fully_verified", outcome.unevaluable.is_empty())?;
                result.set_item(
                    "coordinates",
                    outcome
                        .coords
                        .points
                        .iter()
                        .map(|p| [p.x, p.y, p.z])
                        .collect::<Vec<_>>(),
                )?;
                for repair in outcome.repaired {
                    let row = stereo_element_dict(py, repair.element)?;
                    row.set_item("atoms_moved", repair.atoms_moved)?;
                    row.set_item("max_displacement", repair.max_displacement)?;
                    repaired_rows.append(row)?;
                }
                for element in outcome.unevaluable {
                    unevaluable.append(stereo_element_dict(py, element)?)?;
                }
            }
            Err(outcome) => {
                result.set_item("ok", false)?;
                result.set_item("fully_verified", false)?;
                result.set_item("coordinates", py.None())?;
                for (element, reason) in outcome.failures {
                    let row = stereo_element_dict(py, element)?;
                    row.set_item("reason", format!("{reason:?}"))?;
                    failures.append(row)?;
                }
                for atom in initial.tetrahedral {
                    if matches!(
                        atom.status,
                        chematic_3d::stereo_constraints::StereoStatus::Unevaluable(_)
                    ) {
                        unevaluable.append(stereo_element_dict(
                            py,
                            chematic_3d::stereo_constraints::StereoElement::Tetrahedral(atom.atom),
                        )?)?;
                    }
                }
                for bond in initial.double_bond {
                    if matches!(
                        bond.status,
                        chematic_3d::stereo_constraints::StereoStatus::Unevaluable(_)
                    ) {
                        unevaluable.append(stereo_element_dict(
                            py,
                            chematic_3d::stereo_constraints::StereoElement::DoubleBond(bond.bond),
                        )?)?;
                    }
                }
            }
        }
        result.set_item("repaired", repaired_rows)?;
        result.set_item("unevaluable", unevaluable)?;
        result.set_item("failures", failures)?;
        Ok(result)
    }
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
