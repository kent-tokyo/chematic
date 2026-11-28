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

#[derive(Clone, Copy)]
struct RepairOptions {
    clash_cutoff: Option<f64>,
    max_alternative_attempts: usize,
    require_regular_geometry: bool,
}

#[derive(Default)]
struct RetryCounts {
    attempts: usize,
    selected_variant: Option<usize>,
    stereo_failures: usize,
    geometry_failures: usize,
    clash_rejections: usize,
    tetrahedral_geometry_rejections: usize,
}

struct RepairRun {
    outcome: Result<
        chematic_3d::stereo_constraints::RepairOutcome,
        chematic_3d::stereo_constraints::RepairFailureReport,
    >,
    initial: chematic_3d::stereo_constraints::StereoVerification,
    comparison: Option<chematic_3d::geometry_check::NonbondedClashComparison>,
    distorted: Option<Vec<chematic_core::AtomIdx>>,
    retry: RetryCounts,
}

fn distorted_after_screen(
    mol: &chematic_core::Molecule,
    coords: &chematic_core::Coords3D,
) -> Vec<chematic_core::AtomIdx> {
    chematic_3d::stereo_constraints::verify_stereo_with_geometry_screen(mol, coords)
        .distorted_tetrahedral
}

/// Try deterministic alternative substituent choices after the default repair
/// was rejected, updating `run` in place when one passes every enabled check.
fn try_alternative_repairs(
    mol: &chematic_core::Molecule,
    geometry: &chematic_core::Coords3D,
    options: RepairOptions,
    cutoff: f64,
    run: &mut RepairRun,
) {
    use chematic_3d::geometry_check::compare_nonbonded_clashes;
    use chematic_3d::stereo_constraints::{
        repair_stereo_with_variant, satisfied_tetrahedral_geometry_is_regular,
    };
    for variant in 1..=options.max_alternative_attempts {
        run.retry.attempts += 1;
        let Ok(alternative) = repair_stereo_with_variant(mol, geometry, variant) else {
            run.retry.stereo_failures += 1;
            continue;
        };
        if !satisfied_tetrahedral_geometry_is_regular(mol, &alternative.coords) {
            run.retry.tetrahedral_geometry_rejections += 1;
            continue;
        }
        let alternative_distorted = options
            .require_regular_geometry
            .then(|| distorted_after_screen(mol, &alternative.coords));
        if alternative_distorted
            .as_ref()
            .is_some_and(|atoms| !atoms.is_empty())
        {
            run.retry.tetrahedral_geometry_rejections += 1;
            continue;
        }
        let Ok(checked) = compare_nonbonded_clashes(mol, geometry, &alternative.coords, cutoff)
        else {
            run.retry.geometry_failures += 1;
            continue;
        };
        if !checked.new_clashes.is_empty() {
            run.retry.clash_rejections += 1;
            continue;
        }
        run.outcome = Ok(alternative);
        run.comparison = Some(checked);
        run.distorted = alternative_distorted;
        run.retry.selected_variant = Some(variant);
        break;
    }
}

fn run_stereo_repair(
    mol: &chematic_core::Molecule,
    geometry: &chematic_core::Coords3D,
    options: RepairOptions,
) -> Result<RepairRun, chematic_3d::geometry_check::GeometryCheckError> {
    use chematic_3d::geometry_check::{
        GeometryCheckConfig, check_geometry, compare_nonbonded_clashes,
    };
    use chematic_3d::stereo_constraints::{repair_stereo, verify_stereo};
    if let Some(cutoff) = options.clash_cutoff {
        let config = GeometryCheckConfig {
            nonbonded_min_distance: cutoff,
            ..GeometryCheckConfig::default()
        };
        check_geometry(mol, geometry, &config)?;
    }
    let initial = verify_stereo(mol, geometry);
    let outcome = repair_stereo(mol, geometry);
    let comparison = match (options.clash_cutoff, &outcome) {
        (Some(cutoff), Ok(repaired)) => Some(compare_nonbonded_clashes(
            mol,
            geometry,
            &repaired.coords,
            cutoff,
        )?),
        _ => None,
    };
    let distorted = match (&outcome, options.require_regular_geometry) {
        (Ok(repaired), true) => Some(distorted_after_screen(mol, &repaired.coords)),
        _ => None,
    };
    let mut retry = RetryCounts {
        selected_variant: outcome.is_ok().then_some(0),
        ..RetryCounts::default()
    };
    let has_distortion = distorted.as_ref().is_some_and(|atoms| !atoms.is_empty());
    if has_distortion {
        retry.selected_variant = None;
    }
    let has_new_clashes = comparison
        .as_ref()
        .is_some_and(|c| !c.new_clashes.is_empty());
    let mut run = RepairRun {
        outcome,
        initial,
        comparison,
        distorted,
        retry,
    };
    if let Some(cutoff) = options.clash_cutoff
        && (has_new_clashes || has_distortion)
    {
        run.retry.selected_variant = None;
        try_alternative_repairs(mol, geometry, options, cutoff, &mut run);
    }
    Ok(run)
}

fn regular_geometry_report<'py>(
    py: Python<'py>,
    distorted: Option<&Vec<chematic_core::AtomIdx>>,
    rejected: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let report = PyDict::new(py);
    if let Some(atoms) = distorted {
        report.set_item("status", if rejected { "rejected" } else { "passed" })?;
        report.set_item(
            "distorted_atoms",
            atoms.iter().map(|a| a.0).collect::<Vec<_>>(),
        )?;
    } else {
        report.set_item("status", "not_evaluated_repair_failed")?;
    }
    Ok(report)
}

fn clash_report<'py>(
    py: Python<'py>,
    cutoff: f64,
    comparison: Option<&chematic_3d::geometry_check::NonbondedClashComparison>,
    rejected: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let report = PyDict::new(py);
    report.set_item("cutoff", cutoff)?;
    let Some(comparison) = comparison else {
        report.set_item("status", "not_evaluated_repair_failed")?;
        return Ok(report);
    };
    report.set_item("status", if rejected { "rejected" } else { "passed" })?;
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
    Ok(report)
}

fn append_unevaluable_from_initial<'py>(
    py: Python<'py>,
    initial: &chematic_3d::stereo_constraints::StereoVerification,
    unevaluable: &Bound<'py, pyo3::types::PyList>,
) -> PyResult<()> {
    use chematic_3d::stereo_constraints::{StereoElement, StereoStatus};
    for atom in &initial.tetrahedral {
        if matches!(atom.status, StereoStatus::Unevaluable(_)) {
            unevaluable.append(stereo_element_dict(
                py,
                StereoElement::Tetrahedral(atom.atom),
            )?)?;
        }
    }
    for bond in &initial.double_bond {
        if matches!(bond.status, StereoStatus::Unevaluable(_)) {
            unevaluable.append(stereo_element_dict(
                py,
                StereoElement::DoubleBond(bond.bond),
            )?)?;
        }
    }
    Ok(())
}

struct OutcomeRows<'a, 'py> {
    repaired: &'a Bound<'py, pyo3::types::PyList>,
    failures: &'a Bound<'py, pyo3::types::PyList>,
    unevaluable: &'a Bound<'py, pyo3::types::PyList>,
}

struct RejectionInfo<'a> {
    clashes: bool,
    distorted: &'a [chematic_core::AtomIdx],
}

type RepairResult = Result<
    chematic_3d::stereo_constraints::RepairOutcome,
    chematic_3d::stereo_constraints::RepairFailureReport,
>;

fn fill_rejected<'py>(
    py: Python<'py>,
    result: &Bound<'py, PyDict>,
    outcome: chematic_3d::stereo_constraints::RepairOutcome,
    rejection: RejectionInfo<'_>,
    rows: &OutcomeRows<'_, 'py>,
) -> PyResult<()> {
    use chematic_3d::stereo_constraints::StereoElement;
    result.set_item("ok", false)?;
    result.set_item("fully_verified", false)?;
    result.set_item("coordinates", py.None())?;
    if rejection.clashes {
        let failure = PyDict::new(py);
        failure.set_item("kind", "geometry")?;
        failure.set_item("reason", "NewNonbondedClashes")?;
        rows.failures.append(failure)?;
    }
    for &atom in rejection.distorted {
        let row = stereo_element_dict(py, StereoElement::Tetrahedral(atom))?;
        row.set_item("reason", "DistortedTetrahedralGeometry")?;
        rows.failures.append(row)?;
    }
    for element in outcome.unevaluable {
        rows.unevaluable.append(stereo_element_dict(py, element)?)?;
    }
    Ok(())
}

fn fill_accepted<'py>(
    py: Python<'py>,
    result: &Bound<'py, PyDict>,
    outcome: chematic_3d::stereo_constraints::RepairOutcome,
    rows: &OutcomeRows<'_, 'py>,
) -> PyResult<()> {
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
        rows.repaired.append(row)?;
    }
    for element in outcome.unevaluable {
        rows.unevaluable.append(stereo_element_dict(py, element)?)?;
    }
    Ok(())
}

fn fill_failed<'py>(
    py: Python<'py>,
    result: &Bound<'py, PyDict>,
    report: chematic_3d::stereo_constraints::RepairFailureReport,
    initial: &chematic_3d::stereo_constraints::StereoVerification,
    rows: &OutcomeRows<'_, 'py>,
) -> PyResult<()> {
    result.set_item("ok", false)?;
    result.set_item("fully_verified", false)?;
    result.set_item("coordinates", py.None())?;
    for (element, reason) in report.failures {
        let row = stereo_element_dict(py, element)?;
        row.set_item("reason", format!("{reason:?}"))?;
        rows.failures.append(row)?;
    }
    append_unevaluable_from_initial(py, initial, rows.unevaluable)
}

fn fill_outcome<'py>(
    py: Python<'py>,
    result: &Bound<'py, PyDict>,
    outcome: RepairResult,
    initial: &chematic_3d::stereo_constraints::StereoVerification,
    rejection: RejectionInfo<'_>,
    rows: &OutcomeRows<'_, 'py>,
) -> PyResult<()> {
    let rejected = rejection.clashes || !rejection.distorted.is_empty();
    match outcome {
        Ok(outcome) if rejected => fill_rejected(py, result, outcome, rejection, rows),
        Ok(outcome) => fill_accepted(py, result, outcome, rows),
        Err(report) => fill_failed(py, result, report, initial, rows),
    }
}

fn retry_report_dict<'py>(py: Python<'py>, retry: &RetryCounts) -> PyResult<Bound<'py, PyDict>> {
    let report = PyDict::new(py);
    report.set_item("attempted", retry.attempts)?;
    report.set_item("selected_variant", retry.selected_variant)?;
    report.set_item("stereo_failures", retry.stereo_failures)?;
    report.set_item("geometry_failures", retry.geometry_failures)?;
    report.set_item("clash_rejections", retry.clash_rejections)?;
    report.set_item(
        "tetrahedral_geometry_rejections",
        retry.tetrahedral_geometry_rejections,
    )?;
    Ok(report)
}

fn build_repair_report<'py>(
    py: Python<'py>,
    run: RepairRun,
    options: RepairOptions,
) -> PyResult<Bound<'py, PyDict>> {
    let RepairRun {
        outcome,
        initial,
        comparison,
        distorted,
        retry,
    } = run;
    let result = PyDict::new(py);
    let repaired_rows = pyo3::types::PyList::empty(py);
    let failures = pyo3::types::PyList::empty(py);
    let unevaluable = pyo3::types::PyList::empty(py);
    result.set_item("retry", retry_report_dict(py, &retry)?)?;
    let rejected_clashes = comparison
        .as_ref()
        .is_some_and(|c| !c.new_clashes.is_empty());
    let rejected_distortion = distorted.as_ref().is_some_and(|atoms| !atoms.is_empty());
    if options.require_regular_geometry {
        let report = regular_geometry_report(py, distorted.as_ref(), rejected_distortion)?;
        result.set_item("regular_geometry_check", report)?;
    } else {
        result.set_item("regular_geometry_check", py.None())?;
    }
    if let Some(cutoff) = options.clash_cutoff {
        let report = clash_report(py, cutoff, comparison.as_ref(), rejected_clashes)?;
        result.set_item("clash_check", report)?;
    } else {
        result.set_item("clash_check", py.None())?;
    }
    let rows = OutcomeRows {
        repaired: &repaired_rows,
        failures: &failures,
        unevaluable: &unevaluable,
    };
    let rejection = RejectionInfo {
        clashes: rejected_clashes,
        distorted: distorted.as_deref().unwrap_or(&[]),
    };
    fill_outcome(py, &result, outcome, &initial, rejection, &rows)?;
    result.set_item("repaired", repaired_rows)?;
    result.set_item("unevaluable", unevaluable)?;
    result.set_item("failures", failures)?;
    Ok(result)
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
    /// With ``max_alternative_attempts`` (at most 16), a rejected repair can
    /// try other deterministic substituent choices; a cutoff is required.
    /// With ``require_regular_geometry=True``, reject output with tetrahedral
    /// centers that fail the optional conservative geometry screen.
    #[pyo3(signature = (coords, *, clash_cutoff=None, max_alternative_attempts=0, require_regular_geometry=false))]
    fn repair_stereo_3d<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
        clash_cutoff: Option<f64>,
        max_alternative_attempts: usize,
        require_regular_geometry: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        if max_alternative_attempts > 16 {
            return Err(PyValueError::new_err(
                "max_alternative_attempts must be at most 16",
            ));
        }
        if max_alternative_attempts > 0 && clash_cutoff.is_none() {
            return Err(PyValueError::new_err(
                "clash_cutoff is required for alternative attempts",
            ));
        }
        if coords.len() != self.inner.atom_count()
            || coords.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(PyValueError::new_err(
                "coordinates must be finite and match the atom count",
            ));
        }
        let geometry = flat_to_coords3d(&coords);
        let mol = self.inner.clone();
        let options = RepairOptions {
            clash_cutoff,
            max_alternative_attempts,
            require_regular_geometry,
        };
        let run = py
            .detach(move || run_stereo_repair(&mol, &geometry, options))
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        build_repair_report(py, run, options)
    }
    /// Report declared atom and double-bond stereo checks at supplied coordinates.
    /// Entries the verifier reports as unsupported are unevaluable. The
    /// verifier does not cover every possible stereochemical specification.
    /// With ``require_regular_geometry=True``, tetrahedral results with extreme
    /// ligand geometry are reported as unevaluable instead of satisfied/violated.
    #[pyo3(signature = (coords, *, require_regular_geometry=false))]
    fn stereo_3d_diagnostics<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
        require_regular_geometry: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        use chematic_3d::stereo_constraints::{verify_stereo, verify_stereo_with_geometry_screen};
        if coords.len() != self.inner.atom_count()
            || coords.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(PyValueError::new_err(
                "coordinates must be finite and match the atom count",
            ));
        }
        let geometry = flat_to_coords3d(&coords);
        let mol = self.inner.clone();
        let (report, distorted) = py.detach(move || {
            if require_regular_geometry {
                let screened = verify_stereo_with_geometry_screen(&mol, &geometry);
                (screened.report, screened.distorted_tetrahedral)
            } else {
                (verify_stereo(&mol, &geometry), Vec::new())
            }
        });
        let result = PyDict::new(py);
        result.set_item("reported", report.n_declared())?;
        result.set_item("satisfied", report.n_satisfied())?;
        result.set_item("violated", report.n_violations())?;
        result.set_item("unevaluable", report.n_unevaluable())?;
        let atoms = pyo3::types::PyList::empty(py);
        for item in report.tetrahedral {
            let row = PyDict::new(py);
            let (status, mut reason) = status_parts(item.status);
            if distorted.contains(&item.atom) {
                reason = Some("DistortedTetrahedralGeometry".to_string());
            }
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
