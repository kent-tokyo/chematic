use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Stable, typed reason string for a [`chematic_perception::StereoDiagnostic`]
/// -- never a free-form message (see `docs/rfcs/stereo2d_reader_integration_rfc.md`).
pub(crate) fn stereo_reason_str(
    reason: chematic_perception::StereoRejectionReason,
) -> &'static str {
    use chematic_perception::StereoRejectionReason::*;
    match reason {
        ContradictoryWedges => "contradictory_wedges",
        MissingCoordinate => "missing_coordinate",
        DegenerateGeometry => "degenerate_geometry",
        UnsupportedCoordination => "unsupported_coordination",
    }
}

pub(crate) fn stereo_diagnostics_to_py<'py>(
    py: Python<'py>,
    diagnostics: &[chematic_perception::StereoDiagnostic],
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    diagnostics
        .iter()
        .map(|d| {
            let dict = PyDict::new(py);
            dict.set_item("atom_idx", d.atom.0)?;
            dict.set_item("reason", stereo_reason_str(d.reason))?;
            Ok(dict)
        })
        .collect()
}

pub(crate) fn bitvec2048_to_bytes(fp: &chematic_fp::bitvec::BitVec2048) -> Vec<u8> {
    fp.to_le_bytes().to_vec()
}

pub(crate) fn flat_to_coords3d(coords: &[[f64; 3]]) -> chematic_3d::Coords3D {
    chematic_3d::Coords3D {
        points: coords
            .iter()
            .map(|c| chematic_3d::Point3::new(c[0], c[1], c[2]))
            .collect(),
    }
}

/// Convert caller-provided coordinates after validating the public Python
/// boundary. Most 3D kernels index coordinates by molecule atom index, so a
/// short vector would otherwise surface as a Rust panic instead of a Python
/// input error.
pub(crate) fn checked_coords3d(
    mol: &chematic_core::Molecule,
    coords: &[[f64; 3]],
) -> PyResult<chematic_3d::Coords3D> {
    if coords.len() != mol.atom_count() {
        return Err(PyValueError::new_err(format!(
            "expected {} coordinates, got {}",
            mol.atom_count(),
            coords.len()
        )));
    }
    if coords.iter().flatten().any(|value| !value.is_finite()) {
        return Err(PyValueError::new_err("coordinates must all be finite"));
    }
    Ok(flat_to_coords3d(coords))
}
