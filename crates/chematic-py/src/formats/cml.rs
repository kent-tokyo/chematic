//! CML text readers exposed to Python.
use super::serialized_input::serialized_text;
use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

type CmlWith3d = (Mol, Option<Vec<Vec<f64>>>);

/// Parse one CML molecule and return ``(Mol, coords_or_none)`` in source order.
/// Cartesian Angstrom coordinates must be complete and finite on every atom,
/// or absent on all atoms. No geometry is generated. Requires unique atom IDs.
/// The lightweight XML profile does not interpret metadata or coordinate stereo.
#[pyfunction]
pub(super) fn from_cml_with_3d(
    #[pyo3(from_py_with = serialized_text)] cml_str: &str,
) -> PyResult<CmlWith3d> {
    let limits = chematic_mol::CmlParseLimits {
        max_input_bytes: 16 * 1024 * 1024,
        max_atoms: 100_000,
        ..Default::default()
    };
    chematic_mol::parse_cml_with_3d_with_limits(cml_str, &limits)
        .map(|(mol, coords)| {
            (
                Mol {
                    num_rings_cache: Default::default(),
                    inner: Arc::new(mol),
                    props: Default::default(),
                },
                coords.map(|rows| rows.into_iter().map(|row| row.to_vec()).collect()),
            )
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Chemical Markup Language (CML) string into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("molecule.cml") as f:
///         mol = chematic.from_cml(f.read())
///
/// Accepts ``str`` or UTF-8 ``bytes``; malformed UTF-8 raises ``UnicodeDecodeError``.
#[pyfunction]
pub(super) fn from_cml(#[pyo3(from_py_with = serialized_text)] cml_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_cml(cml_str)
        .map(|(mol, _coords)| Mol {
            num_rings_cache: Default::default(),
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a structurally valid, non-empty CML string into a ``Mol`` object.
///
/// Unlike :func:`from_cml`, this opt-in entry point rejects missing molecules,
/// empty molecules, and mismatched or unclosed XML elements.
///
/// Accepts ``str`` or UTF-8 ``bytes``; malformed UTF-8 raises ``UnicodeDecodeError``.
#[pyfunction]
pub(super) fn from_cml_strict(
    #[pyo3(from_py_with = serialized_text)] cml_str: &str,
) -> PyResult<Mol> {
    chematic_mol::parse_cml_strict(cml_str)
        .map(|(mol, _coords)| Mol {
            num_rings_cache: Default::default(),
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
