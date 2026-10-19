//! CML text readers exposed to Python.
use super::serialized_input::serialized_text;
use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

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
