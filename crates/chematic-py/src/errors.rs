//! Typed input errors for the parsers.
//!
//! `ChematicInputError` subclasses `ValueError` (so `except ValueError` keeps
//! working) and carries the abstention `category` (`malformed`,
//! `unsupported`, `ambiguous` or `resource_limit`), a stable `code`, and the
//! `format` that was being read, so callers need not match message text.

use pyo3::prelude::*;

pyo3::create_exception!(
    chematic,
    ChematicInputError,
    pyo3::exceptions::PyValueError,
    "Input chematic could not read. Subclasses ValueError. Attributes: \
     ``category`` (``\"malformed\"``, ``\"unsupported\"``, ``\"ambiguous\"`` or \
     ``\"resource_limit\"``), ``code`` (stable identifier, e.g. \
     ``\"smiles_parse\"``, ``\"input_too_large\"``) and ``format`` (e.g. \
     ``\"smiles\"``, ``\"mmcif\"``)."
);

/// A [`ChematicInputError`] with its `category`, `code` and `format` set.
pub(crate) fn input_error(category: &str, code: &str, format: &str, message: String) -> PyErr {
    let err = ChematicInputError::new_err(message);
    Python::attach(|py| {
        let value = err.value(py);
        // Setting attributes on a fresh exception instance cannot fail.
        let _ = value.setattr("category", category);
        let _ = value.setattr("code", code);
        let _ = value.setattr("format", format);
    });
    err
}

/// A malformed-input error for `format` (code `<format>_parse`).
pub(crate) fn malformed(format: &str, message: String) -> PyErr {
    input_error("malformed", &format!("{format}_parse"), format, message)
}
