//! Serialized format text has a separate boundary from identifiers and paths.
use pyo3::exceptions::{PyTypeError, PyUnicodeDecodeError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyBytesMethods, PyString, PyStringMethods};

/// Borrow valid UTF-8 directly from str/immutable bytes for the native call.
/// Never invoke __str__, __bytes__ or __fspath__ on an arbitrary object.
pub(super) fn serialized_text<'a>(value: &'a Bound<'_, PyAny>) -> PyResult<&'a str> {
    if let Ok(text) = value.cast::<PyString>() {
        return text.to_str();
    }
    if let Ok(bytes) = value.cast::<PyBytes>() {
        let data = bytes.as_bytes();
        return core::str::from_utf8(data)
            .map_err(|error| PyUnicodeDecodeError::new_err_from_utf8(value.py(), data, error));
    }
    Err(PyTypeError::new_err(
        "serialized input must be str or bytes containing UTF-8 text",
    ))
}
