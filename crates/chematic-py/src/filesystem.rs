//! Python filesystem paths, independent of serialized chemistry text.
use pyo3::prelude::*;
use std::path::PathBuf;

/// Use Python's filesystem decoding so bytes and surrogate-escaped POSIX names
/// stay lossless. `fsdecode` calls `fspath`, rejecting arbitrary objects rather
/// than accepting their `__str__`; PathBuf then preserves native OS encoding.
pub(crate) fn path(value: &Bound<'_, PyAny>) -> PyResult<PathBuf> {
    value
        .py()
        .import("os")?
        .call_method1("fsdecode", (value,))?
        .extract()
}

/// Decode an optional output path with the same contract as required paths.
pub(crate) fn optional_path(value: &Bound<'_, PyAny>) -> PyResult<Option<PathBuf>> {
    if value.is_none() {
        Ok(None)
    } else {
        path(value).map(Some)
    }
}
