use super::*;

/// Parse bounded RDFile v1 records, retaining raw structures and ordered data.
#[pyfunction]
#[pyo3(signature=(text, *, max_input_bytes=None, max_records=None, max_fields_per_record=None, max_line_len=None, max_reaction_components=None))]
pub(super) fn parse_rdfile<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: Option<usize>,
    max_records: Option<usize>,
    max_fields_per_record: Option<usize>,
    max_line_len: Option<usize>,
    max_reaction_components: Option<usize>,
) -> PyResult<Bound<'py, PyAny>> {
    let mut limits = chematic_mol::RdfileParseLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_records {
        limits.max_records = v;
    }
    if let Some(v) = max_fields_per_record {
        limits.max_fields_per_record = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    if let Some(v) = max_reaction_components {
        limits.max_reaction_components = v;
    }
    let document = chematic_mol::parse_rdfile_with_limits(text, &limits)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let json =
        serde_json::to_string(&document).map_err(|e| PyValueError::new_err(e.to_string()))?;
    py.import("json")?.call_method1("loads", (json,))
}

/// Write current typed RDFile records; source_text is provenance, never replayed.
#[pyfunction]
#[pyo3(signature=(document, *, wrap_data=true))]
pub(super) fn write_rdfile(document: &Bound<'_, PyDict>, wrap_data: bool) -> PyResult<String> {
    let json: String = document
        .py()
        .import("json")?
        .call_method1("dumps", (document,))?
        .extract()?;
    let document: chematic_mol::RdfileDocument =
        serde_json::from_str(&json).map_err(|e| PyValueError::new_err(e.to_string()))?;
    chematic_mol::rdfile::write_rdfile_with_options(
        &document,
        chematic_mol::rdfile::RdfileWriteOptions { wrap_data },
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Explicit lossy graph projection in record order. Registry-only records yield
/// None; fields and coordinates remain available through parse_rdfile instead.
#[pyfunction]
pub(super) fn rdfile_to_smiles(text: &str) -> PyResult<Vec<Option<String>>> {
    let document =
        chematic_mol::parse_rdfile(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    document
        .records
        .iter()
        .map(|r| r.smiles().map_err(|e| PyValueError::new_err(e.to_string())))
        .collect()
}
