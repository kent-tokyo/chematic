use super::*;

/// Project MOL2 to CJSON with checked graph/xyz/name/declared partial charges.
/// Incomplete or undeclared charge columns fail; no values are computed/filled.
/// Other MOL2 metadata and stereo are not retained or verified. Use
/// convert_with_report for a complete loss report. Native conversion limits apply.
#[pyfunction]
pub(super) fn project_mol2_to_cjson(
    #[pyo3(from_py_with = serialized_text)] text: &str,
) -> PyResult<String> {
    chematic_mol::project_mol2_to_cjson(text)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Replace only MOL2 ar bond tokens with explicit single/double orders.
/// Bounded neutral C/N/O/S profile; other source bytes stay unchanged.
#[pyfunction]
pub(super) fn kekulize_mol2_record(text: &str) -> PyResult<String> {
    chematic_mol::kekulize_mol2_record(text)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Inspect a source-preserving MOL2 record, including partial charges and types.
#[pyfunction]
pub(super) fn parse_mol2_record_json(text: &str) -> PyResult<String> {
    chematic_mol::parse_mol2_source_record(text)
        .map(|record| record.to_json().to_string())
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Validate one MOL2 record and return its exact source, preserving metadata.
#[pyfunction]
pub(super) fn roundtrip_mol2_record(text: &str) -> PyResult<String> {
    chematic_mol::parse_mol2_source_record(text)
        .map(|record| chematic_mol::write_mol2_source_record(&record))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Inspect a strictly validated immutable CJSON record.
#[pyfunction]
pub(super) fn parse_cjson_record_json(text: &str) -> PyResult<String> {
    chematic_mol::parse_cjson_record(text)
        .map(|record| record.to_json().to_string())
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Validate the local CJSON profile and return its exact UTF-8 source.
#[pyfunction]
pub(super) fn roundtrip_cjson_record(text: &str) -> PyResult<String> {
    chematic_mol::parse_cjson_record(text)
        .map(|record| chematic_mol::write_cjson_record(&record))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Checked MOL2/CJSON/CML conversion with a shared machine-readable report.
#[pyfunction]
#[pyo3(signature = (text, input_format, output_format, policy="strict"))]
pub(super) fn convert_with_report_json(
    text: &str,
    input_format: &str,
    output_format: &str,
    policy: &str,
) -> PyResult<String> {
    chematic_mol::convert_with_report_json(
        text,
        input_format,
        output_format,
        policy,
        &Default::default(),
    )
    .map_err(|error| PyValueError::new_err(error.to_string()))
}
