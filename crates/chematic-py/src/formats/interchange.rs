use super::*;

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

/// Checked MOL2/CJSON conversion with a shared machine-readable report.
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
