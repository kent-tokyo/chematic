use super::*;

/// Project one checked HIN record to a molecular graph. Source charges/types
/// and coordinates are available separately via parse_hin_record_json.
#[pyfunction]
pub(super) fn from_hin(text: &str) -> PyResult<Mol> {
    let record = chematic_mol::parse_hin_record(text)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(Mol {
        inner: Arc::new(record.to_molecule()),
        props: Default::default(),
    })
}

/// Inspect HIN atom identities, coordinates, partial charges and literal bonds.
#[pyfunction]
pub(super) fn parse_hin_record_json(text: &str) -> PyResult<String> {
    chematic_mol::parse_hin_record(text)
        .map(|record| record.to_json().to_string())
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Validate the bounded HIN profile and return exact original UTF-8 source.
#[pyfunction]
pub(super) fn roundtrip_hin_record(text: &str) -> PyResult<String> {
    chematic_mol::parse_hin_record(text)
        .map(|record| chematic_mol::write_hin_record(&record))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
