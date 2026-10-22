use super::*;

/// Validate every local-profile SDF record and retain the exact UTF-8 document.
/// Any failed record/limit returns no output; opaque metadata/stereo is retained
/// without interpretation or editing. Default bounds: 16 MiB input/record,
/// 1 MiB physical line, 10000 records, 100000 atoms per record and in total.
#[pyfunction]
pub(super) fn roundtrip_sdf_source(
    #[pyo3(from_py_with = serialized_text)] text: &str,
) -> PyResult<String> {
    chematic_mol::roundtrip_sdf_source(text)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

/// Extract one unchanged MOL CTAB, omitting SD fields/delimiter. Additional
/// records fail instead of selecting the first. Uses SDF source bounds with
/// a one-record cap; the native core profile is not full chemistry validation.
#[pyfunction]
pub(super) fn single_sdf_mol_block(
    #[pyo3(from_py_with = serialized_text)] text: &str,
) -> PyResult<String> {
    chematic_mol::single_sdf_mol_block(text)
        .map(str::to_owned)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
