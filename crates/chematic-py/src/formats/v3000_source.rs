use super::*;

/// Validate one V3000 MOL block and retain its exact UTF-8 source.
/// Coordinates, IDs and opaque fields survive without editing or interpreting
/// their chemistry/stereo. Multiple CTABs and malformed core structures fail.
/// Accepts str or UTF-8 bytes, using native source-retention limits.
#[pyfunction]
pub(super) fn roundtrip_mol_v3000_source(
    #[pyo3(from_py_with = serialized_text)] text: &str,
) -> PyResult<String> {
    chematic_mol::roundtrip_mol_v3000_source(text)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
