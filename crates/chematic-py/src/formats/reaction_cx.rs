use super::*;

/// Normalize enhanced-atom-stereo reaction CXSMILES without moving group
/// membership to another atom (RDKit #9732). Supports a/&n/on groups only;
/// other CX metadata, invalid references and non-tetrahedral members raise
/// ValueError. Components sort within each role when canonical=True.
#[pyfunction]
#[pyo3(signature = (text, *, canonical = true))]
pub(super) fn normalize_reaction_cxsmiles(
    #[pyo3(from_py_with = serialized_text)] text: &str,
    canonical: bool,
) -> PyResult<String> {
    let reaction = chematic_rxn::parse_reaction_cxsmiles(text)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    chematic_rxn::write_reaction_cxsmiles(&reaction, canonical)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
