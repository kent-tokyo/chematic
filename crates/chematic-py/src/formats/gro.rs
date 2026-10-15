use super::*;

/// Read one GRO frame without guessing elements from ambiguous atom names.
/// Coordinates/box remain in nm and velocities in nm/ps. Optional canonical
/// element symbols must follow source row order; complete original text is kept.
#[pyfunction]
#[pyo3(signature=(text, *, elements=None, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_gro<'py>(
    py: Python<'py>,
    text: &str,
    elements: Option<Vec<String>>,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut limits = chematic_mol::GroParseLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_atoms {
        limits.max_atoms = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    let mut frame = chematic_mol::parse_gro_with_limits(text, &limits)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let policy = if let Some(elements) = elements {
        let elements = elements
            .iter()
            .map(|s| {
                chematic_core::Element::from_symbol(s)
                    .filter(|e| e.symbol() == s)
                    .ok_or_else(|| {
                        PyValueError::new_err(format!("invalid canonical element symbol '{s}'"))
                    })
            })
            .collect::<PyResult<Vec<_>>>()?;
        frame = frame
            .with_elements(&elements)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        "caller_supplied"
    } else {
        "unresolved"
    };
    let atoms = frame
        .atoms
        .iter()
        .map(|atom| {
            let d = PyDict::new(py);
            d.set_item("residue_number", atom.residue_number)?;
            d.set_item("residue_name", &atom.residue_name)?;
            d.set_item("atom_name", &atom.atom_name)?;
            d.set_item("atom_number", atom.atom_number)?;
            d.set_item("position_nm", atom.position_nm)?;
            d.set_item("velocity_nm_per_ps", atom.velocity_nm_per_ps)?;
            d.set_item("element", atom.element.map(|e| e.symbol()))?;
            Ok(d)
        })
        .collect::<PyResult<Vec<_>>>()?;
    let result = PyDict::new(py);
    result.set_item("title", frame.title)?;
    result.set_item("atoms", atoms)?;
    result.set_item("element_policy", policy)?;
    result.set_item("box_values_nm", frame.box_values_nm)?;
    result.set_item("box_vectors_nm", frame.box_vectors_nm)?;
    result.set_item("source_text", frame.source_text)?;
    Ok(result)
}
