use super::*;

/// Read a standalone nonperiodic or 3D TURBOMOLE coordinate document.
/// Returns explicit elements, Cartesian angstrom coordinates, frozen flags,
/// row lattice vectors and original text. Unsupported data groups fail.
#[pyfunction]
#[pyo3(signature=(text, *, max_input_bytes=None, max_atoms=None, max_line_len=None))]
pub(super) fn parse_turbomole<'py>(
    py: Python<'py>,
    text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_line_len: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut limits = chematic_mol::TurbomoleParseLimits::default();
    if let Some(v) = max_input_bytes {
        limits.max_input_bytes = v;
    }
    if let Some(v) = max_atoms {
        limits.max_atoms = v;
    }
    if let Some(v) = max_line_len {
        limits.max_line_len = v;
    }
    let geometry = chematic_mol::parse_turbomole_with_limits(text, &limits)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let result = PyDict::new(py);
    result.set_item("periodicity", geometry.periodicity())?;
    result.set_item(
        "elements",
        geometry
            .atoms
            .iter()
            .map(|a| a.element.symbol())
            .collect::<Vec<_>>(),
    )?;
    result.set_item(
        "positions_angstrom",
        geometry
            .atoms
            .iter()
            .map(|a| a.position_angstrom)
            .collect::<Vec<_>>(),
    )?;
    result.set_item(
        "frozen",
        geometry.atoms.iter().map(|a| a.frozen).collect::<Vec<_>>(),
    )?;
    result.set_item("lattice_angstrom", geometry.lattice_angstrom)?;
    result.set_item("source_text", geometry.source_text)?;
    Ok(result)
}

/// Write explicit geometry in canonical bohr units, including its 3D lattice.
/// No topology or atom perception; arrays must have exactly matching lengths.
#[pyfunction]
#[pyo3(signature=(elements, positions_angstrom, *, lattice_angstrom=None, frozen=None))]
pub(super) fn write_turbomole(
    elements: Vec<String>,
    positions_angstrom: Vec<[f64; 3]>,
    lattice_angstrom: Option<[[f64; 3]; 3]>,
    frozen: Option<Vec<bool>>,
) -> PyResult<String> {
    if elements.len() != positions_angstrom.len()
        || frozen.as_ref().is_some_and(|v| v.len() != elements.len())
    {
        return Err(PyValueError::new_err(
            "element, position and frozen-flag counts must match",
        ));
    }
    let mut atoms = Vec::new();
    for (i, (s, position_angstrom)) in elements.into_iter().zip(positions_angstrom).enumerate() {
        let element = chematic_core::Element::from_symbol(&s)
            .filter(|e| e.symbol() == s)
            .ok_or_else(|| {
                PyValueError::new_err(format!("invalid canonical element symbol '{s}'"))
            })?;
        atoms.push(chematic_mol::TurbomoleAtom {
            element,
            position_angstrom,
            frozen: frozen.as_ref().is_some_and(|v| v[i]),
        });
    }
    let geometry = chematic_mol::TurbomoleGeometry {
        atoms,
        lattice_angstrom,
        source_text: String::new(),
    };
    chematic_mol::write_turbomole(&geometry).map_err(|e| PyValueError::new_err(e.to_string()))
}
